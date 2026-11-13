"""
plate_maps.py
=============

Utilities for building, inspecting and filing compound plate maps for the
SMDC high-throughput screening workflow, backed by the HiTS MySQL database.

WHAT THIS MODULE IS FOR
-----------------------
Screens are run on 96/384/1536-well plates. Two kinds of information have to be
reconciled before results can be loaded into HiTS:

  1. What the scientist *intended*, recorded by hand in an Excel workbook of
     plate maps (one sheet per destination plate, laid out as grids).
  2. What the inventory database *knows*: which compound lot physically sits in
     which well of which source plate.

The functions here parse (1), look up (2), join them, and write out the
per-plate manifests plus the folder structure that the HiTS loader expects.

THE THREE ENTRY POINTS
----------------------
  draw_plate_heatmap(...)
      Visual QC. Renders one plate as an annotated heatmap. Reads from the
      database if you give it a barcode, or plots a DataFrame you already have.

  make_platemaps_from_barcode(...)
      Given a physical plate barcode, pull its full contents from HiTS and
      write an Excel workbook: one sheet per attribute (alias, SMDC ID, common
      name, lot) laid out in plate geometry, plus a long-form sheet.

  make_dr_platemaps_and_folders(...)
      The main pipeline. Consumes a dose-response workbook, resolves every
      destination well to (SMDC_ID, Lot, source plate, source well,
      concentration, assay ID), writes a `*_DR_plate.xlsx` manifest per plate,
      builds the `<assay>/toload_<i>/<date>/` folder tree, and copies the raw
      instrument files in beside it.

      clean_up_folders(...) is its counterpart: run it after HiTS ingest to
      move `toload_*` to `loaded/`.

INPUT WORKBOOK CONTRACT
-----------------------
`make_dr_platemaps_and_folders` expects an .xlsx where:

  * A sheet named `Info` (or `Sheet1`) is the driver table, one row per
    destination plate, with columns Grid_num, Source_plate, Raw_assay_ID,
    Lot_info, Date_of_screen, Raw_data_file, Assay_name. Column names are
    matched case/punctuation-insensitively; see _standardize_info_columns.

  * For every value of Grid_num, a sheet of that name (Excel truncates sheet
    names to 31 chars, hence the `[0:31]` slicing you will see) holding
    stacked grid BLOCKS. See read_platemap_block for the block convention:
    block 0 = SMDC IDs, block 1 = concentrations (uM), block 2 = a variable
    third field whose meaning is set by the Source_plate / Lot_info columns.

  * The Source_plate cell is a small dispatch language, not just a barcode.
    See get_source_plate_info for the accepted values.

ENVIRONMENT ASSUMPTIONS (these will bite a new user first)
----------------------------------------------------------
  * `mysql_cnx.py` sits next to this file (run notebooks from this folder, or
    put this folder on sys.path).
  * HiTS is reachable through mysql_cnx (by default localhost, i.e. an SSH
    tunnel). Lookups read from QB3_HiTS unless `hits_db` is set.
  * Raw instrument data is on the Montreal share, found by default_share_root():
    `/Volumes/Shared/` on macOS, `Z:` (`/mnt/z/` in WSL) on Windows. Pass
    `share_root=` or set `smdc_share_root` for any other mount.
  * It is meant to be driven from a notebook: importing it configures logging,
    and make_dr_platemaps_and_folders calls the notebook builtin `display()`.
    Nothing prompts on stdin, so the pipeline itself runs unattended.
"""

import pandas as pd

# Opt in to pandas 3.0 behaviour now so .replace()/.fillna() downcasting is
# explicit rather than silently changing dtypes.
pd.set_option('future.no_silent_downcasting', True)

import os
import numpy as np
import hashlib
import time
import humanize
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import shutil
import filecmp
import re
import datetime
from sqlalchemy import bindparam, text
from sqlalchemy.exc import DBAPIError

import logging

# Logging is configured at import time. That is a side effect on whoever imports
# this module, which is fine for a notebook-driven tool and would not be in a
# package.
#
# basicConfig() is a no-op when handlers already exist, and Jupyter installs one,
# so the explicit setLevel() below is what actually makes this module's INFO
# output visible. That is the job the old basicConfig(force=True) was doing, by
# tearing out Jupyter's handler; setting the level instead leaves other people's
# logging alone.
#
# The default is WARNING, so new operators only see what needs their attention.
# The pipeline logs its routine decisions at INFO: call set_log_level('INFO')
# (the notebook's LOG_LEVEL setting) to see them while debugging.
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger()  # root logger
logger.setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)  # root is noisy; mute matplotlib



def set_log_level(level='WARNING'):
    """Show pipeline messages at this level and above: 'WARNING' (default,
    only what needs attention) or 'INFO' (every step, for debugging)."""
    logging.getLogger().setLevel(getattr(logging, str(level).upper()) if isinstance(level, str) else level)


sns.set_theme(rc={'figure.dpi': 150, 'figure.figsize': (5, 3.75)})

from mysql_cnx import create_mysql_cnx

# Plate maps only READ from HiTS, so they use production unless the hits_db
# environment variable says otherwise (mysql_cnx's own default is the nightly
# test database, which is right for writes but stale for plate lookups).
READ_DB = os.environ.get('hits_db') or 'QB3_HiTS'

def _is_wsl():
    """True inside Windows Subsystem for Linux."""
    if os.environ.get('WSL_DISTRO_NAME'):
        return True
    try:
        with open('/proc/version') as f:
            return 'microsoft' in f.read().lower()
    except OSError:
        return False


SHARE_ROOT_HELP = ("Set SHARE_ROOT in the notebook's settings cell (or pass share_root=, or set the "
                   "smdc_share_root environment variable) to the folder where the Montreal share is "
                   "mounted, e.g. '/mnt/z/' in WSL, '/Volumes/Shared/' on a Mac, 'Z:/' on Windows.")


def default_share_root():
    """Where the Montreal share (raw instrument data) is mounted on this machine.

    Order: the `smdc_share_root` environment variable if set; else the first
    of the standard mounts that exists: `/Volumes/Shared/` (macOS) or
    `/mnt/z/` (WSL, i.e. Windows drive Z:); else `Z:/` on native Windows.
    The expected setup is WSL on Windows with Z: mounted at /mnt/z. Anywhere
    else, or when /mnt/z isn't mounted, this raises and asks for the mount
    point instead of guessing. Pass `share_root=` to use a different drive.
    """
    env = os.environ.get('smdc_share_root')
    if env:
        return env
    for candidate in ('/Volumes/Shared/', '/mnt/z/'):
        if os.path.isdir(candidate):
            return candidate
    if os.name == 'nt':
        return 'Z:/'
    if _is_wsl():
        raise FileNotFoundError(
            "This is WSL, but the Montreal share is not mounted at /mnt/z. Mount it with "
            "`sudo mkdir -p /mnt/z && sudo mount -t drvfs Z: /mnt/z -o uid=$(id -u),gid=$(id -g)` "
            "(see the WSL setup doc for the /etc/fstab line), or, if it is mounted somewhere else: "
            + SHARE_ROOT_HELP)
    raise FileNotFoundError("Can't find the Montreal share. " + SHARE_ROOT_HELP)


def _connect():
    """Open a read connection to HiTS (see READ_DB) and check it works.

    create_engine is lazy, so a dead SSH tunnel would otherwise only show up
    at the first query, or worse, inside a lookup whose error is handled.
    Fail here instead, with a message that names the likely cause.
    """
    engine = create_mysql_cnx(db=READ_DB)
    try:
        with engine.connect() as conn:
            conn.execute(text("select 1"))
    except DBAPIError as e:
        engine.dispose()
        url = engine.url
        raise ConnectionError(
            f"Cannot reach HiTS ({url.database} on {url.host}:{url.port}). "
            "Is the SSH tunnel up and are the mysql_* environment variables right? "
            f"Original error: {e.orig}"
        ) from e
    return engine

# Excel-style plate row labels, in order: A..Z then AA..AF.
# Slices of this list give the rows for each plate format:
#   96 -> [0:8] (A-H), 384 -> [0:16] (A-P), 1536 -> [0:32] (A-AF).
# Module-level and mutable, so treat it as read-only.
# NOTE: `read_platemap_block` rebinds `plate_rows` as a local int, so this
# global is shadowed inside that function. Confusing, but not a live bug.
plate_rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
              'V', 'W', 'X', 'Y', 'Z', 'AA', 'AB', 'AC', 'AD', 'AE', 'AF']


###############################################################################
# SECTION 1: small generic helpers
###############################################################################

def _stringify_plate_label(value, missing='-'):
    """Convert a plate annotation value to display text."""
    if pd.isna(value):
        return missing
    text = str(value)
    return missing if text == 'nan' else text


def _join_plate_labels(values, sep='\n', missing='-'):
    """Join grouped plate annotations while tolerating mixed dtypes."""
    return sep.join(_stringify_plate_label(value, missing=missing) for value in values)


###############################################################################
def now_str():
    """ A helper function to generate a timestamp string for file naming.
    Args:
        None
    Returns:
        str: A date string in the format YYYYMMDD.
    """
    return time.strftime("%Y%m%d")


###############################################################################
def get_md5_hash(path, size):
    """ Generate the md5 checksum of a file, reading it in chunks.

    Used to verify that a raw data file copied intact. Currently only called
    from commented-out verification code in make_dr_platemaps_and_folders.

    Args:
        path (str): path to file
        size (int): chunk size in bytes to read at a time (65536 is typical);
            chunking keeps large instrument files off the heap
    Returns:
        str: hex digest of the md5 sum.
    """
    m = hashlib.md5()
    with open(path, 'rb') as f:
        b = f.read(size)
        while len(b) > 0:
            m.update(b)
            b = f.read(size)
    return m.hexdigest()


###############################################################################
# SECTION 2: plate geometry (well naming and index -> well assignment)
###############################################################################
def plate_well_list(size=384):
    """ Generate a list of well names for a plate of given size.
    Args:
        size (int): The size of the plate. Default is 384.
    Returns:
        list: Well names in row-major order, e.g. ['A01', 'A02', ... 'P24'].
            Columns are zero-padded to 2 digits, which is the convention HiTS
            and the liquid handlers both use.
    """
    plate_rows_dict = {384: plate_rows[0:16], 1536: plate_rows[0:32], 96: plate_rows[0:8]}
    plate_cols_dict = {384: np.arange(1, 25), 1536: np.arange(1, 49), 96: np.arange(1, 13)}
    rows = plate_rows_dict[size]
    cols = plate_cols_dict[size]

    return [row + str(col).zfill(2) for row in rows for col in cols]


###############################################################################
def assign_wells_by_idx(df, plate_col='Destination Plate', well_col='Destination Well', idx_col='idx', size=384,
                        exclude_cols=['01', '02', '23', '24']):
    """ Assign wells to a dataframe based on the index of the dataframe.
    Args:
        df (DataFrame): The dataframe to assign wells to.
        plate (str): The name of the column containing the plate name(s)
        wells (str): The name of the column to assign wells to.
        idx (str): The name of the column to use as the index. Index must start at 0.
        size (int): The size of the plate. Default is 384.
        exclude_cols (list): Column labels (as zero-padded strings) to keep free,
            typically the outer columns reserved for controls. Default leaves
            the four edge columns of a 384w plate empty.
    Returns:
        DataFrame: The dataframe with wells assigned. Row order follows the
            per-plate groups, so reset_index downstream if order matters.

    The index restarts at 0 for each plate in plate_col, so idx_col must be a
    within-plate index, not a global one.
    """
    df[idx_col] = df[idx_col].astype(int)

    plate_rows_dict = {384: plate_rows[0:16], 1536: plate_rows[0:32], 96: plate_rows[0:8]}
    rows = plate_rows_dict[size]
    plate_cols_dict = {384: np.arange(1, 25), 1536: np.arange(1, 49), 96: np.arange(1, 13)}
    cols = plate_cols_dict[size]

    all_wells = [row + str(col).zfill(2) for row in rows for col in cols]
    assert len(all_wells) == size, f"Error: {len(all_wells)} wells were created, but {size} were expected."

    # Drop the reserved columns. NOTE: this is a substring test against the
    # whole well name, so an exclude value like '01' would also knock out
    # nothing else here, but a single-digit value like '1' would match 'A01',
    # 'A10'..'A19' and 'A21'. Keep exclude_cols zero-padded and 2 chars.
    all_wells = [well for well in all_wells if all(col not in well for col in exclude_cols)]
    # After exclusion len(all_wells) < size, so the highest indices map to NaN.
    # That is intentional: an over-long idx range fails loudly as empty wells.
    wellidx = dict(zip(np.arange(0, size), all_wells))

    dfs = []
    for plate in df[plate_col].unique():
        tmp = df[df[plate_col] == plate]
        # NOTE: chained assignment on a slice. Works because of the .copy()
        # semantics pandas currently applies, but will warn//break under
        # copy-on-write. tmp = tmp.copy() first if you touch this.
        tmp[well_col] = tmp[idx_col].map(wellidx)
        dfs.append(tmp)
    df = pd.concat(dfs)
    return df


###############################################################################
def draw_plate_heatmap(barcode='OTPARP1234-Eco1', df=None, wells='Well', values='SMDC_ID', labels='SMDC_ID',
                       plate_size=None, center=None, vmin=None, vmax=None, width=25, outdir=None, engine=None,
                       cbar=None, cmap=None):
    """
    Draw a heatmap of one plate, for visual QC of a layout or a result.

    Args:
        barcode (str): Plate barcode. Used as the title and output filename,
            and to query HiTS if df is None.
        df (DataFrame, optional): Long-form data, one row per well. If None,
            it is fetched via make_platemaps_from_barcode(barcode).
        wells (str): Name of the column holding well names like 'A01'.
        values (str): Column to colour by. Must be numeric-coercible.
        labels (str): Column to print inside each cell, or None for no labels.
        plate_size (int, optional): 96, 384 or 1536. Inferred from the observed
            rows/cols if None.
        center, vmin, vmax (float, optional): Passed through to sns.heatmap.
        width (int): Figure width in inches; height is width/3.
        outdir (path, optional): If given, save <barcode>.png here and close the
            figure instead of leaving it open.
        cbar (bool, optional): show the colourbar. Default: on, except for SMDC IDs.
        cmap (optional): passed through to sns.heatmap.
    Returns:
        tuple: (df, pl, plabels) - the cleaned long-form frame, the pivoted
            value grid, and the pivoted label grid. Returned mainly for
            debugging what got plotted.
    Side effects:
        Draws (and optionally writes) a plate heatmap.

    Special case: when labels == 'SMDC_ID' the labels are taken from the value
    grid itself and the colourbar is suppressed, since colouring by ID is
    meaningless; you just want to read the IDs off the plate.
    """
    if df is None:
        if engine is not None:
            df = make_platemaps_from_barcode(barcode=barcode, outdir=None, cnx=engine)
        else:
            df = make_platemaps_from_barcode(barcode=barcode, outdir=None)
    df = df[~df[wells].isna()].copy()
    df[wells] = df[wells].astype(str)
    # Split 'A01' / 'AA01' into row letters and column number. Columns are
    # always the last 2 chars by convention, so the remainder is the row.
    df['Col'] = df[wells].str[-2:]
    # NOTE: str.strip(y) strips any of the *characters* in y, not the suffix.
    # It works today only because row letters and digits do not overlap.
    # x.removesuffix(y) would say what is meant.
    df['Row'] = [x.strip(y) for x, y in zip(df[wells], df.Col)]
    df.Col = df.Col.astype(int)

    # 1536 -> 32 rows (A-AF), matching every other size dict in this file.
    plate_rows_dict = {96: plate_rows[0:8], 384: plate_rows[0:16], 1536: plate_rows[0:32]}
    plate_cols_dict = {96: np.arange(1, 13), 384: np.arange(1, 25), 1536: np.arange(1, 49)}
    # plate_size = plate_size or (384 if len(df) > 96 else 96)
    # Infer the plate format from the smallest geometry that contains every
    # observed row and column. A sparse 384w plate using only A-H/1-12 will be
    # misread as 96w; pass plate_size explicitly for partial plates.
    if plate_size is None:
        if set(df.Row.unique()).issubset(set(plate_rows_dict[96])) and set(df.Col.unique()).issubset(
                set(plate_cols_dict[96])):
            plate_size = 96
        elif set(df.Row.unique()).issubset(set(plate_rows_dict[384])) and set(df.Col.unique()).issubset(
                set(plate_cols_dict[384])):
            plate_size = 384
        else:
            plate_size = 1536

    # Categorical + observed=False in the pivot below is what forces empty wells
    # to still appear in the grid, so the heatmap has true plate geometry
    # rather than only the wells that happen to have data.
    df.Row = pd.Categorical(df.Row, categories=plate_rows_dict[plate_size], ordered=True)
    df.Col = pd.Categorical(df.Col, categories=plate_cols_dict[plate_size], ordered=True)

    cbar_override, cbar = cbar, True
    pl = df.pivot_table(columns='Col', index='Row', values=values, dropna=False, observed=False)
    pl = pl.astype(float)
    if labels is not None:
        plabels = df.pivot_table(columns='Col', index='Row', values=labels, aggfunc=_join_plate_labels, dropna=False,
                                 observed=False)
        # -99 is used throughout this module as a "no value" sentinel that
        # survives an int cast, then is blanked out again on the way to display.
        plabels = plabels.fillna(-99).astype(str).replace('-99', '', regex=True)
        # for col in plabels.columns:
        # plabels[col]=plabels[col].str.replace('.0', '', regex=False)
        # For IDs, print the numeric grid itself and drop the colourbar.
        if labels == 'SMDC_ID':
            plabels = pl.copy()
            plabels = plabels.fillna(-99).astype(int).astype(str).replace('-99', '', regex=True)
            cbar = False
    else:
        plabels = None
    if cbar_override is not None:
        cbar = cbar_override

    # Always reserve the colourbar column, hidden when there is no colourbar,
    # so every plate drawn at the same width lines up well for well (e.g. the
    # layout and concentration maps of one DR plate, drawn one above the other).
    fig, (ax, cax) = plt.subplots(1, 2, figsize=(width, width / 3), gridspec_kw={"width_ratios": [1, .01]})
    if not cbar:
        cax.set_visible(False)
    sns.heatmap(pl, annot=plabels, fmt='s', linewidth=1, linecolor='white', ax=ax, cbar_ax=cax if cbar else None, cbar=cbar,
                center=center, vmin=vmin, vmax=vmax, cmap=cmap,
                annot_kws={"size": 7}
                )
    ax.grid(which='major')
    ax.set_title(f"{barcode} {values}");
    # Fixed margins rather than tight_layout, which would size the plate
    # differently with and without colourbar tick labels.
    fig.subplots_adjust(left=0.04, right=0.95, bottom=0.1, top=0.92, wspace=0.02)
    if outdir is not None:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        plt.savefig(os.path.join(outdir, f'{barcode}.png'), bbox_inches='tight')
        plt.close()
    return df, pl, plabels


###############################################################################
# SECTION 3: HiTS database reads
#
# Every query here is a SQLAlchemy text() with bound parameters (:name, and
# bindparam(..., expanding=True) for IN lists). Never build SQL by putting
# values into the string with f-strings, % or +, even for barcodes or IDs
# from the lab's own files: a quote in a value breaks the query, and this
# code is public. Keep new queries parameterised the same way.
###############################################################################
def make_platemaps_from_barcode(barcode='OTPARP1234-Eco1', outdir='./', cnx=None):
    """
    Create an excel file with each sheet showing a plate map with vendor alias,
    smdc_id, and lot for downstream reference.

    Args:
        barcode (str): Unique plate barcode to pull from database
        outdir (path): Place to save excel file. If None, no file will be saved
            and the long-form frame is simply returned. The DR pipeline calls
            it this way.
        cnx (Engine, optional): An open SQLAlchemy connection to reuse. If None,
            one is opened and disposed here. Pass one in when looping over many
            barcodes.
    Returns:
        df: Long form table of all relevant plate information, one row per well
            (a well may appear more than once if it has multiple aliases; see
            the filtering step below).

    Schema note: INV_WELL is the physical well, INV_COMP_LOT the compound lot in
    it, INV_COMPOUND the compound, INV_ALIAS its various names. The joins are
    LEFT so empty wells still come back.

    TODO: decompose further. This does four separable things in one body (query,
    collapse the alias fan-out, add plate geometry, write the workbook), which is
    why it cannot be unit tested without a database. Split it, but only once,
    and delete the monolith when you do.
    """
    # pull data from hits
    close = False
    if cnx is None:
        close = True
        cnx = _connect()

    query = text("""
    select IW_COORDINATES, ICL_ALIAS as Vendor_Alias, ICL_IC_ID as SMDC_ID, IA_VALUE as Common_Name, IA_IAT_ID as Alias_Type,
        ICL_LOT_NUM as Lot_Num, IW_CONC_UM as Conc_uM, IW_FINAL_CONC_NM as Conc_nM, IC_SMILES as SMILES
    from INV_WELL left join INV_COMP_LOT on IW_ICL_ID = ICL_ID left join INV_ALIAS on IA_IC_ID=ICL_IC_ID left join INV_COMPOUND on IC_COMPOUND_ID=ICL_IC_ID
    where IW_IP_ID = 
        (select IP_ID from INV_PLATE where IP_BARCODE = :barcode);
    """)
    df = pd.read_sql(query, cnx, params={'barcode': barcode})
    if close:
        cnx.dispose()
    if df.empty:
        raise ValueError(f"Source plate {barcode!r} has no wells in HiTS. Check the barcode in the "
                         f"Info tab's Source_plate column. {SOURCE_PLATE_HELP}")

    # The INV_ALIAS join fans each well out to one row per alias type. Collapse
    # back to one row per well: prefer the Common_Name alias, else keep an
    # arbitrary row and null the name out.
    wells = []
    for well in df.IW_COORDINATES.unique():
        tmp = df[df.IW_COORDINATES == well]
        tmp = tmp[tmp.Alias_Type == 'Common_Name']  # this will be one row
        if len(tmp) == 0:
            tmp = pd.DataFrame(df[df.IW_COORDINATES == well].iloc[0]).T  # this will be one row
            tmp['Common_Name'] = np.nan
        wells.append(tmp)
    df = pd.concat(wells)
    try:
        df = df.drop(columns='Alias_Type')
    except KeyError:
        pass

    # format data
    # Row is everything before the 2-digit column, so 1536-well rows AA-AF stay
    # separate from row A.
    df['Row'] = df.IW_COORDINATES.str.slice(stop=-2)
    df['Col'] = df.IW_COORDINATES.str.slice(start=-2).astype(int)
    df = df.rename(columns={'IW_COORDINATES': 'Well'})
    df['Plate'] = barcode
    cols = ['Vendor_Alias', 'SMDC_ID', 'Common_Name', 'Lot_Num']

    # return or save
    if outdir is None:
        return df

    # save as plate layouts: one sheet per attribute, each pivoted into plate
    # geometry so it can be read like the physical plate. Text columns are
    # joined, numeric ones averaged (a no-op when each well is unique).
    # The 'None' / -99 sentinels exist only to survive the pivot's aggfunc,
    # and are converted back to NaN immediately afterwards.
    else:
        with pd.ExcelWriter(os.path.join(outdir, f'{barcode}_plate_maps.xlsx')) as writer:
            for col in cols:
                if col in ['Vendor_Alias', 'Common_Name']:
                    aggfunc = lambda x: _join_plate_labels(x, sep=' ', missing='None')
                else:
                    df[col] = df[col].fillna(-99).astype(int)
                    aggfunc = 'mean'
                tmpdf = df.pivot_table(index='Row', columns='Col', values=col, aggfunc=aggfunc)
                # plate order A..Z then AA..AF, not the lexical A, AA, AB, ..., B
                tmpdf = tmpdf.loc[sorted(tmpdf.index, key=lambda r: (len(r), r))]
                tmpdf = tmpdf.replace('None', np.nan).replace(-99, np.nan)
                tmpdf.index.name = col
                tmpdf.to_excel(writer, sheet_name=col)
            df = df.replace('None', np.nan).replace(-99, np.nan)
            df.to_excel(writer, sheet_name='Long_Form')
        return df


###############################################################################
# (lookup_icl_id_by_plate_and_ic_id, further down, is the variant that returns the
# ICL_ID itself, for registering plates.)
# The four lookup_* functions below are the fallback ladder for resolving a
# compound to a lot and a physical source well. get_lot_info tries them in
# roughly this order of confidence:
#   lookup_id_and_lot_by_plate_and_well  - most specific: you know plate + well
#   lookup_lot_by_plate_and_ic_id        - you know the compound and its plate
#   lookup_plate_by_lot_num              - you know the compound and its lot
#   lookup_lot_by_user                   - last resort: newest lot this person made
###############################################################################
def lookup_lot_by_user(ids, users):
    """ Find the most recent (highest) lotnum that the user has created for their compound.

    Last-resort fallback: use only when nothing else identifies the lot, since
    "newest lot registered by this person" is an assumption, not a fact about
    what was pipetted.

    Args:
        ids (list of strings): The list of SMDC IDs to look up. Integers, but in string form.
        users (list of string tuples): The list of users in the form ("FirstName","LastName").
            Required: there is no default person.
    Returns:
        fullresults (DataFrame): the SMDC IDs and Lots as two columns in a df.
            One block of rows per user, concatenated; if two users both made
            lots of the same compound you get duplicate SMDC_IDs.
    """
    if not users:
        raise ValueError(
            "Some compounds have no lot, and the last-resort lookup (newest lot registered by a person) "
            "needs to know who registered them. Pass users=[(\"First\", \"Last\")] (USERS in the notebook), "
            "or give those wells a lot or source plate in the plate map.")
    # pull data from hits
    cnx = _connect()
    res = []
    for fname, lname in users:
        logger.info('looking up users..')
        query = text("select DM_ID from DIRECTORY_MEMBER where DM_FNAME = :fname and DM_LNAME = :lname")
        user = pd.read_sql(query, cnx, params={'fname': fname, 'lname': lname})
        if user.empty:
            raise ValueError(f"No HiTS user named {fname} {lname} in DIRECTORY_MEMBER")
        user = user.DM_ID.iloc[0]

        logger.info('looking up compound lot by user..')
        query = text("""
        select IC_COMPOUND_ID as SMDC_ID, max(ICL_LOT_NUM) as Lot from INV_COMPOUND 
        inner join INV_COMP_LOT on IC_COMPOUND_ID=ICL_IC_ID
        where IC_COMPOUND_ID in :ids and ICL_INS_BY = :user
        group by IC_COMPOUND_ID
        """).bindparams(bindparam('ids', expanding=True))
        fullresults = pd.read_sql(query, cnx, params={'ids': list(ids), 'user': int(user)})
        res.append(fullresults)
    cnx.dispose()
    # create df
    logger.info('creating df..')
    fullresults = pd.concat(res).reset_index(drop=True)
    return fullresults


###############################################################################
def lookup_plate_by_lot_num(ids=['285699', '130958'], lot_nums=['1', '1'], like_filter=['Source plate', "96-MA"],
                            replace_vals=[]):
    """ Find which source plates and wells hold a given (compound, lot) pair.

    Docstring corrected: this does not find lot numbers, it finds plates. ids
    and lot_nums are consumed pairwise, so they must be equal length and
    aligned.

    Args:
        ids (list of strings): The list of SMDC IDs to look up. Integers, but in string form.
        lot_nums (list of strings): Lots, aligned elementwise with ids. The
            literal string 'nan' is skipped, which is how unresolved lots pass
            through silently.
        like_filter (list of strings): [column, substring] pair, applied as a
            substring filter on the result, e.g. keep only mother-assay plates.
            None to skip.
        replace_vals (list): SMDC IDs to skip entirely (controls, DMSO, etc.);
            see _build_replace_metadata.
    Returns:
        fullresults (DataFrame): SMDC_ID, Lot, Source plate, Source well.
            Raises from pd.concat if every id was skipped and dfs is empty.
    """
    # pull data from hits
    cnx = _connect()
    # metadata = MetaData()

    res = []

    logger.info('looking up plates by lot_num..')
    dfs = []
    for (smdc, lot_num) in zip(ids, lot_nums):
        # logger.info(f"Looking up plate for {smdc}, lot {lot_num}")
        smdc = str(smdc)
        lot_num = str(lot_num)
        if smdc not in replace_vals:
            if lot_num != 'nan':
                query = text("""
                select ICL_IC_ID as SMDC_ID, ICL_LOT_NUM as Lot, IP_BARCODE as 'Source plate', IW_COORDINATES as 'Source well'
                from INV_PLATE inner join INV_WELL on INV_PLATE.IP_ID = INV_WELL.IW_IP_ID 
                inner join INV_COMP_LOT on INV_WELL.IW_ICL_ID = INV_COMP_LOT.ICL_ID 
                where ICL_IC_ID = :smdc and ICL_LOT_NUM = :lot_num
                """)
                fullresults = pd.read_sql(query, cnx, params={'smdc': smdc, 'lot_num': lot_num})
                dfs.append(fullresults)
    cnx.dispose()

    # create df
    logger.info('creating df..')
    fullresults = pd.concat(dfs).reset_index(drop=True)
    if like_filter is not None:
        fullresults = fullresults[fullresults[like_filter[0]].str.contains(like_filter[1])]
    return fullresults


###############################################################################
# def lookup_plates_by_ico_id(id=135):
#     cnx = _connect()
#     query=f"""
#     select IP_ID from INV_PLATE where IP_ICO_ID in

def lookup_lot_by_plate_and_ic_id(ids=[], plates=[], replace_vals={}, matched=False):
    """ Find the lot and source well for compounds known to sit on given plates.

    Docstring corrected: this returns whatever lots are on the plate, not the
    highest lot a user made.

    Args:
        ids (list of ints): The list of SMDC IDs to look up.
        plates (list of strings): The list of plates corresponding to each molecule.
        replace_vals (dict): Values (controls/solvents) to skip; membership is
            tested with `in`, so a list works too.
        matched (bool): False = cross product, every id against every plate
            (len(ids) * len(plates) queries; slow and can return a compound
            found on several plates). True = zip ids with plates elementwise.
            get_lot_info always passes True.
    Returns:
        fullresults (DataFrame): SMDC_ID, Lot, Source plate, Source well.
    """
    # pull data from hits
    cnx = _connect()
    # metadata = MetaData()

    res = []

    dfs = []
    if not matched:
        logger.info('looking up lots by plate and ic_id..')
        for smdc in ids:
            for plate in plates:
                if smdc not in replace_vals:
                    query = text("""
                    select ICL_IC_ID as SMDC_ID, ICL_LOT_NUM as Lot, IP_BARCODE as 'Source plate', IW_COORDINATES as 'Source well'
                    from INV_PLATE inner join INV_WELL on INV_PLATE.IP_ID = INV_WELL.IW_IP_ID 
                    inner join INV_COMP_LOT on INV_WELL.IW_ICL_ID = INV_COMP_LOT.ICL_ID 
                    where ICL_IC_ID = :smdc and IP_BARCODE = :plate
                    """)
                    fullresults = pd.read_sql(query, cnx, params={'smdc': smdc, 'plate': plate})
                    dfs.append(fullresults)
    else:
        logger.info('looking up lots by matched plate and ic_id..')
        for smdc, plate in zip(ids, plates):
            if smdc not in replace_vals:
                query = text("""
                select ICL_IC_ID as SMDC_ID, ICL_LOT_NUM as Lot, IP_BARCODE as 'Source plate', IW_COORDINATES as 'Source well'
                from INV_PLATE inner join INV_WELL on INV_PLATE.IP_ID = INV_WELL.IW_IP_ID 
                inner join INV_COMP_LOT on INV_WELL.IW_ICL_ID = INV_COMP_LOT.ICL_ID 
                where ICL_IC_ID = :smdc and IP_BARCODE = :plate
                """)
                fullresults = pd.read_sql(query, cnx, params={'smdc': smdc, 'plate': plate})
                dfs.append(fullresults)
    cnx.dispose()

    # create df
    logger.info('creating df..')
    # Concat FIRST, then normalise case on the combined frame. These two lines
    # used to be the other way round, which applied .upper() to whatever frame
    # the last loop iteration happened to leave bound and then threw it away,
    # so barcodes were never actually normalised and case-sensitive merges
    # downstream silently dropped rows.
    if len(dfs) == 0:
        logger.warning("No source plates matched; returning an empty frame.")
        return pd.DataFrame(columns=['SMDC_ID', 'Lot', 'Source plate', 'Source well'])
    fullresults = pd.concat(dfs).reset_index(drop=True)
    fullresults['Source plate'] = fullresults['Source plate'].str.upper()
    return fullresults


###############################################################################
def _clean_lookup_value(value, numeric=False):
    """Normalise an id / lot / plate / well for matching; None means "not given".

    With numeric=True, 1191231, 1191231.0 and '1191231.0' all become '1191231'
    (pandas turns int columns into floats as soon as one NaN is present). Leave
    numeric=False for plates and wells: a barcode like '0123' must keep its zero.
    """
    if value is None:
        return None
    if isinstance(value, str):
        s = value.strip()
        if s == '' or s.lower() == 'nan':
            return None
    else:
        if pd.isna(value):
            return None
        s = str(value)
    if numeric:
        try:
            f = float(s)
        except ValueError:
            return s
        return str(int(f)) if f.is_integer() else s
    return s


def _upper_or_none(value):
    return value.upper() if value is not None else None


def lookup_icl_id_by_plate_and_ic_id(ids, plates, lots=None, wells=None, replace_vals=(), strict=True, cnx=None):
    """ Look up the inventory lot id (INV_COMP_LOT.ICL_ID, a.k.a. LOT_ID) for compounds on given plates.

    Sibling of lookup_lot_by_plate_and_ic_id, but returns the ICL_ID that
    INV_WELL.IW_ICL_ID needs when registering new plates, and is stricter: a
    compound that is not on its plate, or whose plate holds more than one lot of
    it, is an error rather than a silently missing or duplicated row.

    ids/plates (and lots/wells, if given) are zipped elementwise. One query is
    made per distinct plate, and all queries are parameterised.

    Args:
        ids (sequence): SMDC IDs (INV_COMP_LOT.ICL_IC_ID). Ints, floats or strings.
        plates (sequence of str): Source plate barcodes, aligned with ids. Case
            handling in the query follows the DB collation (case-insensitive on
            MySQL's default); results are matched back to requests ignoring case.
        lots (sequence, optional): ICL_LOT_NUM aligned with ids. Narrows the match;
            use it when a plate can hold several lots of one compound. Missing/NaN
            entries mean "any lot".
        wells (sequence of str, optional): Source wells aligned with ids. Narrows
            the match to one physical well; missing/NaN entries mean "any well".
        replace_vals (iterable): SMDC IDs to skip (controls, DMSO, ...).
        strict (bool): True raises ValueError if any requested compound matches no
            well or more than one distinct ICL_ID. False logs a warning and returns
            whatever matched (every candidate row for an ambiguous request).
        cnx: SQLAlchemy engine. If None, one is created and disposed here.
    Returns:
        DataFrame: SMDC_ID, Lot, Source plate (upper-cased, as in the sibling
        lookups), Source well, LOT_ID (= ICL_ID). One row per matching well.

    Example, attaching LOT_ID to a dose-response transfer file:
        res = lookup_icl_id_by_plate_and_ic_id(
            df['SMDC_ID'], df['Source plate'], lots=df['Lot'], wells=df['Source well'])
        df['_plate'] = df['Source plate'].str.upper()
        df = df.merge(res, left_on=['SMDC_ID', 'Lot', '_plate', 'Source well'],
                      right_on=['SMDC_ID', 'Lot', 'Source plate', 'Source well'], how='left',
                      suffixes=('', '_db')).drop(columns=['_plate', 'Source plate_db'])
    """
    n = len(ids)
    for name, seq in (('plates', plates), ('lots', lots), ('wells', wells)):
        if seq is not None and len(seq) != n:
            raise ValueError(f"{name} has {len(seq)} entries but ids has {n}; they are zipped elementwise")
    ids, plates = list(ids), list(plates)
    lots = list(lots) if lots is not None else [None] * n
    wells = list(wells) if wells is not None else [None] * n

    skip = {_clean_lookup_value(v, numeric=True) for v in replace_vals}
    req = pd.DataFrame({
        'smdc': [_clean_lookup_value(v, numeric=True) for v in ids],
        'plate_raw': [_clean_lookup_value(v) for v in plates],
        'lot': [_clean_lookup_value(v, numeric=True) for v in lots],
        'well': [_upper_or_none(_clean_lookup_value(v)) for v in wells],
    }, dtype=object)
    req['plate'] = req['plate_raw'].map(_upper_or_none)
    req = req[req['smdc'].notna() & ~req['smdc'].isin(skip)].drop_duplicates().reset_index(drop=True)

    cols = ['SMDC_ID', 'Lot', 'Source plate', 'Source well', 'LOT_ID']
    sql = text("""
        select ICL_IC_ID as SMDC_ID, ICL_LOT_NUM as Lot, IP_BARCODE as 'Source plate',
               IW_COORDINATES as 'Source well', ICL_ID as LOT_ID
        from INV_PLATE inner join INV_WELL on INV_PLATE.IP_ID = INV_WELL.IW_IP_ID
        inner join INV_COMP_LOT on INV_WELL.IW_ICL_ID = INV_COMP_LOT.ICL_ID
        where IP_BARCODE in :plates and ICL_IC_ID in :ids
    """).bindparams(bindparam('plates', expanding=True), bindparam('ids', expanding=True))

    own_cnx = cnx is None
    if own_cnx:
        cnx = _connect()
    frames = []
    try:
        logger.info('looking up ICL_IDs by plate and ic_id..')
        for _, grp in req.dropna(subset=['plate']).groupby('plate'):
            params = {'plates': list(grp['plate_raw'].unique()),
                      'ids': [int(i) if i.isdigit() else i for i in grp['smdc'].unique()]}
            frames.append(pd.read_sql(sql, cnx, params=params))
    finally:
        if own_cnx:
            cnx.dispose()

    found = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)
    found = found[cols].copy()
    found['Source plate'] = found['Source plate'].str.upper()
    key = pd.DataFrame({
        '_fid': found.index,
        'smdc': found['SMDC_ID'].map(lambda v: _clean_lookup_value(v, numeric=True)),
        'plate': found['Source plate'],
        '_lot': found['Lot'].map(lambda v: _clean_lookup_value(v, numeric=True)),
        '_well': found['Source well'].map(lambda v: _upper_or_none(_clean_lookup_value(v))),
    })
    m = req.reset_index(names='_req').merge(key, on=['smdc', 'plate'], how='left')
    keep = m['_fid'].notna() & (m['lot'].isna() | (m['lot'] == m['_lot'])) & (m['well'].isna() | (m['well'] == m['_well']))
    matched = m[keep].copy()
    matched['_icl'] = found.loc[matched['_fid'].astype(int), 'LOT_ID'].to_numpy()

    n_icl = matched.groupby('_req')['_icl'].nunique().reindex(req.index, fill_value=0)
    problems = []
    for label, bad in (('no matching well', n_icl == 0), ('more than one lot (pass lots=/wells= to disambiguate)', n_icl > 1)):
        if bad.any():
            sample = req.loc[bad, ['smdc', 'plate_raw', 'lot', 'well']].head(10).to_dict('records')
            problems.append(f"{int(bad.sum())} request(s) with {label}, e.g. {sample}")
    if problems:
        msg = 'lookup_icl_id_by_plate_and_ic_id: ' + '; '.join(problems)
        if strict:
            raise ValueError(msg)
        logger.warning(msg)

    fids = sorted(matched['_fid'].astype(int).unique())
    return found.loc[fids, cols].reset_index(drop=True)


###############################################################################
def lookup_id_and_lot_by_plate_and_well(plates=[], wells=[]):
    """ Look up the compound ID and lot physically present at given plate/well pairs.

    The most trustworthy of the four lookups: plate + well identifies exactly
    one physical well. Used for Biomek-style long-format worklists, where the
    source well is already stated.

    Args:
        plates (list of strings): Plate barcodes.
        wells (list of strings): Wells, aligned elementwise with plates.
    Returns:
        fullresults (DataFrame): SMDC_ID, Lot, Source plate, Source well.
            One query per pair, so this is O(n) round trips.
    """
    cnx = _connect()
    dfs = []
    for plate, well in zip(plates, wells):
        query = text("""
        select ICL_IC_ID as SMDC_ID, ICL_LOT_NUM as Lot, IP_BARCODE as 'Source plate', IW_COORDINATES as 'Source well'
        from INV_PLATE inner join INV_WELL on INV_PLATE.IP_ID = INV_WELL.IW_IP_ID 
        inner join INV_COMP_LOT on INV_WELL.IW_ICL_ID = INV_COMP_LOT.ICL_ID 
        where IP_BARCODE = :plate and IW_COORDINATES = :well
        """)
        fullresults = pd.read_sql(query, cnx, params={'plate': plate, 'well': well})
        dfs.append(fullresults)
    cnx.dispose()

    # create df
    logger.info('creating df..')
    fullresults = pd.concat(dfs).reset_index(drop=True)
    fullresults['Source plate'] = fullresults['Source plate'].str.upper()
    return fullresults


###############################################################################
# SECTION 4: parsing the hand-written Excel plate maps
###############################################################################
def _row_labels(n):
    """ Plate row labels: A..P for 16, A..Z then AA..AF for 32. """
    letters = [chr(ord('A') + i) for i in range(26)]
    return (letters + ['A' + x for x in letters])[:n]


def _infer_platemap_size(plate_map):
    """ 96, 384 or 1536 from the row labels in a plate map sheet's first column.

    The largest label present decides: AF only exists on 1536, P on 384 and
    up, H on 96 and up. Falls back to 384, which then fails with
    MissingBlockError if the sheet has no complete A-P block.
    """
    labels = set(plate_map.iloc[:, 0].dropna().astype(str).str.strip().str.upper())
    if 'AF' in labels:
        return 1536
    if 'P' in labels:
        return 384
    if 'H' in labels:
        return 96
    return 384


class MissingBlockError(ValueError):
    """The plate map sheet has fewer complete blocks than the one asked for."""


def read_platemap_block(plate_map=None, block=0, value_name='SMDC_ID', size=None):
    """ Read one grid block out of a block-formatted plate map sheet.

    THE BLOCK CONVENTION (this is the thing to understand about the workbook):
    a single sheet stacks several plate-shaped grids vertically. Each block has
    a header row of column numbers (1, 2, ... N) and then one row per plate
    row, with the row label (A, B, ...) in the first column. Completely blank
    rows are dropped first, so blank lines between blocks don't matter. Blocks
    are found by their A..P (A..AF for 1536) label runs, so a sheet without
    the per-block header rows also reads correctly.

    By convention in this pipeline:
        block 0 -> SMDC IDs
        block 1 -> concentrations in uM (optional: may be missing, or hold a
                   plate suffix on cherry-pick maps)
        block 2 -> a context-dependent third field (Lot, Source plate, or
                   Assay ID), decided by the Info sheet's Source_plate and
                   Lot_info columns. See get_source_plate_info / get_lot_info.

    Args:
        plate_map (DataFrame): sheet read with header=None-ish semantics; the
            first column holds row letters, the rest hold values.
        block (int): which block, 0 at the top.
        size (int, optional): 96, 384 or 1536. Default None: infer it from the
            row labels (AF means 1536, else P means 384, else H means 96).
    Returns:
        df (DataFrame): long-format block data with columns row, col, value_name.
            col is a zero-padded 2-char string so `row + col` reconstructs the
            well name.

    NOTE: `plate_rows` is rebound here as an int, shadowing the module-level
    list of row labels.
    """
    if size is None:
        size = _infer_platemap_size(plate_map)
    if size == 384:
        plate_cols = 24
        plate_rows = 16
    elif size == 96:
        plate_cols = 12
        plate_rows = 8
    elif size == 1536:
        plate_cols = 48
        plate_rows = 32
    else:
        raise ValueError(f"size must be 96, 384 or 1536; got {size}")
    cols = ['row']
    cols.extend(range(1, plate_cols + 1))
    # Find each block by its row labels (A, B, ... in the first column) rather
    # than by a fixed stride. The real maps have a header row per block, which
    # the old stride of plate_rows + 1 read correctly. A sheet that left out
    # those header rows would have been shifted by one row per block under the
    # stride; label matching reads both.
    sheet = plate_map.dropna(how='all', axis='index').reset_index(drop=True)
    labels = sheet.iloc[:, 0].astype(str).str.strip().str.upper().tolist()
    expected = _row_labels(plate_rows)
    starts = [i for i in range(len(labels) - plate_rows + 1)
              if labels[i:i + plate_rows] == expected]
    if block >= len(starts):
        raise MissingBlockError(
            f"Plate map has {len(starts)} complete block(s) with rows A to {expected[-1]} in the "
            f"first column, so block {block} does not exist. Check the row labels and that no "
            "populated row was added inside a block."
        )
    df = sheet.iloc[starts[block]:starts[block] + plate_rows, 0:plate_cols + 1].copy()
    df.columns = cols
    df['row'] = expected
    df = df.melt(id_vars=['row'], value_vars=range(1, plate_cols + 1), ignore_index=True,
                 value_name=value_name, var_name='col')
    # Zero-pad so row + col == 'A01', matching HiTS well coordinates.
    df.col = df.col.astype(str).str.rjust(2, '0')
    return df


###############################################################################
# SECTION 5: column-name normalisation
#
# These exist because the Info sheet is hand-typed. "Grid num", "grid_num" and
# "GridNum" should all mean Grid_num. Everything is compared on a normalised
# key (lowercase, alphanumerics only) and renamed to one canonical spelling, so
# the rest of the pipeline can assume exact column names.
###############################################################################
def _norm_key(s):
    """Normalise a column name to lowercase alphanumerics for fuzzy matching.

    'Grid num' / 'grid_num' / 'GridNum' all collapse to 'gridnum'.
    """
    return re.sub(r'[^a-z0-9]+', '', str(s).strip().lower())


def _safe_lower(v):
    """Lowercase a possibly-NaN cell, mapping NaN to ''.

    Lets the dispatch logic use `'map' in _safe_lower(x)` without guarding for
    NaN at every call site.
    """
    if pd.isna(v):
        return ''
    return str(v).strip().lower()


def _series_get_ci(row, key, default=np.nan):
    """Fetch a value from a Series by case/punctuation-insensitive key.

    Args:
        row (Series): typically one row of the Info sheet.
        key (str): canonical column name, e.g. 'Grid_num'.
        default: returned when no column normalises to key.
    Returns:
        The cell value, or default.
    """
    key_norm = _norm_key(key)
    idx_map = {_norm_key(k): k for k in row.index}
    if key_norm in idx_map:
        return row[idx_map[key_norm]]
    return default


def _canonicalize_columns(df, canonical_map):
    """Rename columns to canonical spellings, matching on normalised keys.

    Args:
        df (DataFrame): frame whose columns may be spelled loosely.
        canonical_map (dict): {normalised_key: canonical_name}.
    Returns:
        DataFrame: renamed copy. Columns not in the map are left alone.
    """
    norm_to_current = {_norm_key(c): c for c in df.columns}
    rename = {}
    for norm_key, canonical in canonical_map.items():
        if norm_key in norm_to_current and norm_to_current[norm_key] != canonical:
            rename[norm_to_current[norm_key]] = canonical
    return df.rename(columns=rename)


def _standardize_info_columns(info):
    """Canonicalise the Info sheet's columns. This is the Info sheet's schema."""
    return _canonicalize_columns(info, {
        'gridnum': 'Grid_num',
        'sourceplate': 'Source_plate',
        'rawassayid': 'Raw_assay_ID',
        'lotinfo': 'Lot_info',
        'dateofscreen': 'Date_of_screen',
        'rawdatafile': 'Raw_data_file',
        'assayname': 'Assay_name',
    })


def _standardize_dids_columns(dids):
    """Canonicalise the destination-plate frame's columns.

    Called defensively after nearly every merge in the DR pipeline, because
    each lookup contributes its own spelling of the concentration column.
    Note several distinct source spellings all fold into '[compound] uM'.
    """
    return _canonicalize_columns(dids, {
        'smdcid': 'SMDC_ID',
        'lot': 'Lot',
        'sourceplate': 'Source plate',
        'sourcewell': 'Source well',
        'destinationplate': 'Destination plate',
        'destinationwell': 'Destination well',
        'assayid': 'Assay ID',
        'conc': '[compound] uM',
        'compoundum': '[compound] uM',
        'destinationplateconcentration': '[compound] uM',
    })


def _build_replace_metadata(replace_dict=None):
    """Build the map of non-compound well contents, and the list of values to skip.

    Wells labelled DMSO / Water / Background are not compounds, so their
    SMDC_ID is nulled out and they are skipped by the database lookups.
    Caller-supplied entries in replace_dict are merged first, so the built-in
    solvent/control entries take precedence over any conflicting override.

    Args:
        replace_dict (dict, optional): extra {label: replacement} pairs, e.g.
            mapping a project's control name to np.nan or to a real SMDC ID.
    Returns:
        tuple: (replace_dict, replace_vals) where replace_vals is the
            deduplicated set of replacement values, with NaN rendered as the
            string 'nan' so it can be compared against stringified IDs.
    """
    replace_dict = dict(replace_dict or {})
    replace_dict.update({'DMSO': np.nan, 'dmso': np.nan})
    replace_dict.update({'Water': np.nan, 'water': np.nan})
    replace_dict.update({'Background': np.nan, 'background': np.nan})
    replace_dict.update({'nan': np.nan, u'\xa0': np.nan, '': np.nan})

    replace_vals = list(set(replace_dict.values()))
    replace_vals = [x if str(x) != 'nan' else 'nan' for x in replace_vals]
    logger.info(f"SMDC ID replace values are {replace_vals}")
    return replace_dict, replace_vals


###############################################################################
# SECTION 5b: shared plumbing for the folder pipeline
#
# make_dr_platemaps_and_folders and clean_up_folders MUST agree about where a
# folder lives, or clean-up silently fails to find what the pipeline built.
# These four helpers are the single source of truth for that agreement. Change
# them here, not at the call sites.
###############################################################################
def _normalize_share_path(raw):
    r"""Normalise a hand-typed Raw_data_file value to a path relative to the share root.

    The value is typed in by hand from whichever machine acquired the data, so it
    may arrive as a Windows drive letter (Z:\...), a UNC path, an smb:// URL, a
    macOS /Volumes mount, or a raw IP. Each branch strips one known prefix; the
    magic numbers are len(prefix).

    Fragile and append-only by nature: expect to add a branch the first time
    someone uses a new machine. That is exactly why this lives in one place now.

    Args:
        raw (str): the Raw_data_file cell.
    Returns:
        str: path relative to the share root, with no leading slash, suitable
            for os.path.join(share_root, result).

    Every branch slices exactly len(prefix). The 'Montreal/Shared/' branch used
    to slice [17:] against a 16-character prefix, eating the first character of
    the remaining path ('SMDC/run.txt' -> 'MDC/run.txt'), so that branch could
    never resolve to a real file. Fixed; test_plate_maps.py covers it.
    """
    p = str(raw)
    p = (p.replace('\\\\', '/').replace('\\', '/')
         .replace('Z:', '').replace('Y:', '').replace('A:', '')
         )
    p = re.sub(r'^/*\d{1,3}(?:\.\d{1,3}){3}/shared/', '', p, flags=re.IGNORECASE)
    if p.startswith('/shared/'):
        p = p[8:]
    elif p.startswith('Montreal/Shared/'):
        p = p[16:]  # len('Montreal/Shared/') == 16
    elif p.startswith('shared/'):
        p = p[7:]
    elif p.startswith('/Volumes/Shared/'):
        p = p[16:]
    elif re.match(r'^/mnt/[A-Za-z]/', p):
        p = p[7:]  # WSL view of a Windows drive, e.g. /mnt/z/ for Z:
    elif p.lower().startswith('smb://montreal.ucsf.edu/shared/'):
        p = p[31:]
    elif p.startswith('smb://montreal/Shared/'):
        p = p[22:]
    elif p.startswith('mb://montreal.ucsf.edu/shared/'):
        p = p[30:]
    elif p.lower().startswith("/montreal.ucsf.edu/shared/"):
        p = p[26:]
    elif re.match(r'^/*montreal/shared/', p, flags=re.IGNORECASE):
        p = re.sub(r'^/*montreal/shared/', '', p, flags=re.IGNORECASE)  # UNC \\montreal\Shared\
    elif p.startswith('/'):
        # Must stay last: a leading slash would make os.path.join ignore the
        # share root entirely and return an absolute path.
        p = p[1:]
    return p


_SHARE_PREFIX = re.compile(
    r'^(?:[AYZayz]:[\\/]|\\\\|//|/Volumes/Shared/|/mnt/[A-Za-z]/|smb://|/shared/|shared/|Montreal/Shared/)',
    re.IGNORECASE)


def resolve_share_path(path, share_root=None):
    r"""Turn a share path typed on any machine into one that works on this one.

    A path that already exists is returned unchanged. Otherwise, if it looks
    like a Montreal share path (Z:\..., \\server\..., /Volumes/Shared/...,
    /mnt/z/..., smb://...), it is re-rooted on this machine's share mount, so
    a Z:\ path copied from File Explorer works in WSL (/mnt/z/...) and a Mac
    path works on Windows. Anything else (relative paths, other folders) is
    returned unchanged.
    """
    if path is None:
        return path
    p = str(path)
    if os.path.exists(p) or not _SHARE_PREFIX.match(p):
        return p
    return os.path.join(share_root or default_share_root(), _normalize_share_path(p))


def share_path(path, share_root=None, what='path'):
    """resolve_share_path, then check the result exists.

    For notebook settings: a wrong path or an unmounted share fails here,
    with both the typed and the mapped path, instead of several cells later.
    Args:
        path (str): a path as typed or copied (Z:\\..., /Volumes/Shared/..., /mnt/z/...).
        share_root (str, optional): where the share is mounted. Default: default_share_root().
        what (str): the setting's name, for the message (e.g. 'INFILE').
    Returns:
        str: the path on this machine.
    """
    resolved = resolve_share_path(path, share_root)
    if not os.path.exists(resolved):
        mapped = f" (looked for {resolved!r})" if resolved != str(path) else ""
        raise FileNotFoundError(f"{what} not found: {str(path)!r}{mapped}. Check the path for typos. "
                                f"If the share is mounted somewhere else: {SHARE_ROOT_HELP}")
    return resolved



class _FileInUse:
    """with _FileInUse(path): ... turns a PermissionError while writing path
    into one that says to close the file in Excel (or wherever it is open)."""
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        return self

    def __exit__(self, kind, err, tb):
        if kind is not None and issubclass(kind, PermissionError):
            raise PermissionError(
                f"Can't write {self.path}: it is open in another program (Excel, a text editor, Pipeline Pilot) "
                "or you don't have write permission there. Close it everywhere else and run this cell again."
            ) from err
        return False


def _copy_file(src, dst):
    """Copy a file's contents, then its timestamps if the filesystem allows it.

    shutil.copy2 fails with PermissionError on some network mounts (a drvfs
    mount of Z: in WSL owned by root, for example): the data copies fine but
    setting the modification time is refused. Keep the copy and log a
    warning instead of stopping the run.
    """
    with _FileInUse(dst):
        shutil.copyfile(src, dst)
    try:
        shutil.copystat(src, dst)
    except OSError as e:
        logger.warning(f"Copied {dst} but could not keep the original timestamps ({e}). "
                       "In WSL, mount Z: with your uid/gid to avoid this; see the setup notes.")
    return dst


def _copy_tree(src, dst, ignore=None):
    """copytree(..., dirs_exist_ok=True) that uses _copy_file for each file.

    shutil.copytree also sets directory timestamps and raises shutil.Error at
    the end if that is refused, which a WSL network mount does. Directory
    timestamps don't matter here, so they aren't copied.
    """
    def fail(err):  # an unreadable folder must stop the copy, not be skipped
        raise err

    for root, dirs, files in os.walk(src, onerror=fail):
        skip = set(ignore(root, dirs + files)) if ignore else set()
        dirs[:] = [d for d in dirs if d not in skip]
        target = os.path.join(dst, os.path.relpath(root, src))
        os.makedirs(target, exist_ok=True)
        for name in files:
            if name not in skip:
                _copy_file(os.path.join(root, name), os.path.join(target, name))


def _missing_from_copy(src, dst, ignore=None):
    """Files under src (minus ignored ones) that are absent from dst or differ
    in size there. Empty means dst holds a full copy."""
    missing = []
    for root, dirs, files in os.walk(src, onerror=lambda err: missing.append(str(err))):
        skip = set(ignore(root, dirs + files)) if ignore else set()
        dirs[:] = [d for d in dirs if d not in skip]
        for name in files:
            if name in skip:
                continue
            a = os.path.join(root, name)
            b = os.path.join(dst, os.path.relpath(a, src))
            if not os.path.isfile(b) or os.path.getsize(a) != os.path.getsize(b):
                missing.append(a)
    return missing


def _read_info_sheet(infile):
    """Read and normalise the workbook's Info sheet.

    Single source of truth for how the Info sheet is parsed, so
    make_dr_platemaps_and_folders and clean_up_folders cannot drift apart in
    which rows they keep or how they parse dates. They previously did both.

    Args:
        infile (path): .xlsx (all sheets) or .csv (Info sheet only).
    Returns:
        tuple: (dp, info) where dp is {sheet_name: DataFrame} for the whole
            workbook and info is the normalised Info sheet. Date_of_screen is an
            int like 20260717.
    """
    if 'xlsx' in infile:
        dp = pd.read_excel(infile, engine='openpyxl', sheet_name=None)
    else:
        dp = {}
        try:
            dp['Info'] = pd.read_csv(infile)
        except Exception:
            dp['Sheet1'] = pd.read_csv(infile)
    info_sheet = 'Info' if 'Info' in dp else 'Sheet1'
    info = dp[info_sheet]
    # pandas reads cells like 'NA' and 'n/a' as empty. In Source_plate those
    # mean "no source plate", so keep that column as typed (only truly empty
    # cells become NaN). Same shape and index as `info`.
    if 'xlsx' in infile:
        raw = pd.read_excel(infile, engine='openpyxl', sheet_name=info_sheet,
                            keep_default_na=False, na_values=[''])
    else:
        raw = pd.read_csv(infile, keep_default_na=False, na_values=[''])

    info = _standardize_info_columns(info)
    # Some workbooks carry a title row above the real header.
    if "Grid_num" not in info.columns:
        info.columns = info.iloc[0]
        info = info.drop(index=0)
        info = _standardize_info_columns(info)
    if 'Source_plate' in info.columns and raw.shape == dp[info_sheet].shape:
        pos = list(info.columns).index('Source_plate')
        info['Source_plate'] = raw.iloc[:, pos].reindex(info.index).values

    info = info.dropna(how='all', axis='index').reset_index(drop=True)
    info = info.dropna(subset='Date_of_screen', how='all', axis='index').reset_index(drop=True)
    try:
        info.Date_of_screen = pd.to_datetime(info.Date_of_screen.astype(int).astype(str), format='mixed')
    except Exception:
        info.Date_of_screen = pd.to_datetime(info.Date_of_screen)
    info.Date_of_screen = info.Date_of_screen.astype(str).str.replace('-', '').astype(int)
    info = info.dropna(how='all', axis='columns')
    return dp, info


SOURCE_PLATE_HELP = (
    'Fill in Source_plate for every plate in the Info tab. Accepted values (type them exactly, '
    'case does not matter):\n'
    '  "use lots"   lots are drawn in block 2 (bare "lots" or "lot" also works)\n'
    '  "use map"    source plate barcodes are drawn in block 2 (bare "map" also works)\n'
    '  "use biomek" lots come from a Biomek worklist tab (bare "biomek" also works)\n'
    '  "use my name" (or "use users", "NA", "n/a", "n.a.", "n.a", "none")   no source plate; '
    'each compound gets the newest lot registered by USERS; source plate and well are written as NA\n'
    '  one or more source plate barcodes separated by ";", e.g. "PLATE1-MA01; PLATE2-MA01"')


NO_SOURCE_PLATE_WORDS = {'na', 'n/a', 'n.a', 'n.a.', 'none',
                         # same rule, said as "look the lot up by my name"
                         'use my name', 'use users', 'use user', 'use name', 'my name', 'users'}


OLD_SCREEN_DAYS = 92   # about 3 months

_ENVISION_STARTED = re.compile(r'^\s*Assay Started:\s*,+\s*([^,]+)', re.IGNORECASE)


def envision_run_date(path):
    """The date an EnVision export says the assay ran ('Assay Started:'), or
    None if the file isn't an EnVision CSV or has no readable date."""
    if not str(path).lower().endswith('.csv') or not os.path.isfile(path):
        return None
    try:
        with open(path, newline='', errors='replace') as f:
            for line in f:
                m = _ENVISION_STARTED.match(line)
                if m:
                    value = m.group(1).strip()
                    for fmt in ('%m/%d/%Y %I:%M:%S %p', '%m/%d/%Y %H:%M:%S', '%m/%d/%Y'):
                        try:
                            return pd.Timestamp(datetime.datetime.strptime(value, fmt)).normalize()
                        except ValueError:
                            pass
                    stamp = pd.to_datetime(value, errors='coerce')
                    return None if pd.isna(stamp) else stamp.normalize()
    except OSError:
        return None
    return None


def check_screen_dates(info, max_age_days=OLD_SCREEN_DAYS, today=None, share_root=None):
    """Warn about Date_of_screen values that look wrong.

    For each Info row whose Raw_data_file is an EnVision CSV with an
    'Assay Started:' date, Date_of_screen is compared with that date and any
    difference is flagged. Rows without one (other instruments, file not
    found) fall back to the age rule: more than ~3 months old, or in the
    future, is flagged and the operator is asked whether it is historical
    data or a wrong date. Warnings only, so historical data still goes
    through.

    Args:
        info (DataFrame): the parsed Info sheet (Date_of_screen as YYYYMMDD ints).
        max_age_days (int): the fallback age limit. Default about 3 months.
        today (date, optional): for tests. Default: today.
        share_root (str, optional): where the share is mounted, to find
            Raw_data_file. Default: default_share_root(); if the share
            can't be found, every row uses the age rule.
    Returns:
        DataFrame: the flagged rows with a 'reason' column (empty if none).
    """
    cols = ['Raw_assay_ID', 'Assay_name', 'Date_of_screen', 'Grid_num', 'reason']
    if 'Date_of_screen' not in info.columns:
        return pd.DataFrame(columns=cols)
    today = pd.Timestamp(today or pd.Timestamp.today().normalize())
    dates = pd.to_datetime(info['Date_of_screen'].astype(str), format='%Y%m%d', errors='coerce')

    def run_date(raw):
        if not isinstance(raw, str) or not raw.strip():
            return None
        try:
            root = share_root or default_share_root()
            path = raw if os.path.exists(raw) else os.path.join(root, _normalize_share_path(raw))
        except (FileNotFoundError, ValueError):
            path = raw
        return envision_run_date(path)

    raw_col = info['Raw_data_file'] if 'Raw_data_file' in info.columns else pd.Series(None, index=info.index)
    reasons = []
    for date, raw in zip(dates, raw_col):
        ran = run_date(raw)
        if pd.isna(date):
            reasons.append('unreadable Date_of_screen')
        elif ran is not None:
            reasons.append('' if ran == date else f'the EnVision file says it ran {ran:%Y-%m-%d}')
        elif (today - date).days < 0:
            reasons.append('in the future')
        elif (today - date).days > max_age_days:
            reasons.append(f'{(today - date).days} days old')
        else:
            reasons.append('')
    flagged = info.assign(_date=dates, reason=reasons)[[r != '' for r in reasons]]
    if flagged.empty:
        return pd.DataFrame(columns=cols)
    keys = [c for c in ('Raw_assay_ID', 'Assay_name', '_date', 'reason') if c in flagged.columns]
    plates = flagged['Grid_num'].astype(str) if 'Grid_num' in flagged.columns else pd.Series('?', index=flagged.index)
    versions = flagged.assign(_plates=plates).groupby(keys, dropna=False)['_plates'].agg(', '.join).reset_index()
    lines = [f"  assay {v.get('Raw_assay_ID', '?')} {v.get('Assay_name', '')}: Date_of_screen "
             f"{v._date:%Y-%m-%d} - {v.reason}; plates {v._plates}" if not pd.isna(v._date) else
             f"  assay {v.get('Raw_assay_ID', '?')} {v.get('Assay_name', '')}: {v.reason}; plates {v._plates}"
             for _, v in versions.iterrows()]
    logger.warning(
        f"Check Date_of_screen for {len(versions)} assay version(s):\n" + '\n'.join(lines) +
        f"\nWhere the EnVision file has a different date, the Info tab is probably wrong. Otherwise, dates over "
        f"{max_age_days} days old are fine for historical data; if this isn't historical, fix the date "
        "(often the year) and rerun step 2.")
    have = [c for c in cols if c in flagged.columns]
    return flagged[have].reset_index(drop=True)


def check_destination_plates(info, cnx=None):
    """Warn when a destination plate (Grid_num) is listed twice in the Info
    tab, or already exists in HiTS inventory as an IP_BARCODE or IP_GRID_NUM.

    Uploading creates the destination plate in inventory, so an existing one
    usually means the barcode was reused or an earlier upload of it got part
    way. It's a warning, not an error: a deliberate rerun is fine.

    Returns:
        dict: {'duplicated': [grid_nums listed more than once],
               'in_hits': DataFrame of matching INV_PLATE rows}.
    """
    grids = [str(g).strip() for g in info.get('Grid_num', pd.Series(dtype=str)).dropna()]
    grids = [g for part in grids for g in part.split(';') if g.strip()]
    grids = [g.strip() for g in grids]
    dup = sorted({g for g in grids if grids.count(g) > 1})
    if dup:
        logger.warning(f"Destination plate(s) listed more than once in the Info tab: {', '.join(dup)}. "
                       "Each Grid_num must be a different plate.")
    in_hits = pd.DataFrame(columns=['IP_BARCODE', 'IP_GRID_NUM', 'IP_INS_DATE'])
    if grids:
        close = cnx is None
        cnx = cnx or _connect()
        try:
            query = text("select IP_BARCODE, IP_GRID_NUM, IP_INS_DATE from INV_PLATE "
                         "where IP_BARCODE in :g or IP_GRID_NUM in :g").bindparams(
                bindparam('g', expanding=True))
            in_hits = pd.read_sql(query, cnx, params={'g': sorted(set(grids))})
        finally:
            if close:
                cnx.dispose()
        if len(in_hits):
            rows = '\n'.join(f"  {r.IP_BARCODE} (grid_num {r.IP_GRID_NUM}, added {r.IP_INS_DATE})"
                             for r in in_hits.itertuples())
            logger.warning(f"{len(in_hits)} destination plate(s) already exist in HiTS inventory ({READ_DB}):\n"
                           f"{rows}\nIs this a rerun of an earlier upload? If not, the barcode is already "
                           "taken: give the plate a new Grid_num.")
    return {'duplicated': dup, 'in_hits': in_hits}


def _source_plate_mode(value):
    """Which Source_plate rule a cell selects.

    Returns 'blank', 'na', 'biomek', 'lot', 'map' or 'barcodes'. Keywords only
    count when the cell starts with 'use ' or is the bare word ('lots',
    'map', ...), so a real barcode that happens to contain 'MAP' or 'lot' is
    still read as a barcode. 'NA', 'n/a', 'n.a.' and 'none' mean no source
    plate.
    """
    v = _safe_lower(value).strip()
    if v in ('', 'nan'):
        return 'blank'
    if v in NO_SOURCE_PLATE_WORDS or v.startswith('n.a.'):
        return 'na'
    if v.startswith('use ') or v in ('lot', 'lots', 'map', 'biomek'):
        if 'biomek' in v:
            return 'biomek'
        if 'lot' in v:
            return 'lot'
        if 'map' in v:
            return 'map'
        raise ValueError(f"Unknown Source_plate keyword {value!r}. {SOURCE_PLATE_HELP}")
    return 'barcodes'


def _blank_source_plate_message(grid_nums):
    return f"Source_plate is blank for {', '.join(grid_nums)}. {SOURCE_PLATE_HELP}"


def _check_source_plates(info):
    """Raise before anything is written if any Info row has no Source_plate."""
    if 'Source_plate' in info.columns:
        values = info['Source_plate']
    else:  # the column was empty in every row, so _read_info_sheet dropped it
        values = pd.Series([np.nan] * len(info), index=info.index)
    modes = values.map(_source_plate_mode)  # also rejects unknown 'use ...' keywords up front
    blank = modes == 'blank'
    if blank.any():
        raise ValueError(_blank_source_plate_message(info.loc[blank, 'Grid_num'].astype(str).tolist()))


def _assay_dir_name(raw_assay_id, assay_name):
    """Build the <assay_id>_<assay_name> folder name for one Info row.

    Single source of truth: clean_up_folders used str(raw_assay_id) while the
    pipeline used str(int(raw_assay_id)) and had an extra 'map' branch, so an ID
    that stringified as '123.0' in one and '123' in the other produced folders
    that could never be found again.

    Args:
        raw_assay_id: the Raw_assay_ID cell. The literal 'map' means assay IDs
            vary per well, in which case the folder is named for the assay alone.
        assay_name (str): the Assay_name cell.
    Returns:
        str: folder name.
    """
    if 'map' in str(raw_assay_id):
        return str(assay_name)
    return str(int(float(raw_assay_id))) + '_' + str(assay_name)


def _iter_load_groups(info):
    """Yield (grouped_df, rpt) for each raw data file in the Info sheet.

    Single source of truth for the toload_<i> numbering. Both the pipeline and
    clean_up_folders iterate through this, so `i` means the same thing to both.
    They previously grouped differently (by Raw_data_file vs by assay+date),
    which is how clean-up could fail to find a folder that plainly existed.

    Args:
        info (DataFrame): normalised Info sheet from _read_info_sheet.
    Yields:
        tuple: (df, rpt) where df is one file's rows collapsed to one row per
            (assay, name, date, file) with the other fields ';'-joined, and rpt
            is True when the raw file is a .rpt (one report covering the whole
            run, so all its plates share toload_0).
    """
    for file in info.Raw_data_file.unique():
        df = info[(info.Raw_data_file == file)].reset_index(drop=True)
        rpt = str(file).endswith('.rpt')
        cols = ['Raw_assay_ID', 'Assay_name', 'Date_of_screen', 'Raw_data_file']
        if 'Assay_name' not in df.columns:
            df['Assay_name'] = ''
        # dict.fromkeys dedupes while preserving order, so Grid_num becomes the
        # ';'-separated plate list the caller unpacks.
        df = df.astype(str)
        df = df.groupby(cols).agg(func=lambda x: ';'.join(list(dict.fromkeys(x)))).reset_index()
        yield df, rpt


###############################################################################
# SECTION 6: the dose-response pipeline
#
# make_dr_platemaps resolves ONE destination plate. It is orchestration only:
# each get_*_info helper below owns one column of the output and is called in
# dependency order, because later steps read columns earlier ones added.
###############################################################################
def make_dr_platemaps(dp, row, users, basedir=None, replace_dict=None, control_dict=None):
    """ Resolve one destination plate to a full, loadable well manifest.

    Args:
        dp (dict): {sheet_name: DataFrame} for the whole workbook, as returned
            by pd.read_excel(..., sheet_name=None).
        row (Series): one row of the Info sheet, naming this plate's Grid_num,
            Source_plate, Raw_assay_ID and Lot_info.
        users (list of tuples): [(FirstName, LastName)], used only by the
            last-resort lot lookup.
        basedir: accepted and ignored; kept for call-site compatibility.
        replace_dict (dict, optional): extra control/solvent labels, see
            _build_replace_metadata.
        control_dict (dict, optional): labels used in block 0 for scientist controls.
            Use {'Sci_PosControl': '<positive label>', 'Sci_NegControl': '<negative label>'}.
            Matching wells are emitted as Sci_PosControl / Sci_NegControl booleans,
            and the labels are automatically treated as non-compound SMDC_ID values.
    Returns:
        DataFrame: one row per destination well with columns SMDC_ID, Lot,
            Source plate, Source well, Destination plate, Destination well,
            [compound] uM, Assay ID. Missing values are filled with the string
            'NA', which is what the HiTS loader expects. Wells with zero
            concentration are dropped, i.e. treated as empty wells.

    Non-interactive.
    """
    _ = basedir  # accepted but unused; silences linters about the dead param
    # Control labels are annotations, not SMDC IDs. Add them to the replacement
    # map automatically so callers do not have to repeat them in replace_dict.
    # Keys are normalised only to recognise them ('sci_poscontrol' -> 
    # 'Sci_PosControl'); labels are compared normalised in
    # get_destination_plate_info but replaced using the text as typed.
    normalized_controls = _normalize_control_dict(control_dict)

    replace_dict = dict(replace_dict or {})
    for label in (control_dict or {}).values():
        if label is not None and str(label).strip():
            replace_dict[str(label)] = np.nan
            replace_dict[str(label).strip()] = np.nan
    replace_dict, replace_vals = _build_replace_metadata(replace_dict)
    logger.info(f"Labels replaced in block 0: {replace_dict}")

    # Open one connection and thread it through the lookups. A missing
    # credential raises here, next to its cause, instead of failing later.
    cnx = _connect()

    grid_num = str(_series_get_ci(row, 'Grid_num', ''))
    source_plate_val = _safe_lower(_series_get_ci(row, 'Source_plate', ''))
    logger.info(grid_num)

    # Two shapes of input. Biomek worklists are already long format and already
    # name their source wells, so they skip block parsing entirely and only need
    # lots looked up. Everything else is a hand-drawn block grid.
    if _source_plate_mode(source_plate_val) == "biomek":
        logger.info("Source plate is long format")
        dids = dp[grid_num[0:31]].copy()  # [0:31]: Excel truncates sheet names
        dids = _standardize_dids_columns(dids)
        dids['row'] = dids['Destination well'].str.extract(r'([A-Za-z]+)')[0]
        dids['col'] = dids['Destination well'].str.extract(r'(\d+)')[0].astype(int)
        maps = None
    else:
        # Order matters: each step consumes columns the previous one added.
        dids, maps = get_destination_plate_info(dp, row, replace_dict,
                                                normalized_controls)  # block 0 -> SMDC_ID + controls
        dids = get_concentration_info(dids, maps)  # block 1 -> [compound] uM
        dids = get_source_plate_info(row, dids, maps, cnx, replace_vals)  # -> Source plate/well
        dids = _standardize_dids_columns(dids)

    dids['Source plate'] = dids['Source plate'].str.upper()
    dids = get_assay_id_info(dids, maps, raw_assay_id=row['Raw_assay_ID'])
    dids = get_lot_info(row, dids, maps, users, replace_vals)
    dids = dids.rename(columns={'Concentration': '[compound] uM', 'Concentration_uM': '[compound] uM',
                                'concentration': '[compound] uM', 'conc': '[compound] uM'})
    # Keep the scientist control annotations in the HiTS manifest. They are true
    # only for wells whose block-0 text matches a control_dict label (ignoring
    # case, spaces and punctuation).
    # Biomek worklists skip get_destination_plate_info, which is what normally
    # adds these columns, so default them to False when they are missing.
    dids = _ensure_control_flags(dids)
    dids = dids[['SMDC_ID', 'Lot', 'Source plate', 'Source well', 'Destination plate',
                 'Destination well', '[compound] uM', 'Assay ID',
                 'Sci_PosControl', 'Sci_NegControl']]

    if cnx is not None:
        cnx.dispose()

    dids = _standardize_dids_columns(dids)
    dids = dids[dids['Destination well'].notna()]
    dids = dids.drop_duplicates()

    # Sanity check on row count. A plate should have exactly one row per well
    # (320 = a 384w plate with the 4 control columns excluded). Anything else
    # usually means a merge fanned out because a compound was found on more
    # than one source plate.
    if len(dids) not in [96, 384, 1536, 320]:
        logging.warning(f"Your plate has {len(dids)} wells. Was a compound pulled from two plates?")

    # Only now drop the empty (zero-concentration) wells. Doing this after the
    # count above means the count sees the full plate, so the warning stays
    # meaningful instead of firing on every plate with empty wells.
    dids = handle_zero_concentrations(dids)
    dids = _round_concentrations(dids)

    # HiTS wants the literal string 'NA' for empty fields, not blanks or NaN.
    dids = dids.fillna('NA').replace("-", "NA").replace('water', 'NA').replace('Water', 'NA')
    return dids


def get_source_plate_info(row, dids, maps, cnx, replace_vals):
    """ Resolve Source plate / Source well for every destination well.

    THE Source_plate CELL IS A DISPATCH LANGUAGE (see _source_plate_mode).
    Keywords count only as 'use ...' or the bare word, so a barcode that
    contains 'MAP' or 'lot' is still a barcode. Case-insensitive:

        'NA', 'n/a', 'n.a.', 'none'
                        no source plate exists (compounds added by hand).
                        Source plate/well are literally 'NA'; lots come from
                        the newest lot registered by `users`.
        'use lots'      lots are drawn in block 2 of the map; plates are then
                        looked up from (SMDC_ID, Lot).
        'use map'       source plate barcodes are drawn in block 2 of the map,
                        i.e. wells came from different plates.
        'use biomek'    handled upstream in make_dr_platemaps, not here.
        blank           an error. Blank used to mean 'use lots', which hid a
                        forgotten Source_plate until much later.
        anything else   treated as one or more literal barcodes, separated by
                        ';' and/or newlines. Every listed plate is pulled from
                        HiTS and joined on SMDC_ID, so a compound present on
                        two listed plates fans out into duplicate rows (this is
                        what the well-count warning catches).

    Args:
        row (Series): one Info row.
        dids (DataFrame): destination wells so far; must have SMDC_ID.
        maps (DataFrame): this plate's raw map sheet, for block reads.
        cnx: open connection, reused across barcode lookups.
        replace_vals (list): control/solvent values to skip.
    Returns:
        DataFrame: dids with Source plate and Source well (and sometimes Lot).

    NOTE: `source_plates` is only bound in some branches, so the
    `if len(source_plates) > 0` below relies on the branch above having set it.
    The 'n.a.' branch escapes this only because it sets 'Source plate' first.
    """
    dids = _standardize_dids_columns(dids)
    source_plate_raw = _series_get_ci(row, 'Source_plate', '')
    source_plate_val = _safe_lower(source_plate_raw)
    mode = _source_plate_mode(source_plate_raw)

    if mode == 'na':
        logger.info("Source plate is N.A.")
        dids['Source plate'] = 'NA'
        dids['Source well'] = 'NA'
    elif mode == 'blank':
        raise ValueError(_blank_source_plate_message(
            [str(_series_get_ci(row, 'Grid_num', '?'))]))
    elif mode == 'lot':
        logger.info("Will get source plate info from lot lookup")
        try:
            lots = read_platemap_block(plate_map=maps, block=2, value_name='Lot')
        except MissingBlockError:
            # Fine when every compound well was pinned as '12346-1' in block 0.
            if 'Lot' in dids.columns and not dids.loc[dids.SMDC_ID.notna(), 'Lot'].isna().any():
                lots = None
            else:
                raise
        if lots is None:
            pass
        elif 'Lot' in dids.columns:
            lots['Destination well'] = lots.row + lots.col
            # Some lots were pinned with '12346-1' in block 0. Those win; the
            # other wells take their lot from block 2.
            block2 = dids['Destination well'].map(lots.set_index('Destination well')['Lot'])
            dids['Lot'] = dids['Lot'].where(dids['Lot'].notna(), block2)
        else:
            lots['Destination well'] = lots.row + lots.col
            dids = dids.merge(lots[['Destination well', 'Lot']], how='left')
        source_plates = []
    elif mode == 'map':
        logger.info("Source plate is map")
        sps = read_platemap_block(plate_map=maps, block=2, value_name='Source plate')
        sps['Destination well'] = sps.row + sps.col
        dids = dids.merge(sps[['Destination well', 'Source plate']], how='left')
        source_plates = sps[~sps['Source plate'].isna()]['Source plate'].astype(str).str.strip().tolist()
    else:
        logger.info("Getting source plates from listed barcodes")
        # logger.debug(f"Raw source plate value: {source_plate_raw}")
        # Barcodes may be ';'-separated, newline-separated, or both, because
        # they are pasted in by hand. Split on both, strip, dedupe.
        source_plates = str(source_plate_raw).split(';')
        final_source_plates = []
        for x in source_plates:
            final_source_plates.extend(x.splitlines())
        source_plates = [x.strip() for x in final_source_plates if x.strip()]
        source_plates = list(set(source_plates))
        # logger.debug(f"Extracted source plates: {source_plates}")

    if 'Source plate' not in dids.columns:
        if len(source_plates) > 0:
            logger.info("Getting source plate info from listed barcodes")
            dfs = []
            for plate in source_plates:
                logger.info(f"Processing source plate {plate}")
                df = make_platemaps_from_barcode(barcode=plate, outdir=None, cnx=cnx)
                df = _standardize_dids_columns(df)
                df = df[~df.SMDC_ID.isna()]
                df.SMDC_ID = df.SMDC_ID.astype(float)
                dfs.append(df)
            df = pd.concat(dfs)
            # Rename by name, not by position. The previous
            # `df.columns = [...]` assignment silently mislabelled every column
            # if make_platemaps_from_barcode's SELECT list changed order. This
            # raises instead, and names the column that went missing.
            rename_map = {
                'Well': 'Source well',
                'Vendor_Alias': 'alias',
                'Common_Name': 'Common Name',
                'Lot_Num': 'Lot',
                'Plate': 'Source plate',
            }
            missing_cols = set(rename_map) - set(df.columns)
            if missing_cols:
                raise KeyError(
                    f"make_platemaps_from_barcode returned unexpected columns: "
                    f"missing {sorted(missing_cols)}. Got {list(df.columns)}. "
                    f"Did its SELECT list change?"
                )
            df = df.rename(columns=rename_map)
            dids = dids.merge(df[['SMDC_ID', 'Lot', 'Source plate', 'Source well']], how='left')
        elif 'Lot' in dids.columns:
            logger.info("Getting source plate info from lot lookup")
            lookups = dids[['SMDC_ID', 'Lot']].drop_duplicates()
            try:
                plates = lookup_plate_by_lot_num(
                    ids=lookups.SMDC_ID.astype(str).tolist(),
                    lot_nums=lookups.Lot.astype(str).tolist(),
                    like_filter=['Source plate', "MA01"],
                    replace_vals=replace_vals
                )
                dids = dids.merge(plates, how='left')
            except (DBAPIError, ConnectionError):
                raise  # a database/connection error must stop the run, not become 'NA'
            except Exception as e:
                logger.warning(f"Error looking up plates by lot num: {e}")
                dids['Source plate'] = 'NA'
                dids['Source well'] = 'NA'
        else:
            logger.info("No source plate info available, will lookup after finding lots")
            dids['Source plate'] = 'NA'
            dids['Source well'] = 'NA'

    return _standardize_dids_columns(dids)


def get_lot_info(row, dids, maps, users, replace_vals):
    """ Resolve a Lot for every well that has an SMDC_ID, via a fallback ladder.

    The lot is what makes a result reproducible: same compound, different
    synthesis batch. It is also the field most often missing from a hand-drawn
    map, hence the cascade. Each stage only fills wells the previous ones left
    empty, in decreasing order of confidence:

      1. If Lot is already present (drawn in the map, or split from a
         'SMDC-Lot' string), keep it.
      2. Lot_info == 'use map'  -> read block 2 as Lot.
         Source_plate ~ 'map'   -> look up by (SMDC_ID, source plate).
         Source_plate ~ 'biomek'-> look up by (source plate, source well); this
                                   branch also REPLACES SMDC_ID with the
                                   database's answer, since the worklist's own
                                   ID column is not trusted.
      3. Still missing -> look up by (SMDC_ID, any known source plate).
      4. Still missing -> newest lot registered by `users`. An assumption.
      5. Still missing -> leave NaN, warn, and let it become 'NA' downstream.

    Args:
        row (Series): one Info row.
        dids (DataFrame): destination wells, must have SMDC_ID.
        maps (DataFrame): this plate's raw map sheet.
        users (list of tuples): [(FirstName, LastName)] for stage 4.
        replace_vals (list): control/solvent values to skip.
    Returns:
        DataFrame: dids with a Lot column.

    NOTE: stage 4 maps by SMDC_ID with no plate constraint, so a compound whose
    lot could not be traced gets this person's newest lot attached regardless of
    whether that lot was ever on the plate. Silent, and only a warning if even
    that fails. If lot provenance matters for downstream QC, this stage is the
    one to interrogate.
    """
    dids = _standardize_dids_columns(dids)
    source_plate_val = _safe_lower(_series_get_ci(row, 'Source_plate', ''))
    lot_info_val = _safe_lower(_series_get_ci(row, 'Lot_info', ''))

    if 'Lot' in dids.columns and dids.Lot.isna().any() and maps is not None and (
            lot_info_val == 'use map' or _source_plate_mode(source_plate_val) == 'map'):
        # Some wells were pinned with '12346-1' in block 0; fill the rest from
        # stage 2 as if nothing had been pinned. Pinned lots always win.
        logger.info("Some lots pinned in block 0; filling the rest from the map settings")
        pinned = dids.dropna(subset=['Lot']).set_index('Destination well')['Lot']
        rest = get_lot_info(row, dids.drop(columns='Lot'), maps, users, replace_vals)
        from_pin = rest['Destination well'].map(pinned)
        rest['Lot'] = from_pin.where(from_pin.notna(), rest['Lot'])
        return rest

    if 'Lot' not in dids.columns:
        if (lot_info_val == 'use map') and (maps is not None):
            logger.info("Lot info is use map")
            aids = read_platemap_block(plate_map=maps, block=2, value_name='Lot')
            aids['Destination well'] = aids.row + aids.col
            dids = dids.merge(aids[['Destination well', 'Lot']], how='left')
        elif _source_plate_mode(source_plate_val) == 'map':
            logger.info("Lot info from plate map")
            ids, plates = _id_plate_pairs(dids)
            lots = lookup_lot_by_plate_and_ic_id(ids=ids, plates=plates, replace_vals=replace_vals, matched=True)
            lots.SMDC_ID = lots.SMDC_ID.astype(float)
            dids = dids.merge(lots, how='left')
        elif _source_plate_mode(source_plate_val) == 'biomek':
            logger.info("Lot info from source plate and well")
            lots = lookup_id_and_lot_by_plate_and_well(plates=dids['Source plate'].tolist(),
                                                       wells=dids['Source well'].tolist())
            lots.SMDC_ID = lots.SMDC_ID.astype(float)
            dids = dids.drop(columns='SMDC_ID').merge(lots, how='left')
        else:
            dids['Lot'] = np.nan

    # Stage 3: recover missing lots from source plate + SMDC_ID.
    # `missing_mask` is recomputed from scratch before each stage; that is the
    # loop-free way this cascade narrows down.
    missing_mask = (~dids.SMDC_ID.isna()) & ((dids.Lot.isna()) | (dids.Lot == ''))
    if len(dids[missing_mask]) > 0:
        logger.info('Getting lot info from source plates and SMDC IDs')
        ids, plates = _id_plate_pairs(dids[missing_mask])
        if len(plates) == 0:
            logger.info("No source plate info available, will lookup lots by user")
        else:
            lots = lookup_lot_by_plate_and_ic_id(ids=ids, plates=plates, replace_vals=replace_vals, matched=True)
            lots.SMDC_ID = lots.SMDC_ID.astype(float)
            # Preserve non-compound rows (controls/solvents). Their SMDC_ID was
            # deliberately replaced with NaN earlier, but they still belong in the
            # destination manifest -- especially scientist control wells carrying
            # Sci_PosControl/Sci_NegControl=True.  The old split only retained rows
            # with a non-null SMDC_ID, silently deleting every annotated control.
            non_compounds = dids[dids.SMDC_ID.isna()]
            nolots = dids[(dids.SMDC_ID.notna()) & (dids.Lot.isna())]
            haslots = dids[(dids.SMDC_ID.notna()) & (dids.Lot.notna())]
            # Split, fill only the empty compound half, then recombine with the
            # untouched non-compound rows. NOTE: the lookup merge can still fan
            # out if a compound has more than one lot on the plate.
            nolots = nolots.drop(columns=['Lot']).merge(lots, how='left')
            dids = pd.concat([haslots, nolots, non_compounds], ignore_index=True)
    else:
        logger.info("There are no missing lots or SMDC IDs are all NaN")

    # Stage 4: last resort, assume the newest lot this user registered.
    missing_mask = (~dids.SMDC_ID.isna()) & ((dids.Lot.isna()) | (dids.Lot == ''))
    if len(dids[missing_mask]) > 0:
        ids = dids[missing_mask].SMDC_ID.astype(int).astype(str).unique().tolist()
        # Expected with 'NA' / 'use my name'; anywhere else it's a guess the operator should check.
        level = (logging.INFO if _source_plate_mode(_series_get_ci(row, 'Source_plate', '')) == 'na'
                 else logging.WARNING)
        logger.log(level, f"Getting lot info from user for {len(ids)} compound(s) on "
                          f"{_series_get_ci(row, 'Grid_num', '')}: the newest lot registered by USERS. "
                          "Check those wells.")
        lot_map = _unique_user_lots(lookup_lot_by_user(ids=ids, users=users))
        mask = dids["SMDC_ID"].notna() & dids["Lot"].isna()
        dids.loc[mask, "Lot"] = dids.loc[mask, "SMDC_ID"].map(lot_map)

    # Stage 5: give up. Normalise '' to NaN so the fillna('NA') at the end of
    # make_dr_platemaps catches it.
    missing_mask = (~dids.SMDC_ID.isna()) & ((dids.Lot.isna()) | (dids.Lot == ''))
    if len(dids[missing_mask]) > 0:
        logger.warning("no Lot from plate map")
        dids.loc[(~dids.SMDC_ID.isna()) & (dids.Lot.isna()), 'Lot'] = np.nan
        dids.loc[(~dids.SMDC_ID.isna()) & (dids.Lot == ''), 'Lot'] = np.nan

    return _standardize_dids_columns(dids)


def get_assay_id_info(dids, maps, raw_assay_id):
    """Resolve assay ID field from map or single value.

    Args:
        dids (DataFrame): destination wells, must have 'Destination well'.
        maps (DataFrame): this plate's raw map sheet.
        raw_assay_id: the Info sheet's Raw_assay_ID cell. If it contains 'map',
            assay IDs vary by well and are read from block 2; otherwise it is a
            single ID applied to the whole plate.
    Returns:
        DataFrame: dids with an 'Assay ID' column.

    NOTE: block 2 is shared between Lot, Source plate and Assay ID, so a plate
    cannot have both per-well assay IDs and a per-well third field of another
    kind. That constraint is implicit in the workbook format.
    """
    if pd.isna(raw_assay_id):
        dids['Assay ID'] = np.nan
    elif 'map' not in _safe_lower(raw_assay_id):
        logger.info("Assay ID is not multiple")
        dids['Assay ID'] = int(float(raw_assay_id))
    else:
        logger.info("Assay ID is multiple")
        aids = read_platemap_block(plate_map=maps, block=2, value_name='Assay ID')
        aids['Destination well'] = aids.row + aids.col
        dids = dids.merge(aids[['Destination well', 'Assay ID']], how='left')

    return dids


def get_concentration_info(dids, maps):
    """Resolve concentration field (block 1) and attach to destination wells.

    Args:
        dids (DataFrame): destination wells, must have 'Destination well'.
        maps (DataFrame): this plate's raw map sheet.
    Returns:
        DataFrame: dids with a float '[compound] uM' column.

    Rounded to 7 dp so that concentrations produced by repeated serial-dilution
    arithmetic in Excel compare equal across plates instead of differing in the
    last float bit.
    """
    # Block 1 is optional. Cherry-pick maps may have no block 1, or use it for
    # a plate suffix instead of concentrations; the suffix and a whole-plate
    # concentration are then passed in during plate registration. In both cases
    # leave the concentration blank here, as the old fixed-stride reader did.
    try:
        dums = read_platemap_block(plate_map=maps, block=1, value_name='[compound] uM')
    except MissingBlockError:
        logger.info("No block 1 in the plate map; leaving concentrations blank")
        dids['[compound] uM'] = np.nan
        return dids
    values = dums['[compound] uM'].dropna()
    values = values[values.astype(str).str.strip() != '']
    numeric = pd.to_numeric(values, errors='coerce')
    if len(values) and numeric.isna().all():
        logger.info("Block 1 holds text (for example a plate suffix), not "
                    "concentrations; leaving concentrations blank")
        dids['[compound] uM'] = np.nan
        return dids
    if numeric.isna().any():
        bad = sorted(values[numeric.isna()].astype(str).unique())[:5]
        raise ValueError(f"Block 1 mixes concentrations with text values {bad}. "
                         "Fix those wells in the plate map.")
    dums['Destination well'] = dums.row + dums.col
    dids = dids.merge(dums[['Destination well', '[compound] uM']], how='left')
    dids['[compound] uM'] = pd.to_numeric(dids['[compound] uM'], errors='coerce').round(7)
    logger.info("Read concentrations and rounded to 7 decimal places")
    return dids


def get_destination_plate_info(dp, row, replace_dict, control_dict=None):
    """Read block 0 (SMDC IDs) and set up the destination-well frame.

    Args:
        dp (dict): {sheet_name: DataFrame} for the workbook.
        row (Series): one Info row; supplies Grid_num (= destination barcode
            AND sheet name).
        replace_dict (dict): control/solvent labels to null out.
        control_dict (dict, optional): {'Sci_PosControl': label, 'Sci_NegControl': label}
            with labels already normalised by _normalize_control_dict (as
            make_dr_platemaps passes them).
    Returns:
        tuple: (dids, maps) - the long-form destination frame with SMDC_ID,
            Destination well and Destination plate, plus the raw map sheet for
            later block reads.

    Handles the 'SMDC-Lot' shorthand: if any ID cell contains a dash, the whole
    column is split into SMDC_ID and Lot, which lets a scientist pin an exact
    lot in the map itself.
    """
    grid_num = str(_series_get_ci(row, 'Grid_num', ''))
    raw_assay_id = _series_get_ci(row, 'Raw_assay_ID', np.nan)

    # Sheet name == destination plate barcode, truncated to Excel's 31-char cap.
    maps = dp[grid_num[0:31]]
    maps = maps.dropna(how='all')

    dids = read_platemap_block(plate_map=maps, block=0, value_name='SMDC_ID')

    # Annotate controls BEFORE replacing text labels with NaN. Exact string
    # matching is intentional: the caller supplies the text used in the map.
    control_dict = dict(control_dict or {})
    pos_label = control_dict.get('Sci_PosControl')
    neg_label = control_dict.get('Sci_NegControl')
    # Compare normalized text rather than raw Python strings. This handles harmless
    # differences such as case and surrounding whitespace in hand-entered maps.
    map_labels = dids.SMDC_ID.map(_norm_key)
    dids['Sci_PosControl'] = False if pos_label is None else map_labels.eq(pos_label)
    dids['Sci_NegControl'] = False if neg_label is None else map_labels.eq(neg_label)

    for kind, label, col in [
        ('positive', pos_label, 'Sci_PosControl'),
        ('negative', neg_label, 'Sci_NegControl'),
    ]:
        if label is not None:
            n_matches = int(dids[col].sum())
            if n_matches == 0:
                observed_text = sorted({
                    str(v).strip() for v in dids.SMDC_ID.dropna()
                    if not str(v).strip().replace('.', '', 1).isdigit()
                })
                logger.warning(
                    f"Configured {kind} control label {label!r} matched 0 wells on "
                    f"{grid_num}. Text annotations seen in block 0: {observed_text}"
                )
            else:
                logger.info(
                    f"Annotated {n_matches} {kind} control well(s) on {grid_num} "
                    f"from label {label!r}."
                )
    # Control wells are annotations, not compounds. Null them here using the
    # same normalised comparison as the flags, so 'no enzyme ' or 'NO ENZYME'
    # in the map does not reach the float cast below.
    control_labels = {v for v in (pos_label, neg_label) if v}
    if control_labels:
        dids.loc[map_labels.isin(control_labels), 'SMDC_ID'] = np.nan

    # Null out (or remap) control/solvent labels so the column can become float.
    # Done before the ID-lot split so a label with a dash in it ('No-Enzyme')
    # is not mistaken for one.
    for key, value in replace_dict.items():
        dids.SMDC_ID = dids.SMDC_ID.replace(key, value)

    # '12346-1' pins SMDC ID 12346 to lot 1. Only cells that are exactly
    # <digits>-<digits> (spaces allowed) are split; other wells keep Lot empty
    # and get their lot from the usual lookups in get_lot_info.
    id_lot = dids.SMDC_ID.astype(str).str.extract(r'^\s*(\d+)\s*-\s*(\d+)\s*$')
    pinned = id_lot[0].notna()
    if pinned.any():
        logger.info(f"Read {int(pinned.sum())} SMDC ID-lot cell(s) like 12346-1")
        dids['Lot'] = np.nan
        dids.loc[pinned, 'SMDC_ID'] = id_lot.loc[pinned, 0]
        dids.loc[pinned, 'Lot'] = id_lot.loc[pinned, 1].astype(float)

    # An unrecognised text label (or a malformed ID like '12346-1-2') reaches
    # this cast and raises a ValueError naming the offending string, which is
    # the intended failure mode: fix the map or add the label to replace_dict.
    dids.SMDC_ID = dids.SMDC_ID.astype(float)

    dids['Destination well'] = dids.row + dids.col
    dids['Destination plate'] = grid_num
    logger.info("Read SMDC ID's")

    return _standardize_dids_columns(dids), maps


CONC_SIG_FIGS = 10


def _round_concentrations(dids, sig_figs=CONC_SIG_FIGS):
    """Round '[compound] uM' to sig_figs significant figures.

    The same dilution often reaches the map in two forms: a formula result
    (83.33333333333333) in some cells and a typed or pasted value
    (83.3333333333333) in others. They print the same but are different
    numbers, so Pipeline Pilot's 'Insert Conc Variable' sees two values for
    one concentration and fails with "Cannot insert multiple values for
    HLD_VALUE". Ten significant figures is far beyond any real dilution.
    """
    col = '[compound] uM'
    if col not in dids.columns:
        return dids
    conc = pd.to_numeric(dids[col], errors='coerce')
    ok = conc.notna() & np.isfinite(conc) & (conc != 0)
    if ok.any():
        dids = dids.copy()
        dids.loc[ok, col] = conc[ok].map(lambda v: float(f'{v:.{sig_figs}g}'))
    return dids


def handle_zero_concentrations(dids):
    """Drop zero-concentration wells, treating them as empty wells.

    A zero in the concentration column is not loadable: dose-response curves are
    fit on a log scale, so zero has no place on the axis. In practice a zero
    means the well held no compound, so the well is simply dropped from the
    manifest, which is how an empty well is represented.

    This used to prompt on stdin for a replacement value, which is why the
    pipeline could not run unattended. Pressing Enter or typing NA at that prompt
    produced exactly this behaviour, and that was the recommended answer, so the
    prompt is gone and the recommended answer is now the only one.

    If you ever have zeros that genuinely mean "below the lowest dilution"
    rather than "empty", fix them in the workbook. Do not reintroduce the prompt.

    Args:
        dids (DataFrame): destination wells with a concentration column.
    Returns:
        DataFrame: with zero-concentration rows removed.

    Call this AFTER the plate row-count check in make_dr_platemaps, so the count
    is taken on the full plate. Dropping rows first made the check fire a false
    "was a compound pulled from two plates?" warning on every plate that had any
    empty wells.
    """
    # Support both 'conc' and '[compound] uM' column names
    conc_col = '[compound] uM' if '[compound] uM' in dids.columns else 'conc'

    if conc_col not in dids.columns:
        return dids

    # Control wells (Sci_PosControl / Sci_NegControl) hold no compound, so they
    # have no concentration: blank it (written as NA), whatever the map says,
    # and keep the well, or the control flags never reach the manifest.
    controls = pd.Series(False, index=dids.index)
    for flag in ('Sci_PosControl', 'Sci_NegControl'):
        if flag in dids.columns:
            controls |= dids[flag].fillna(False).astype(bool)
    if controls.any():
        dids = dids.copy()
        dids.loc[controls, conc_col] = np.nan
    zero_mask = (dids[conc_col] == 0) & ~controls
    n_zero = int(zero_mask.sum())
    if n_zero == 0:
        return dids

    logger.info(
        f"Dropping {n_zero} well(s) with zero concentration; treating them as empty wells."
    )
    dids = dids[~zero_mask]
    if len(dids) > 0:
        logger.info(f"Min conc of plate is now {dids[conc_col].min()}")
    return dids


###############################################################################
def make_dr_platemaps_and_folders(infile='', users=None, basedir='./', dose_response=True,
                                  replace_dict=None, control_dict=None, share_root=None):
    """
    Create an excel file with dose response plate information for DR assay upload
    to HiTS. Each sheet in the input file should be labeled with the destination
    plate barcode and have SMDC IDs and Concs. There should be a sheet labeled 'Info'
    with additional data including source plate barcodes, assay IDs and dates.
    Args:
        infile (str): Path to the input file containing plate information.
        users (list of tuples): List of user names in the form of (FirstName, LastName),
            used only for the last-resort lot lookup. No default: if that lookup is
            needed and users is empty, the run stops with a ValueError.
        basedir (str): Place to move data to and create folders in.
        dose_response (bool): Whether to create dose response platemaps. If
            False, folders are still created and raw data still copied; only
            the per-plate manifests are skipped.
        replace_dict (dict, optional): Extra control/solvent labels; see
            _build_replace_metadata.
        control_dict (dict, optional): Control labels written in block 0, e.g.
            {'positive': 'Positive Control', 'negative': 'Negative Control'}.
            These labels are converted to Sci_PosControl / Sci_NegControl boolean
            columns in each *_DR_plate.xlsx output.
        share_root (str, optional): where the Montreal share is mounted, e.g.
            'Z:/', '/mnt/z/' or '/Volumes/Shared/'. Default: default_share_root().
    Returns:
        info (DataFrame): the parsed Info sheet, as a record of what was done.
            load_folders(info) lists the toload_<i>/<date> folders it wrote.

    Side effects, which are the actual point of this function:
        * creates <basedir>/<assay_id>_<assay_name>/toload_<i>/<date>/
        * writes <grid_num>_DR_plate.xlsx into it (dose_response only)
        * writes/updates grid_num.txt listing the plates in that folder
        * copies the raw instrument file in from the Montreal share
          (share_root, default: default_share_root())

    Interactive: make_dr_platemaps may prompt on stdin, and display() assumes a
    notebook. Grouping is per Raw_data_file, since one instrument file may cover
    several plates.

    NOTE the `i` semantics: for .rpt files i is pinned to 0 so everything lands
    in toload_0, because one .rpt covers the whole run. For other formats i is
    the enumerate index, giving toload_0, toload_1, ... per row. The numbering
    and the folder name both come from _iter_load_groups / _assay_dir_name, which
    clean_up_folders also uses, so the two cannot drift apart.
    """

    # read plate layout file (shared with clean_up_folders)
    share_root = share_root or default_share_root()
    infile = resolve_share_path(infile, share_root)
    basedir = resolve_share_path(basedir, share_root)
    dp, info = _read_info_sheet(infile)
    check_screen_dates(info, share_root=share_root)
    if dose_response:
        _check_source_plates(info)  # fail now, not after folders are made and files copied
    check_destination_plates(info)
    load_folders = []
    dr_plates = []  # (grid_num, path) of every *_DR_plate.xlsx written, for check_dr_plates
    raw_files = []  # raw instrument files copied into the toload folders, for tidying

    # One instrument file may back several plates, so group by file.
    for df, rpt in _iter_load_groups(info):
        # display(df)  # NOTE: notebook-only builtin; NameError if run as a script

        for i, row in df.iterrows():
            if rpt:
                i = 0
            row.Date_of_screen = str(int(row.Date_of_screen))
            # create folders
            outdir = os.path.join(basedir, _assay_dir_name(row.Raw_assay_ID, row.Assay_name),
                                  f'toload_{i}', row.Date_of_screen)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
                logger.info(f'create folder {outdir}')
            if outdir not in load_folders:
                load_folders.append(outdir)

            # make and save platemaps if dose response
            gridnum = pd.DataFrame(row.Grid_num.split(';'), columns=['grid_num'])
            if dose_response:
                for grid_num in gridnum.grid_num:
                    # The group row joins every plate's Source_plate/Lot_info
                    # with ';'. Use this plate's own Info row instead, so plates
                    # sharing a raw file can still use different settings.
                    plate_row = _plate_info_row(info, row, grid_num)
                    dids = make_dr_platemaps(dp, plate_row, users, replace_dict=replace_dict, control_dict=control_dict)
                    dids = dids.drop_duplicates()

                    # Apply scientist-control flags LAST, after all plate/lot/source logic.
                    dids = _apply_control_flags_from_map(dids, dp, plate_row, control_dict)

                    dr_path = os.path.join(outdir, f'{grid_num}_DR_plate.xlsx')
                    with _FileInUse(dr_path):
                        dids.to_excel(dr_path, index=False)
                    dr_plates.append((grid_num, dr_path))
                row.Grid_num = ';'.join(gridnum.grid_num)

            # create or modify grid_num.txt: the manifest the HiTS loader reads
            # to know which destination plates this folder's data covers.
            # Appended to rather than overwritten, since a folder may be filled
            # across several runs.
            _append_grid_num(outdir, gridnum)

            # move data files to folders
            logger.info('Copying data...')
            starttime = time.time()
            # Path normalisation is shared with clean_up_folders.
            row.Raw_data_file = _normalize_share_path(row.Raw_data_file)
            logger.info(row.Raw_data_file)

            origfile = os.path.join(share_root, row.Raw_data_file)
            fname = row.Raw_data_file.rsplit('/', 1)[-1]
            outfile = os.path.join(outdir, fname)
            logger.info(f'{origfile} -> {outfile}')

            # _copy_file rather than os.system("cp ..."): no shell quoting to
            # get wrong, and a failed copy raises instead of scrolling past.
            # It keeps mtime where the mount allows, which matters for provenance.
            if not os.path.exists(origfile):
                raise FileNotFoundError(
                    f"Raw_data_file for {row.Grid_num} not found: {origfile!r}. Check the Info tab's "
                    f"Raw_data_file path. If the share is mounted somewhere else: {SHARE_ROOT_HELP}")
            _copy_file(origfile, outfile)
            if outfile not in raw_files:
                raw_files.append(outfile)
            _warn_if_envision_blocks_out_of_order(outfile)
            elapsed = time.time() - starttime
            logger.info(f'Done copying in {humanize.precisedelta(elapsed, minimum_unit="seconds", format="%d")}')
    # The toload_<i>/<date> folders this run wrote to, so later steps (the
    # EnVision check) look only at them and not at older assays in basedir.
    info.attrs['load_folders'] = load_folders
    info.attrs['dr_plates'] = dr_plates
    info.attrs['raw_files'] = raw_files
    info.attrs['infile'] = infile
    if dr_plates:
        check_dr_concentrations(info)
    return info


def check_dr_concentrations(info, rel_tol=1e-3):
    """Warn about concentrations that look like the same dilution written two ways.

    Pipeline Pilot's DR upload stores one concentration row per assay and
    value ("Insert Conc Variable"), shared by every plate and compound in that
    assay, after rounding IW_FINAL_CONC_NM to the nearest 0.001 nM (1 fM):
    SUBST_CONC_UM := Floor(IW_FINAL_CONC_NM * 1000 + 0.5) / 1000000.
    Values equal at 1 fM merge cleanly. Values that still differ at 1 fM but
    by less than rel_tol (0.1%) are almost certainly one dilution typed two
    ways (e.g. 83.3333 and 83.33333), and would be stored as two separate
    concentrations. Real dilution steps differ by far more than 0.1%.

    Reads the DR files make_dr_platemaps_and_folders wrote, grouped by Assay ID.

    Returns:
        DataFrame, one row per suspicious pair: Assay ID, the two values in
        uM, and the plates and wells holding each. Empty = nothing to fix.
    """
    dr_plates = info.attrs.get('dr_plates') or []
    frames = []
    for grid_num, path in dr_plates:
        d = pd.read_excel(path)
        if '[compound] uM' not in d:
            continue
        d['[compound] uM'] = pd.to_numeric(d['[compound] uM'], errors='coerce')
        d = d[d['[compound] uM'].notna() & (d['[compound] uM'] > 0)]
        frames.append(d.assign(**{'Destination plate': grid_num}))
    cols = ['Assay ID', 'conc 1 uM', 'conc 2 uM', 'wells with conc 1', 'wells with conc 2']
    if not frames:
        return pd.DataFrame(columns=cols)
    d = pd.concat(frames, ignore_index=True)
    if 'Assay ID' not in d:
        d['Assay ID'] = 'NA'
    d['fM'] = np.floor(d['[compound] uM'] * 1e9 + 0.5).astype('int64')   # what PP stores, in fM
    rows = []
    for assay, g in d.groupby('Assay ID'):
        levels = sorted(g.fM.unique())
        for a, b in zip(levels, levels[1:]):
            if (b - a) / a < rel_tol:
                def where(v):
                    w = g[g.fM == v]
                    return ', '.join(f"{p} {x}" for p, x in zip(w['Destination plate'], w['Destination well']))
                rows.append({'Assay ID': assay, 'conc 1 uM': a / 1e9, 'conc 2 uM': b / 1e9,
                             'wells with conc 1': where(a), 'wells with conc 2': where(b)})
    out = pd.DataFrame(rows, columns=cols)
    if len(out):
        logger.warning(f"{len(out)} pair(s) of concentrations differ by less than {rel_tol:.1%} but are not equal "
                       "to 1 fM, so Pipeline Pilot would store them as separate concentrations. They are probably "
                       "one dilution typed two ways: make them identical in the plate map and run step 2 again.\n"
                       + '\n'.join(f"  assay {r['Assay ID']}: {r['conc 1 uM']:.12g} vs {r['conc 2 uM']:.12g} uM "
                                   f"({r['wells with conc 1']} / {r['wells with conc 2']})" for _, r in out.iterrows()))
    return out


def load_folders(info):
    """The toload_<i>/<date> folders make_dr_platemaps_and_folders wrote to.

    Args:
        info (DataFrame): what make_dr_platemaps_and_folders returned.
    Returns:
        list[str]: one path per folder, in the order they were created. Each
            holds a grid_num.txt, the DR plate files and the copied raw data.
    """
    folders = info.attrs.get('load_folders')
    if folders is None:
        raise ValueError("info has no record of its load folders. Pass the DataFrame that "
                         "make_dr_platemaps_and_folders returned, not a copy or a re-read Info sheet.")
    return list(folders)


def _warn_if_envision_blocks_out_of_order(path):
    """Warn if a copied raw CSV's blocks aren't laid out as Plate information,
    Background information, then results with row A on line 12 of each block
    (raw_data_parsers._block_layout_problems). Best fixed by re-exporting from
    EnVision; otherwise the tidy step re-arranges the copy. Other problems
    (hand-pasted values, reads that differ) are reported by the tidy step.
    Non-CSV files and files with no 'Plate information' section are skipped."""
    if not str(path).lower().endswith('.csv'):
        return []
    try:
        with open(path, newline='') as f:
            text = f.read()
        if 'Plate information' not in text:
            return []
        try:
            from raw_data_parsers import _block_layout_problems, _envision_ncols, plate_dict
        except ImportError:
            return []
        import csv
        nrows = {v[2] + 1: v[1] for v in plate_dict.values()}.get(_envision_ncols(path), 16)
        lines = list(csv.reader(text.splitlines()))
        problems = [p for p in _block_layout_problems(lines, nrows)
                    if 'block starting here' in p or 'before any Plate information' in p]
    except Exception as e:
        logger.warning(f"Could not check the block order of {path}: {e}")
        return []
    if problems:
        logger.warning(f"{path}: the EnVision export is ordered wrong (each block should be Plate information, "
                       f"Background information, then the results, with row A on line 12). Ideally re-export it "
                       f"from EnVision with Plate and Background information placed before the results. If you "
                       f"can't, the tidy step re-arranges this copy.\n  " + '\n  '.join(problems[:5]))
    return problems


def raw_files(info):
    """The raw instrument files make_dr_platemaps_and_folders copied into
    the toload folders (one per Raw_data_file), in the order it copied them."""
    files = info.attrs.get('raw_files')
    if files is None:
        raise ValueError("info has no record of its raw files. Pass the DataFrame that "
                         "make_dr_platemaps_and_folders returned, not a copy or a re-read Info sheet.")
    return list(files)


def _draw_dr_layout(grid_num, dr):
    """Draw a DR file as a plate: SMDC ID in compound wells, POS / NEG in
    control wells, coloured by well type (compound, positive, negative)."""
    from matplotlib.colors import ListedColormap
    d = dr[['Destination well', 'SMDC_ID']].copy()
    pos = dr['Sci_PosControl'].astype(bool) if 'Sci_PosControl' in dr else False
    neg = dr['Sci_NegControl'].astype(bool) if 'Sci_NegControl' in dr else False
    ids = pd.to_numeric(d.SMDC_ID, errors='coerce')
    d['layout: compound / POS / NEG'] = np.where(pos, 2.0, np.where(neg, 3.0, np.where(ids.notna(), 1.0, np.nan)))
    d['label'] = np.where(pos, 'POS', np.where(neg, 'NEG', ids.map(lambda v: '' if pd.isna(v) else str(int(v)))))
    d['label'] = d['label'].replace('', np.nan)
    draw_plate_heatmap(barcode=grid_num, df=d, wells='Destination well',
                       values='layout: compound / POS / NEG', labels='label', plate_size=384, vmin=1, vmax=3, cbar=False,
                       cmap=ListedColormap(['#cfe3f5', '#9be3a3', '#f5a3a3']))


_ID_CELL = re.compile(r'^\s*(\d+)(?:\.0+)?(?:\s*-\s*\d+)?\s*$')


def _expected_smdc_id(cell):
    """The SMDC ID a block-0 cell asks for: 12346, 12346.0 and 12346-1 all give
    12346. Anything else (blank, a control or solvent label) gives NaN."""
    m = _ID_CELL.match(str(cell)) if pd.notna(cell) else None
    return float(m.group(1)) if m else np.nan


def _plates_to_draw(draw, n_plates):
    """How many plates check_dr_plates draws: an int, 'all', or True (= all)."""
    if draw is True or (isinstance(draw, str) and draw.strip().lower() == 'all'):
        return n_plates
    if draw is None or draw is False:
        return 0
    try:
        n = int(draw)
    except (TypeError, ValueError):
        raise ValueError(f"draw must be a number of plates or 'all', not {draw!r}") from None
    if n < 0:
        raise ValueError(f"draw must be 0 or more, not {n}")
    return min(n, n_plates)


def check_dr_plates(info, draw=1):
    """Check each *_DR_plate.xlsx step 2 wrote against the user's input map.

    For every destination plate, the SMDC ID in each well of the DR file is
    compared with block 0 of that plate's tab in the input workbook. Numbers
    and ID-lot cells like 12346-1 are compounds; text cells (controls, DMSO,
    anything in REPLACE_DICT) are expected to be empty in the DR file.

    Args:
        info (DataFrame): what make_dr_platemaps_and_folders returned.
        draw (int or 'all'): how many DR files to draw (layout with controls,
            plus concentrations), to compare by eye with the input tab.
            1 (default) draws the first plate, 'all' draws every plate,
            0 or None draws none. True means 'all'.
    Returns:
        (summary, mismatches): summary has one row per plate (compound wells
        in the input, compound wells in the DR file, mismatched wells, wells
        with no Lot, no source plate or no concentration). mismatches lists
        every well where the two disagree. Both empty of problems = file matches.
    """
    dr_plates = info.attrs.get('dr_plates')
    infile = info.attrs.get('infile')
    if dr_plates is None or infile is None:
        raise ValueError("info has no record of its DR plate files. Pass the DataFrame that "
                         "make_dr_platemaps_and_folders returned (with dose_response=True).")
    n_draw = _plates_to_draw(draw, len(dr_plates))
    dp, _ = _read_info_sheet(infile)
    rows, bad = [], []
    for i, (grid_num, path) in enumerate(dr_plates):
        dr = pd.read_excel(path)
        want = read_platemap_block(plate_map=dp[grid_num[0:31]], block=0, value_name='Input')
        want['Destination well'] = want.row + want.col
        want['Expected'] = want.Input.map(_expected_smdc_id)
        got = dr.groupby('Destination well', as_index=False).SMDC_ID.first()
        both = want[['Destination well', 'Input', 'Expected']].merge(got, how='outer')
        differs = ~((both.Expected == both.SMDC_ID) | (both.Expected.isna() & both.SMDC_ID.isna()))
        mism = both[differs].assign(**{'Destination plate': grid_num})
        bad.append(mism[['Destination plate', 'Destination well', 'Input', 'Expected', 'SMDC_ID']])
        compounds = dr[dr.SMDC_ID.notna()]
        rows.append({
            'Destination plate': grid_num,
            'compound wells (input)': int(both.Expected.notna().sum()),
            'compound wells (DR file)': int(compounds['Destination well'].nunique()),
            'mismatched wells': int(differs.sum()),
            'no Lot': int(compounds.Lot.isna().sum()) if 'Lot' in dr else len(compounds),
            'no source plate': (int(compounds['Source plate'].isna().sum())
                                if 'Source plate' in dr else len(compounds)),
            'no conc': (int(compounds['[compound] uM'].isna().sum())
                        if '[compound] uM' in dr else len(compounds)),
            'positive control wells': int(dr['Sci_PosControl'].astype(bool).sum()) if 'Sci_PosControl' in dr else 0,
            'negative control wells': int(dr['Sci_NegControl'].astype(bool).sum()) if 'Sci_NegControl' in dr else 0,
            'file': path,
        })
        if i < n_draw:
            _draw_dr_layout(grid_num, dr)
            if '[compound] uM' in dr and dr['[compound] uM'].notna().any():
                # No colourbar: the printed concentrations are the key, and the
                # colour change between dilution steps is enough to spot a slip.
                draw_plate_heatmap(barcode=grid_num, df=dr, wells='Destination well',
                                   values='[compound] uM', labels='[compound] uM', plate_size=384,
                                   cbar=False)
    if 0 < n_draw < len(dr_plates):
        logger.warning(f"Drew {n_draw} of {len(dr_plates)} DR plates. Use draw='all' to draw every plate.")
    summary = pd.DataFrame(rows)
    mismatches = pd.concat(bad, ignore_index=True) if bad else pd.DataFrame()
    n = int(summary['mismatched wells'].sum()) if len(summary) else 0
    if n:
        logger.warning(f"{n} well(s) in the DR plate files don't match the input map. See mismatches.")
    else:
        logger.info(f"All {len(summary)} DR plate file(s) match the input map.")
    return summary, mismatches


# ###############################################################################
# inv_well_cols=['IW_SOURCE_ID','IW_IP_ID','IW_BARCODE','IW_WELL_INDEX','IW_COORDINATES','IW_ICL_ID','IW_VOL_UL','IW_CONC_UM','IW_FINAL_CONC_NM','IW_DILUTION_FACTOR','IW_FREEZETHAW','IW_INS_BY','IW_INS_DATE','IW_CONTAM','IW_NOTES']
# # def make_well_upload_list(pm):

###############################################################################
###############################################################################
# SECTION 7: post-ingest housekeeping
###############################################################################
def clean_up_folders(infile=None, basedir='./', delete_orig=False, share_root=None):
    """ After data is uploaded into HiTS, delete files from all_data (optional)
    and move each dated folder from toload to loaded.

    Run this only once the HiTS ingest has been confirmed. It is destructive and
    has no dry-run mode.

    Args:
        infile (path): path to the plate map info file. Must be the SAME file
            passed to make_dr_platemaps_and_folders, since the toload_<i> paths
            are recomputed from it by re-running the same enumerate logic.
        basedir (path): path to the directory for which to move data
        delete_orig (bool): whether to delete the original data from its
            location on the share. Irreversible; leave False unless the copy in
            loaded/ is confirmed good.
        share_root (str, optional): where the Montreal share is mounted; only
            used with delete_orig=True. Default: default_share_root().

    Folder locations are computed with the same helpers the pipeline uses
    (_read_info_sheet, _iter_load_groups, _assay_dir_name), so the two cannot
    disagree about where a folder lives. They used to, which is why clean-up
    could fail to find a folder that plainly existed.
    """
    # read data (shared with make_dr_platemaps_and_folders)
    infile = resolve_share_path(infile, share_root)
    basedir = resolve_share_path(basedir, share_root)
    dp, info = _read_info_sheet(infile)

    # Work out every move first and refuse to start if any would overwrite a
    # different file already in loaded/ (e.g. an earlier batch with the same
    # screen date). Identical files are fine: that is a rerun after an
    # interrupted copy.
    moves = _plan_loaded_moves(info, basedir)
    conflicts = _loaded_conflicts(moves)
    if conflicts:
        raise FileExistsError(
            "These files are already in loaded/ with different contents, so nothing "
            "was moved or deleted. Rename the folder in toload_<n> or move it by hand:\n  "
            + "\n  ".join(conflicts[:20])
        )

    # move each date from toload to loaded folders. Copy then remove, rather
    # than shutil.move, so an interrupted run leaves the source intact.
    for oldpath, newpath in moves:
        logger.info(f'{oldpath} -> {newpath}')
        os.makedirs(os.path.dirname(newpath), exist_ok=True)
        _copy_tree(oldpath, newpath, ignore=ignore_files)
        missing = _missing_from_copy(oldpath, newpath, ignore=ignore_files)
        if missing:
            raise OSError(f"Copy of {oldpath} to {newpath} is incomplete, so nothing was removed from "
                          "toload and no original was deleted. Missing or different: " + ", ".join(missing[:10]))
        shutil.rmtree(oldpath)  # raises on failure; the full copy is already in loaded/

    # Delete the originals only after every move above has finished, and only
    # when this original's own copy is in loaded/ (the same assay folder, date
    # and file name the pipeline copied it to). Any failure above raises
    # before this point, so the originals stay put. Contents aren't compared:
    # the tidy step may have edited the toload copy on purpose.
    if delete_orig:
        root = share_root or default_share_root()
        copies = {}
        for df, _rpt in _iter_load_groups(info):  # the rows the pipeline built folders from
            for _, r in df.iterrows():
                rel = _normalize_share_path(r.Raw_data_file)
                copies.setdefault(rel, set()).add(os.path.join(
                    basedir, _assay_dir_name(r.Raw_assay_ID, r.Assay_name),
                    'loaded', str(int(float(r.Date_of_screen))), rel.rsplit('/', 1)[-1]))
        for rel, expected in copies.items():
            # Same path normalisation the pipeline used to find the file.
            rdf = os.path.join(root, rel)
            if not os.path.exists(rdf):
                continue
            if not any(os.path.isfile(c) for c in expected):
                logger.warning(f"Not deleting {rdf}: its copy isn't in loaded/ (looked for {sorted(expected)})")
                continue
            logger.info(f"deleting original {rdf}")
            # os.remove rather than os.system("rm ..."): no shell involved,
            # so spaces and quotes in the path cannot make this delete
            # something else.
            os.remove(rdf)


CONTROL_COLUMNS = {'sciposcontrol': 'Sci_PosControl', 'scinegcontrol': 'Sci_NegControl'}


GROUP_COLUMNS = ['Raw_assay_ID', 'Assay_name', 'Date_of_screen', 'Raw_data_file']


def _plate_info_row(info, group_row, grid_num):
    """ The Info settings for one destination plate inside a load group.

    _iter_load_groups joins every non-key column across the plates that share a
    raw file, assay and date (e.g. Source_plate 'use map;use lots'). Folders
    and copies are shared, but the per-plate settings must not be, so take
    them from the plate's own Info row. Falls back to the group row if the
    plate cannot be found exactly once.
    """
    out = group_row.copy()
    out['Grid_num'] = grid_num
    keys = [c for c in GROUP_COLUMNS if c in info.columns and c in group_row.index]
    same_group = pd.Series(True, index=info.index)
    for c in keys:
        same_group &= info[c].astype(str) == str(group_row[c])
    has_plate = info['Grid_num'].astype(str).map(
        lambda v: grid_num in [x.strip() for x in v.split(';')])
    match = info[same_group & has_plate]
    if len(match) != 1:
        if len(match) > 1:
            logger.warning(f"{grid_num} is on {len(match)} Info rows; using the joined settings")
        return out
    own = match.iloc[0]
    for c in own.index:
        if c not in keys and c != 'Grid_num':
            out[c] = own[c]
    return out


def _apply_control_flags_from_map(dids, dp, row, control_dict):
    """ Re-read block 0 of the original plate map and set the control flags.

    Runs after all plate/lot/source logic, so the scientist's annotations are
    read before they were replaced during SMDC_ID processing. Labels match the
    same normalised way as in get_destination_plate_info. Biomek worklists are
    long format with no block 0, so their flags are left as they are.
    """
    if not control_dict:
        return dids
    if _source_plate_mode(_series_get_ci(row, 'Source_plate', '')) == 'biomek':
        logger.info(f"{row.Grid_num}: Biomek worklist, no plate-map annotations to read for controls")
        return dids
    dids = dids.copy()
    map_sheet = dp[str(row.Grid_num)[0:31]].dropna(how='all')
    annotations = read_platemap_block(plate_map=map_sheet, block=0, value_name='annotation')
    annotations['Destination well'] = annotations['row'] + annotations['col']
    annotation_keys = annotations['annotation'].map(_norm_key)
    for control_col, control_value in _normalize_control_dict(control_dict).items():
        control_wells = annotations.loc[annotation_keys == control_value, 'Destination well']
        dids[control_col] = dids['Destination well'].isin(control_wells)
    logger.info(
        f"Final control flags for {row.Grid_num}: "
        f"Sci_PosControl={int(dids['Sci_PosControl'].sum())}, "
        f"Sci_NegControl={int(dids['Sci_NegControl'].sum())}"
    )
    return dids


def _normalize_control_dict(control_dict):
    """ {'Sci_PosControl': 'No Enzyme'} -> {'Sci_PosControl': 'noenzyme'}.

    Output keys are the canonical column names whatever their spelling in the
    notebook; labels are lowercased alphanumerics for comparison only.
    """
    out = {}
    for key, label in (control_dict or {}).items():
        if label is None:
            continue
        out[CONTROL_COLUMNS.get(_norm_key(key), key)] = _norm_key(label)
    return out


def _id_plate_pairs(dids):
    """ Unique (SMDC_ID, Source plate) pairs taken from the same rows.

    lookup_lot_by_plate_and_ic_id(matched=True) zips ids with plates, so both
    lists must come from the same wells in the same order.
    """
    if 'Source plate' not in dids.columns:
        return [], []
    plates = dids['Source plate'].astype(str).str.strip()
    ok = dids.SMDC_ID.notna() & dids['Source plate'].notna() & ~plates.str.upper().isin(['', 'NA', 'NAN'])
    pairs = list(dict.fromkeys(zip(dids.loc[ok, 'SMDC_ID'].astype(float).astype(int).astype(str), plates[ok])))
    return [i for i, _ in pairs], [p for _, p in pairs]


def _ensure_control_flags(dids):
    """ Sci_PosControl / Sci_NegControl as booleans, False where missing.

    Biomek worklists skip get_destination_plate_info, which is what normally
    adds these columns.
    """
    dids = dids.copy()
    for col in ['Sci_PosControl', 'Sci_NegControl']:
        if col not in dids.columns:
            dids[col] = False
        dids[col] = dids[col].astype(object).where(dids[col].notna(), False).astype(bool)
    return dids


def _unique_user_lots(lots):
    """ SMDC_ID -> Lot from lookup_lot_by_user, refusing to guess between lots.

    Two users can each have registered a lot of the same compound. The same lot
    twice is fine. Different lots mean this last-resort fallback cannot know
    which one was pipetted, so it raises instead of picking one.
    """
    lots = lots.drop_duplicates()
    conflicts = lots[lots.SMDC_ID.duplicated(keep=False)]
    if not conflicts.empty:
        detail = conflicts.groupby('SMDC_ID')['Lot'].apply(lambda x: sorted(x.tolist())).to_dict()
        raise ValueError(
            "More than one user registered a lot for these compounds, so the lot "
            f"cannot be guessed: {detail}. Put the lot in the Info sheet "
            "(Lot_info) or pass only the user who made this plate."
        )
    return lots.set_index("SMDC_ID")["Lot"]


def _append_grid_num(outdir, gridnum):
    """ Create or extend grid_num.txt, the manifest the HiTS loader reads to
    know which destination plates this folder's data covers.

    Appended to rather than overwritten, since a folder may be filled across
    several runs. Read as text so an all-digit barcode keeps its leading zeros.
    """
    path = os.path.join(outdir, 'grid_num.txt')
    gridnum = gridnum.astype(str)
    if os.path.exists(path):
        gn = pd.read_csv(path, dtype=str, keep_default_na=False)
        gn = pd.concat([gn, gridnum]).drop_duplicates(keep='first')
    else:
        gn = gridnum
    gn['grid_num'] = [str(x)[0:50] for x in gn['grid_num']]  # HiTS barcode field caps at 50
    with _FileInUse(path):
        gn.to_csv(path, index=False)
    return gn


def _plan_loaded_moves(info, basedir):
    """ (toload_<i>/<datadir>, loaded/<datadir>) pairs for clean_up_folders. """
    moves = []
    for df, rpt in _iter_load_groups(info):
        for i, row in df.iterrows():
            if rpt:
                i = 0
            outdir = os.path.join(basedir, _assay_dir_name(row.Raw_assay_ID, row.Assay_name))
            toload = os.path.join(outdir, f'toload_{i}')
            loaded = os.path.join(outdir, 'loaded')
            if not os.path.exists(toload):
                logger.warning(f'No such folder, skipping: {toload}')
                continue
            for datadir in sorted(os.listdir(toload)):
                # skip macOS metadata that is not a real directory
                if datadir == ".DS_Store" or datadir.startswith('._'):
                    continue
                pair = (os.path.join(toload, datadir), os.path.join(loaded, datadir))
                if pair not in moves:
                    moves.append(pair)
    return moves


def _loaded_conflicts(moves):
    """ Files a move would overwrite with different contents. """
    conflicts = []
    for oldpath, newpath in moves:
        if not os.path.isdir(oldpath):
            continue
        for root, _, files in os.walk(oldpath):
            for f in files:
                if f == '.DS_Store' or f.startswith('._'):
                    continue
                src = os.path.join(root, f)
                dst = os.path.join(newpath, os.path.relpath(src, oldpath))
                if os.path.exists(dst) and not filecmp.cmp(src, dst, shallow=False):
                    conflicts.append(dst)
    return conflicts


def ignore_files(dir, files):
    """copytree ignore callback: skip .DS_Store and macOS ._ resource forks.

    Args:
        dir: directory being visited (required by the copytree API, unused).
        files (list): names in that directory.
    Returns:
        list: names for copytree to skip.
    """
    return [f for f in files if f == '.DS_Store' or f.startswith('._')]

###############################################################################
# SECTION 7: fixing up *_DR_plate.xlsx files before upload
#
# Used by prepare_hits_upload.ipynb steps 4 and 5. Everything that edits files
# is a dry run unless write=True.
###############################################################################
def dr_plate_files(outdir, assay_ids):
    """ Yield every *_DR_plate.xlsx under outdir for the given assay IDs.

    Folders are matched by name prefix (`<assay id>_<assay name>`). Anything
    under a `loaded` folder and macOS `._` sidecar files are skipped.

    Args:
        outdir (str): the upload folder (the notebook's OUTDIR).
        assay_ids (list of str): assay IDs, e.g. ['2863', '2865'].
    """
    outdir = resolve_share_path(outdir)
    prefixes = tuple(f"{a}_" for a in assay_ids)
    exact = {str(a) for a in assay_ids}
    for d in sorted(os.listdir(outdir)):
        # '2863' must not pick up '28630_<name>', so match on '<id>_' or the bare id.
        if not (d.startswith(prefixes) or d in exact):
            continue
        for root, _, files in os.walk(os.path.join(outdir, d)):
            if 'loaded' in root:
                continue
            for f in sorted(files):
                if f.endswith('DR_plate.xlsx') and not f.startswith('._'):
                    yield os.path.join(root, f)


def find_duplicate_source_plates(outdir, assay_ids):
    """ List compounds that appear with more than one source plate/well.

    Args:
        outdir (str): the upload folder.
        assay_ids (list of str): assay IDs to check.
    Returns:
        DataFrame with columns Plate (the DR file), SMDC_ID and SPlateWell, a
        list of "<source plate> <source well>" entries for that compound.
    """
    frames = []
    for path in dr_plate_files(outdir, assay_ids):
        t = pd.read_excel(path, engine='openpyxl', keep_default_na=False)
        t = t[['SMDC_ID', 'Source plate', 'Source well']].drop_duplicates()
        t = t[~t.SMDC_ID.astype(str).isin(['', 'NA'])]  # empty and control wells are not compounds
        t = t[t.SMDC_ID.duplicated(keep=False)].copy()
        t['Plate'] = path
        frames.append(t)
    cols = ['Plate', 'SMDC_ID', 'SPlateWell']
    if not frames or all(f.empty for f in frames):
        return pd.DataFrame(columns=cols)
    dups = pd.concat(frames)
    dups['SPlateWell'] = dups['Source plate'].astype(str) + ' ' + dups['Source well'].astype(str)
    dups = dups.pivot_table(index=['Plate', 'SMDC_ID'], values='SPlateWell', aggfunc=list)
    return dups.sort_values('Plate').reset_index()[cols]


def remove_duplicate_source_plates(dups, keep=None, write=False):
    """ Keep one source plate/well per duplicated compound and drop the others.

    Args:
        dups (DataFrame): output of find_duplicate_source_plates.
        keep (list of int): for each row of dups, in order, the position of the
            RIGHT entry in its SPlateWell list. Defaults to 0 for every row.
        write (bool): False (default) only reports; True saves the files.
    Returns:
        dict: DR file path -> (rows before, rows after).
    """
    dups = dups.reset_index(drop=True)
    keep = [0] * len(dups) if keep is None else list(keep)
    if len(keep) != len(dups):
        raise ValueError(f"keep has {len(keep)} entries but dups has {len(dups)} rows")
    # Check every index before touching any file: an out-of-range index would
    # otherwise keep nothing and delete every row for that compound.
    for i, (k, entries) in enumerate(zip(keep, dups.SPlateWell)):
        if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or not 0 <= k < len(entries):
            raise ValueError(f"keep[{i}] = {k!r} but row {i} has {len(entries)} entries (use 0 to {len(entries) - 1})")
    changes = {}
    for plate in dups.Plate.unique():
        # keep_default_na=False keeps the loader's literal 'NA' cells as 'NA'
        # (and blanks blank) when the file is saved again.
        t = pd.read_excel(plate, engine='openpyxl', keep_default_na=False)
        before = len(t)
        for i, row in dups[dups.Plate == plate].iterrows():
            for j, entry in enumerate(row.SPlateWell):
                if j == keep[i]:
                    continue
                src_plate, src_well = entry.rsplit(' ', 1)
                # Compare as text: SPlateWell entries are strings, but a numeric
                # barcode is read from Excel as a number.
                t = t[~((t.SMDC_ID.astype(str) == str(row.SMDC_ID))
                        & (t['Source plate'].astype(str) == src_plate)
                        & (t['Source well'].astype(str) == src_well))]
        changes[plate] = (before, len(t))
        logger.info(f"{plate}: {before} -> {len(t)} rows")
        if write:
            with _FileInUse(plate):
                t.to_excel(plate, index=False)
    if not write:
        logger.info("Dry run: pass write=True to save.")
    return changes


def fill_control_lots(outdir, assay_ids, control_lots, write=False):
    """ Set the Lot of control compounds in every DR file for these assays.

    Empty cells are written as the literal 'NA' the loader expects.

    Args:
        outdir (str): the upload folder.
        assay_ids (list of str): assay IDs to fix.
        control_lots (dict): {SMDC_ID: lot}.
        write (bool): False (default) only reports; True saves the files.
    Returns:
        list of str: the DR files that were (or would be) changed.
    """
    paths = []
    for path in dr_plate_files(outdir, assay_ids):
        t = pd.read_excel(path, engine='openpyxl')
        for smdc, lot in control_lots.items():
            t.loc[t.SMDC_ID == smdc, 'Lot'] = lot
        t = t.fillna('NA')
        paths.append(path)
        logger.info(path)
        if write:
            with _FileInUse(path):
                t.to_excel(path, index=False)
    if not write:
        logger.info("Dry run: pass write=True to save.")
    return paths



def source_plate_barcodes(info):
    """ Real plate barcodes listed in an Info tab's Source_plate column.

    Splits on ';' and drops the dispatch keywords ("use lots", "use map",
    "use echo/biomek", "n.a."). Order is preserved, duplicates removed.
    """
    out = []
    if 'Source_plate' not in info.columns:
        return out
    for v in info['Source_plate'].dropna().astype(str):
        if _source_plate_mode(v) != 'barcodes':
            continue
        for b in v.replace('\n', ';').split(';'):
            b = b.strip()
            if b and b not in out:
                out.append(b)
    return out
