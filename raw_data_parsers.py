#!/usr/bin/env python
# 20260210 Start replicating Pipeline Pilot data parsing functions in Python. This is a work in progress, and will be updated as I go through the PP protocol and figure out what functions I need to replicate.

import pandas as pd
import os
import numpy as np
import logging
import re

logger = logging.getLogger(__name__)

plate_dict = {'384': ['P',16,24], '96': ['H',8,12], '1536': ['AF',32,48]}


def _natural_sort_key(value):
    """Sort filenames in human order (for example, plate2 before plate10)."""
    return [(0, int(part)) if part.isdigit() else (1, part.casefold())
            for part in re.split(r'(\d+)', str(value))]

######### GENERAL FILE PARSING FUNCTIONS #######


def load_csv_files(load_dir="./toload"):
    """
    Reads in all CSV files from a to-load folder, and returns a list of dataframes and corresponding filenames.
    Args:
        load_dir: Directory containing the CSV files to load
    Returns:
        dfs: List of dataframes containing the data from each CSV file
        filenames: List of filenames corresponding to each dataframe
    """
    load_dir = _existing_dir(load_dir)
    filenames=[]
    dfs=[]
    for root, dirs, files in os.walk(load_dir):
        for file in files:
            if file.endswith(".csv"):
                logger.info(f"Loading file: {file}")
                filenames.append(file)
                dfs.append(pd.read_csv(os.path.join(root, file), header=None))
    return dfs, filenames


def make_wells(df, row, column):
    """ fix column ids for well plate
    Args:
        df: dataframe of data blocks
    Returns:
        df: dataframe of data blocks
    """ 
    df[column]=df[column].astype(int)
    df['well']=df[row]+df[column].astype(int).astype(str).str.zfill(2)
    return df


def make_well_type(df, ctrls, neg_ctrl_locs=[1,2], pos_ctrl_locs=[23,24]):
    """ mark control wells in dataframe of data blocks
    Args:
        df: dataframe of data blocks
        ctrls: type of controls to mark (cols or rows)
        neg_ctrl_locs: list of column or row numbers corresponding to negative controls
        pos_ctrl_locs: list of column or row numbers corresponding to positive controls
    Returns:
        df: dataframe of data blocks with control wells marked
    """ 
    df=df.copy()
    if ctrls=='cols':
        df['Well_Type']=np.where(df['col'].isin(neg_ctrl_locs), 'negativeControl', "compound")
        df['Well_Type']=np.where(df['col'].isin(pos_ctrl_locs), 'positiveControl', df['Well_Type'])
    elif ctrls=='rows':
        df['Well_Type']=np.where(df['row'].isin(neg_ctrl_locs), 'negativeControl', "compound")
        df['Well_Type']=np.where(df['row'].isin(pos_ctrl_locs), 'positiveControl', df['Well_Type'])
    elif ctrls=='wells':
        df['Well_Type']=np.where(df['well'].isin(neg_ctrl_locs), 'negativeControl', "compound")
        df['Well_Type']=np.where(df['well'].isin(pos_ctrl_locs), 'positiveControl', df['Well_Type'])  
    return df


def prepare_data_for_upload(df,ctrls='cols', neg_ctrl_locs=[1,2], pos_ctrl_locs=[23,24],meas_dict={"Results for Jez US LUM 384 (cps)(1) - channel 0 (CPS)":"intensity"}):
    """ prepare dataframe of data blocks for upload by marking control wells and reordering columns
    Args:
        df: dataframe of data blocks
        ctrls: type of controls to mark ('cols', 'rows', or 'wells')
        neg_ctrl_locs: list of column or row numbers corresponding to negative controls
        pos_ctrl_locs: list of column or row numbers corresponding to positive controls
    Returns:
        df: dataframe of data blocks with control wells marked and columns reordered for upload
    """
    # add Well_Type column based on control locations
    df=make_well_type(df, ctrls, neg_ctrl_locs, pos_ctrl_locs)

    # reorder and rename columns for upload with specified assay result type names
    df.columns=[x if x not in meas_dict.keys() else meas_dict[x] for x in df.columns]

    # make sure measurement columns are numeric for upload
    try:
        df.loc[:,list(meas_dict.values())]=df.loc[:,list(meas_dict.values())].astype(float)
    except:
        logger.warning("Column names were not successfully translated")
    return df


####### ENVISION FILE PARSING FUNCTIONS #######


def _existing_dir(load_dir):
    """load_dir as a folder that exists, or a FileNotFoundError naming it.

    os.walk on a missing folder yields nothing, so a wrong path used to look
    like "0 files, nothing to fix". A Z:\\ or /Volumes/Shared/ path is mapped
    to where the share is mounted here (plate_maps.resolve_share_path) when
    plate_maps can be imported, so a path copied on Windows works in WSL.
    """
    path = str(load_dir)
    try:
        from plate_maps import resolve_share_path
    except ImportError:
        try:
            from database_interactions.plate_maps import resolve_share_path
        except ImportError:
            resolve_share_path = None
    if resolve_share_path is not None:
        path = resolve_share_path(path)
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Folder not found: {load_dir!r}"
                                + (f" (looked for {path!r})" if path != str(load_dir) else "")
                                + ". Check the path, and that the share is mounted.")
    return path


def parse_envision_files(plate_type='384', load_dir="./toload"):
    """ main function to parse envision files. Reads in all CSV files from a to-load folder, extracts plate and background info, and extracts data blocks from each file.
    Args:        load_dir: Directory containing the CSV files to load
    Returns:         dats: List of dataframes containing the data from each CSV file, with plate and background info extracted and associated with each data block
    """
    dfs, filenames = _envision_load_csv_files(load_dir, plate_type=plate_type)

    dats=[]
    for filename, df in zip(filenames, dfs):
        # get rid of basic assay info section, but keep for warnings checks
        assay_info_idx=df[df[0].astype(str).str.startswith("Basic assay information")].index
        if len(assay_info_idx) == 0:
            raise _envision_format_error(filename, "missing the Basic assay information section")
        assay_info=df.loc[assay_info_idx[0]:]
        df=df.loc[:assay_info_idx[0]-1]

        # get indexes of plate, bg, and result sections
        plate_idxs=df[df[0].astype(str).str.startswith("Plate information")].index
        bg_idxs=df[df[0].astype(str).str.startswith("Background information")].index
        res_idxs=df[df[0].astype(str).str.upper().str.contains("RESULT")].index

        _validate_envision_section_order(df, filename, plate_idxs, bg_idxs, res_idxs)
        
        # extract plate and bg rows
        plate_info=pd.concat([_envision_extract_info_block(df, idx+1, idx+3) for idx in plate_idxs]).reset_index(drop=True)
        bg_info=pd.concat([_envision_extract_info_block(df, idx+1, idx2+1) for idx,idx2 in zip(bg_idxs,res_idxs)]).dropna(how='all').reset_index(drop=True)

        # extract data index and data name
        dat_idxs=df[df[df.columns[0]].astype(str).str.upper().str.contains('RESULTS')].reset_index(names='idx')
        dat_idxs=dat_idxs.iloc[:,0:2]
        # extract all data
        dat=_envision_extract_all_data(df, dat_idxs, filename, plate_info, plate_type)
        dat=make_wells(dat, 'row','col')
        dat=dat.sort_values(['plate','row','col']).reset_index(drop=True)
        dats.append(dat)

    dat=pd.concat(dats, ignore_index=True)
    return dat


def _envision_format_error(filename, detail):
    return ValueError(
        f"Unsupported EnVision CSV format in {filename!r}: {detail}. "
        "Re-export the data with Plate information and Background information "
        "enabled and placed at the beginning of each plate, before measurement data."
    )


def _validate_envision_section_order(df, filename, plate_idxs, bg_idxs, res_idxs):
    """Reject exports whose plate/background metadata cannot be paired safely."""
    if len(res_idxs) == 0:
        raise _envision_format_error(filename, "no measurement data sections were found")
    previous_result_idx = -1
    for number, result_idx in enumerate(res_idxs, start=1):
        plate_sections = plate_idxs[(plate_idxs > previous_result_idx) & (plate_idxs < result_idx)]
        background_sections = bg_idxs[
            (bg_idxs > previous_result_idx) & (bg_idxs < result_idx)
        ]
        if len(plate_sections) != 1 or len(background_sections) != 1:
            raise _envision_format_error(
                filename,
                f"measurement section {number} must have exactly one Plate information "
                "and one Background information section after the prior measurement "
                "and before its data; found "
                f"{len(plate_sections)} plate and {len(background_sections)} background sections",
            )
        headers = df.iloc[plate_sections[0] + 1, :].astype(str).str.strip().tolist()
        if not {"Plate", "Repeat", "Barcode"}.issubset(headers):
            raise _envision_format_error(
                filename,
                f"Plate information section {number} is missing Plate, Repeat, or Barcode fields",
            )
        previous_result_idx = result_idx

def _envision_load_csv_files(load_dir="./toload", plate_type='384'):
    """
    Reads in all envision CSV files from a to-load folder, and returns a list of dataframes and corresponding filenames.
    Assumes that Plate Info and Background Info sections are present and exported at the beginning of the data block.
    Args:
        load_dir: Directory containing the CSV files to load
    Returns:
        dfs: List of dataframes containing the data from each CSV file
        filenames: List of filenames corresponding to each dataframe
    """
    load_dir = _existing_dir(load_dir)
    filenames=[]
    dfs=[]
    for root, dirs, files in os.walk(load_dir):
        dirs[:] = sorted(
            (directory for directory in dirs if directory.casefold() != 'loaded'),
            key=_natural_sort_key,
        )
        for file in sorted(files, key=_natural_sort_key):
            if file.endswith(".csv") and not file.startswith("._"): # get rid of hidden files that sometimes come with Mac CSV exports
                logger.info(f"Loading file: {file}")
                filenames.append(file)
                try:
                    df=pd.read_csv(os.path.join(root, file), header=None)
                except Exception as e:
                    logger.error(f"Error reading {file}: {e}. Falling back to envision-specific CSV loading for {file}")
                    ncols=plate_dict[plate_type][2]+2
                    df=pd.read_csv(os.path.join(root, file),  header=None, names=np.arange(0,ncols), usecols=np.arange(0,ncols), skip_blank_lines=False).dropna(how='all', axis=1)
                    
                dfs.append(df)
    return dfs, filenames

def _envision_extract_info_block(df, idx, idx2):
    """ extract info associated with each data block in an envision file.
    Args:
        df: dataframe of full envision file
        idx: start line index of info block
        idx2: end line index of info block
    Returns:
        info: dataframe of info block
    """    
    info=df.iloc[idx:idx2,:]
    info=info.T.reset_index(drop=True)
    info=info.dropna(how='all')
    info=info.set_index(info.columns[0]).T.reset_index(names='idx')
    return info


def _envision_extract_data_block(df, dat_idxs, i, idx2, filename, plate_info):
    """ extract data from a single matrix location in an envision file
    Args:
        df: dataframe of full envision file
        idx: start line index of data block
        idx2: end line index of data block
    Returns:
        dat: dataframe of data block
    """    
    dat=df.iloc[dat_idxs.idx.iloc[i]+1:dat_idxs.idx.iloc[i]+idx2+1+1,:].copy()
    dat.iloc[0,0]='col'
    dat=dat.rename(columns={df.columns[0]:'row'}).set_index('row').T.set_index('col').T.melt(ignore_index=False).reset_index()
    dat['plate']=f"{filename}_{plate_info.Plate.iloc[i]}_{plate_info.Repeat.iloc[i]}_{plate_info.Barcode.iloc[i]}"
    dat['meas']=dat_idxs[0].iloc[i]
    return dat


def _envision_extract_all_data(df, dat_idxs, filename, plate_info, plate_type):
    """ extract data from a list of matrix locations in an envision file
    Args:
        df: dataframe of full envision file
        dat_idxs: dataframe of locations of data blocks and measurements from extract_info()
        filename: name of file
        plate_info: dataframe of plate info
    Returns:
        dat: dataframe of data blocks
    """ 
    idx2=plate_dict[plate_type][1]
    dats=[]
    for i in range(len(dat_idxs)):
        dat=_envision_extract_data_block(df, dat_idxs, i, idx2, filename, plate_info)
        dats.append(dat)
    dat=pd.concat(dats)
    dat=dat.pivot_table(index=['plate','row','col'], columns='meas', values='value', aggfunc=lambda x: x.iloc[0] if len(x) == 1 else Exception()).reset_index()
    # pivot_table drops wells whose measurements are all blank. PP's plate
    # reader still emits those wells, so restore the full plate geometry.
    row_count, column_count = plate_dict[plate_type][1:]
    row_labels = []
    for number in range(1, row_count + 1):
        label = ''
        while number:
            number, remainder = divmod(number - 1, 26)
            label = chr(65 + remainder) + label
        row_labels.append(label)
    wells = pd.MultiIndex.from_product(
        [row_labels, range(1, column_count + 1)], names=['row', 'col']
    )
    dat['col'] = pd.to_numeric(dat['col'], errors='coerce')
    restored = []
    for plate in dat['plate'].drop_duplicates():
        plate_data = dat.loc[dat['plate'] == plate].drop(columns='plate').set_index(['row', 'col'])
        plate_data = plate_data.reindex(wells).reset_index()
        plate_data.insert(0, 'plate', plate)
        restored.append(plate_data)
    return pd.concat(restored, ignore_index=True)



######### ENVISION CSV MAINTENANCE UTILITIES #########

import re
import csv
import io
from typing import Tuple

def _detect_newline_style_from_bytes(b: bytes) -> str:
    """Detect dominant newline style in a byte string.
    Returns one of "\r\n", "\n", or "\r". Defaults to "\n" if unknown."""
    if b.count(b"\r\n") > 0:
        return "\r\n"
    # If file contains CR without LF, treat as classic Mac
    if b.count(b"\r") > 0 and b.count(b"\n") == 0:
        return "\r"
    return "\n"


def _is_csv_blank_line(s: str) -> bool:
    """Return True if the line is visually blank or contains only empty CSV fields.

    A CSV line counts as blank if all tokens separated by commas are either:
      - empty (after trimming whitespace), or
      - empty quoted strings like "" (optionally with internal whitespace only)
    """
    # Fast path: if no comma present, just check if empty or empty quoted
    t = s.strip()
    if ',' not in s:
        if t == '':
            return True
        if (t.startswith('"') and t.endswith('"') and t[1:-1].strip() == ''):
            return True
        return False

    # Split by commas and validate each token
    for token in s.split(','):
        tt = token.strip()
        if tt == '':
            continue
        if tt.startswith('"') and tt.endswith('"') and tt[1:-1].strip() == '':
            continue
        return False
    return True


def ensure_calc_spacing_in_file(csv_path: str, required_empty_lines: int = 3) -> Tuple[bool, int]:
    """Ensure there are N empty lines immediately above any line that starts a
    calculation block labeled like "Calc X: Z'" where X is an integer.

    A line above counts as "empty" if it is either a true empty line or a line
    consisting only of commas (and optional whitespace) representing empty CSV fields.

    The function edits the file in place only if changes are needed.

    Args:
        csv_path: Path to the Envision CSV file.
        required_empty_lines: Number of blank lines that must precede the label line (default 3).

    Returns:
        (changed, edits):
            changed (bool): True if the file was modified.
            edits (int): Number of blank lines inserted in this file.
    """

    # Read bytes to detect newline style robustly
    try:
        with open(csv_path, 'rb') as f:
            raw = f.read()
    except FileNotFoundError:
        logger.error(f"File not found: {csv_path}")
        return (False, 0)
    except Exception as e:
        logger.error(f"Error reading file {csv_path}: {e}")
        return (False, 0)

    newline_style = _detect_newline_style_from_bytes(raw)
    # Decode text (assume UTF-8, fall back to latin-1 to avoid data loss on odd encodings)
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        text = raw.decode('latin-1')

    # Split into logical lines (no newline tokens)
    lines = text.splitlines()

    # Pattern matches a calc block label anywhere in the row (after optional empty leading cells),
    # accepts optional "Calculated results:" prefix, optional surrounding quotes for the cell,
    # and handles ASCII or Unicode right single quote after Z. Allows trailing columns.
    # Use a permissive search: look for the sequence anywhere in the row since the cell may continue after Z'
    # Two patterns: prefer rows with "Calculated results:", fallback to any "Calc X: Z'"
    pattern_calc_results = re.compile(r"Calculated\s*results\s*:.*?Calc\s+\d+\s*:\s*Z['\u2019]", re.IGNORECASE)
    pattern_fallback = re.compile(r"\bCalc\s+\d+\s*:\s*Z['\u2019]", re.IGNORECASE)

    i = 0
    inserted = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        # Match only if line contains Calculated results: ... Calc X: Z',
        # otherwise fall back to any Calc X: Z' lines that appear at top-of-section.
        if (pattern_calc_results.search(line) is not None) or \
           (pattern_fallback.search(line) is not None and (i+1 < n and lines[i+1].strip()=='')):
            # Count how many empty (or comma-only) lines are immediately above
            k = 0
            j = i - 1
            # Track style of the nearest blank: commas count if comma-only is found
            nearest_blank_commas = None
            while j >= 0 and _is_csv_blank_line(lines[j]):
                if nearest_blank_commas is None and lines[j].strip() != '':
                    # comma-only blank line -> record how many commas to replicate style
                    nearest_blank_commas = lines[j].count(',')
                k += 1
                j -= 1
            if k < required_empty_lines:
                blanks_needed = required_empty_lines - k
                insertion_index = i - k  # insert before existing blank block (or right above calc)
                # Choose insertion style:
                # - If a comma-only blank was seen just above, mirror its comma count
                # - Else insert truly empty lines (safer for parsers that skip blank lines)
                if nearest_blank_commas is not None and nearest_blank_commas > 0:
                    insert_str = ',' * nearest_blank_commas
                else:
                    insert_str = ''
                lines[insertion_index:insertion_index] = [insert_str] * blanks_needed
                inserted += blanks_needed
                i += blanks_needed  # move past inserted lines
                n += blanks_needed
        i += 1

    if inserted > 0:
        # Re-join with original newline style; append a trailing newline for POSIX-friendliness
        new_text = (newline_style.join(lines)) + newline_style
        try:
            with open(csv_path, 'w', newline='') as f:
                f.write(new_text)
            logger.info(f"Updated {csv_path}: inserted {inserted} blank line(s)")
            return (True, inserted)
        except Exception as e:
            logger.error(f"Failed to write updates to {csv_path}: {e}")
            return (False, 0)
    else:
        logger.debug(f"No changes needed for {csv_path}")
        return (False, 0)


def ensure_calc_spacing_in_dir(load_dir: str, required_empty_lines: int = 3) -> Tuple[int, int, int]:
    """Walk a directory and ensure calc spacing for all CSVs (skips '._' files).

    Args:
        load_dir: Directory to search for CSV files.
        required_empty_lines: Number of blank lines to enforce above each calc label.

    Returns:
        (files_seen, files_changed, total_inserted)
    """
    load_dir = _existing_dir(load_dir)
    files_seen = 0
    files_changed = 0
    total_inserted = 0

    for root, dirs, files in os.walk(load_dir):
        for file in files:
            if not file.lower().endswith('.csv'):
                continue
            if file.startswith('._'):
                continue
            csv_path = os.path.join(root, file)
            files_seen += 1
            changed, inserted = ensure_calc_spacing_in_file(csv_path, required_empty_lines)
            if changed:
                files_changed += 1
                total_inserted += inserted

    if files_seen == 0:
        logger.warning(f"No CSV files found under '{load_dir}'. Is this the right folder?")
    elif files_changed == 0:
        logger.info(f"Found {files_seen} CSV file(s) in '{load_dir}': spacing already OK, nothing changed")
    else:
        logger.info(f"Found {files_seen} CSV file(s) in '{load_dir}': fixed spacing in {files_changed}, "
                    f"inserted {total_inserted} blank line(s)")
    return files_seen, files_changed, total_inserted


def fill_zprime_blocks_in_file(csv_path: str, require_prefix: bool = True) -> Tuple[bool, int]:
    """Fill Z' data blocks so the full 16x24 matrix (A1..P24) contains the Z' value from A1.

    After a line that begins with "Calculated results:" and contains "Calc X: Z'",
    this function searches ahead for the rigid header row with an empty first cell
    and columns labeled 1..24 in columns 2..25, followed by 16 data rows labeled
    A..P in the first column.

    It writes the A1 value to all cells in columns 1..24 for rows A..P (i.e., the
    entire matrix), regardless of what is already present.

    Args:
        csv_path: Path to the Envision CSV file.
        require_prefix: If True, only consider blocks with the explicit
                        "Calculated results:" prefix; if False, also consider
                        generic "Calc X: Z'" lines.

    Returns:
        (changed, blocks_filled): whether the file was modified and how many Z' blocks were filled.
    """
    try:
        with open(csv_path, 'rb') as f:
            raw = f.read()
    except FileNotFoundError:
        logger.error(f"File not found: {csv_path}")
        return (False, 0)
    except Exception as e:
        logger.error(f"Error reading file {csv_path}: {e}")
        return (False, 0)

    newline_style = _detect_newline_style_from_bytes(raw)
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        text = raw.decode('latin-1')

    lines = text.splitlines()

    # Patterns to locate the Z' block label rows
    pat_calc_results = re.compile(r"Calculated\s*results\s*:.*?Calc\s+\d+\s*:\s*Z['\u2019\u2032]", re.IGNORECASE)
    pat_fallback = re.compile(r"\bCalc\s+\d+\s*:\s*Z['\u2019\u2032]", re.IGNORECASE)

    def parse_csv_row(s: str):
        return next(csv.reader([s]))

    def to_csv_row(cells: list) -> str:
        buf = io.StringIO()
        w = csv.writer(buf, lineterminator='')
        w.writerow(cells)
        return buf.getvalue()

    def is_rigid_header(cells: list) -> bool:
        """Return True if header has empty first cell then columns 1..24 in positions 2..25.
        Accepts leading zeros (e.g., '01'..'09') as valid labels.
        """
        if len(cells) < 25:
            return False
        if str(cells[0]).strip() != '':
            return False
        for i in range(1, 25):
            v = str(cells[i]).strip()
            if v == '':
                return False
            try:
                if int(v) != i:
                    return False
            except Exception:
                return False
        return True

    blocks_filled = 0
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]
        if (require_prefix and pat_calc_results.search(line)) or (not require_prefix and (pat_calc_results.search(line) or pat_fallback.search(line))):
            # Advance to next non-blank line for the header
            j = i + 1
            while j < n and _is_csv_blank_line(lines[j]):
                j += 1
            if j >= n:
                i += 1
                continue
            try:
                header_cells = parse_csv_row(lines[j])
            except Exception:
                i += 1
                continue
            # Rigid header: empty first cell then 1..24 in columns 2..25 (accepts '01' etc.)
            if not is_rigid_header(header_cells):
                i += 1
                continue

            # Expect 16 rows A..P starting at j+1
            data_start = j + 1
            data_end = data_start + 16
            if data_end > n:
                i += 1
                continue

            # Parse rows A..P
            rows = []
            valid_block = True
            for r in range(data_start, data_end):
                try:
                    cells = parse_csv_row(lines[r])
                except Exception:
                    valid_block = False
                    break
                rows.append(cells)
            if not valid_block:
                i += 1
                continue

            # Validate labels A..P in first column
            expected_rows = [chr(ord('A') + k) for k in range(16)]
            for k, cells in enumerate(rows):
                if len(cells) == 0 or str(cells[0]).strip().upper() != expected_rows[k]:
                    valid_block = False
                    break
            if not valid_block:
                i += 1
                continue

            # Get A1 value (row A, column index 1)
            # Ensure cell width for all rows preserves original trailing length
            rowA = rows[0]
            # Determine the max width among these rows to preserve trailing empties
            max_len = max(len(c) for c in rows)
            if max_len < 25:
                max_len = 25
            # Normalize widths and capture original lengths
            orig_lengths = [len(c) for c in rows]
            rows_norm = []
            for c in rows:
                cc = c[:] + [''] * (max_len - len(c))
                rows_norm.append(cc)

            a1_value = rows_norm[0][1]

            # Fill columns 3..22 inclusive with A1 value for all rows A..P
            for r in range(16):
                for col_idx in range(3, 23):
                    rows_norm[r][col_idx] = a1_value

            # Rebuild the 16 data lines using original lengths where possible
            new_data_lines = []
            for r in range(16):
                target_len = max(orig_lengths[r], 25)
                cells_out = rows_norm[r][:target_len]
                new_data_lines.append(to_csv_row(cells_out))

            # Write back into lines
            for offset, new_line in enumerate(new_data_lines):
                lines[data_start + offset] = new_line

            blocks_filled += 1
            # Move i past this block to avoid re-matching within the same block
            i = data_end
            continue
        i += 1

    if blocks_filled > 0:
        new_text = (newline_style.join(lines)) + newline_style
        try:
            with open(csv_path, 'w', newline='') as f:
                f.write(new_text)
            logger.info(f"Updated {csv_path}: filled {blocks_filled} Z' block(s)")
            return (True, blocks_filled)
        except Exception as e:
            logger.error(f"Failed to write updates to {csv_path}: {e}")
            return (False, 0)
    else:
        logger.debug(f"No Z' blocks found/updated for {csv_path}")
        return (False, 0)


def fill_zprime_blocks_in_dir(load_dir: str, require_prefix: bool = True) -> Tuple[int, int]:
    """Fill Z' blocks across all CSVs in a directory tree.

    Args:
        load_dir: Directory to search for CSV files.
        require_prefix: If True, only consider blocks with the explicit
                        "Calculated results:" prefix.
    Returns:
        (files_changed, blocks_filled_total)
    """
    load_dir = _existing_dir(load_dir)
    files_seen = 0
    files_changed = 0
    blocks_total = 0
    for root, dirs, files in os.walk(load_dir):
        for file in files:
            if not file.lower().endswith('.csv'):
                continue
            if file.startswith('._'):
                continue
            csv_path = os.path.join(root, file)
            files_seen += 1
            changed, blocks = fill_zprime_blocks_in_file(csv_path, require_prefix=require_prefix)
            if changed:
                files_changed += 1
                blocks_total += blocks
    if files_seen == 0:
        logger.warning(f"No CSV files found under '{load_dir}'. Is this the right folder?")
    elif files_changed == 0:
        logger.info(f"Found {files_seen} CSV file(s) in '{load_dir}': no empty Z' blocks, nothing changed")
    else:
        logger.info(f"Found {files_seen} CSV file(s) in '{load_dir}': filled {blocks_total} empty Z' "
                    f"block(s) in {files_changed} file(s)")
    return files_changed, blocks_total


import shutil

def extract_rpt_data(infile, samples_only=False, outdir="./"):
    """Extracts data from .rpt file and saves individual sample data to ./samples/ directory.
    
    Args:
        infile (str): Path to .rpt file.
        samples_only (bool): If True, only extract individual sample files.
        outdir (str): Output directory to save extracted data.

    Returns:
        pd.DataFrame: DataFrame containing extracted sample information."""
        
    file=infile
    print('\n\n',file)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # extract key line numbers using grep
    cmd=f"grep -n '^Sample\t[0-9]' '{file}' > {os.path.join(outdir, 'samples.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^Well\t' '{file}' > {os.path.join(outdir, 'wells.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^FileName\t' '{file}' > {os.path.join(outdir, 'fnames.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^BPM\t' '{file}' > {os.path.join(outdir, 'bpms.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^BPI\t' '{file}' > {os.path.join(outdir, 'bpis.txt')}"
    os.system(cmd)

    cmd=f"wc -l < '{file}' > {os.path.join(outdir, 'linenum.txt')}"
    os.system(cmd)

    cmd=f"grep -n '^;Mass' '{file}' > {os.path.join(outdir, 'massstart.txt')}"
    os.system(cmd)

    search='}'
    cmd=f"grep -n {search} '{file}' > {os.path.join(outdir, 'massend.txt')}"
    os.system(cmd)

    # read in key line numbers
    samples = pd.read_csv(os.path.join(outdir, 'samples.txt'), sep=':', header=None)
    wells = pd.read_csv(os.path.join(outdir, 'wells.txt'), sep=':', header=None)
    fnames = pd.read_csv(os.path.join(outdir, 'fnames.txt'), sep=':', header=None)
    bpms = pd.read_csv(os.path.join(outdir, 'bpms.txt'), sep=':', header=None)
    bpis = pd.read_csv(os.path.join(outdir, 'bpis.txt'), sep=':', header=None)

    samples.columns=['linenum','Sample']
    wells.columns=['linenum','Plate','Well']
    fnames.columns=['linenum','Fname']
    bpms.columns=['linenum','BPM']
    bpis.columns=['linenum','BPI']

    samples.Sample=samples.Sample.str.replace('Sample\t','').astype(int)
    wells=wells[wells.Well.notna()]
    wells.Plate=wells.Plate.str.replace('Well\t','').astype(int)

    fnames.Fname=fnames.Fname.str.replace('FileName\t','')
    bpms.BPM=bpms.BPM.str.replace('BPM\t','').astype(float)
    bpis.BPI=bpis.BPI.str.replace('BPI\t','').astype(float)

    wells['Sample']=np.nan
    fnames['Sample']=np.nan
    bpms['Sample']=np.nan
    bpis['Sample']=np.nan

    starts=pd.read_csv(os.path.join(outdir, 'massstart.txt'), sep=':', header=None)
    starts.columns=['linenum','headers']

    ends=pd.read_csv(os.path.join(outdir, 'massend.txt'), sep=':', header=None)
    ends.columns=['linenum','value']

    starts['Sample']=np.nan
    ends['Sample']=np.nan

    nl=pd.read_csv(os.path.join(outdir, 'linenum.txt'), header=None)
    linenum=nl.iloc[0][0]+1

    nsamples=len(samples)

    # assign samples to wells, fnames, bpms, bpis, starts, ends
    for i in range(nsamples):
        if i==nsamples-1:
            linerange=range(samples.linenum.loc[i], linenum)
        else:
            linerange=range(samples.linenum.loc[i], samples.linenum.loc[i+1]+1)
        sample=samples.Sample.loc[i]
        wells.loc[wells.linenum.isin(linerange), 'Sample']=sample
        fnames.loc[fnames.linenum.isin(linerange), 'Sample']=sample
        bpis.loc[bpis.linenum.isin(linerange), 'Sample']=sample
        bpms.loc[bpms.linenum.isin(linerange), 'Sample']=sample
        starts.loc[starts.linenum.isin(linerange), 'Sample']=sample
        ends.loc[ends.linenum.isin(linerange), 'Sample']=sample

    # get end rows for each sample, account for multiple end row matches
    endrows=[]
    for i in range(len(starts)):
        if i==len(starts)-1:
            linerange=range(starts.linenum.loc[i], linenum)
        else:
            linerange=range(starts.linenum.loc[i], starts.linenum.loc[i+1]+1)
        tmp=pd.DataFrame(ends[ends.linenum.isin(linerange)].iloc[0]).T
        endrows.append(tmp)
    endrows=pd.concat(endrows)

    # merge all together
    df=samples.merge(wells, how='left', on='Sample', suffixes=['_sample','_well'])
    df=df.merge(fnames, how='left', on='Sample')
    df=df.merge(bpis, how='left', on='Sample', suffixes=['_fname','_bpi'])
    df=df.merge(bpms, how='left', on='Sample')
    df=df.merge(starts, how='left', on='Sample', suffixes=['_bpm','_start'])
    df=df.merge(endrows, how='left', on='Sample')
    df=df[['Sample', 'Plate', 'Well', 'Fname', 'BPI', 'BPM', 'linenum_start', 'linenum']]
    df.columns=['Sample', 'Plate', 'Well', 'Fname', 'intensitymax', 'massatmax', 'linenum_start', 'linenum']
    df.linenum_start=df.linenum_start.replace(np.nan, -99).astype(int)

    # extract individual sample data files
    if not os.path.exists(os.path.join(outdir, 'samples')):
        os.makedirs(os.path.join(outdir, 'samples'))
    
    file=infile
    
    with open(file, "rb") as f:
        content=f.readlines()

    for i, row in df.iterrows():
        if row.linenum_start==-99:
            continue
        if samples_only:
            samplefile=os.path.join(outdir, 'samples', f'{row.Fname}.txt')
        else:
            samplefile=os.path.join(outdir, 'samples', f'Sample_{str(row.Sample).zfill(3)}_data.txt')
        with open(samplefile, "wb") as f:
            f.writelines(content[row.linenum_start-1:row.linenum-1])
    
    # cleanup intermediate files
    if samples_only:
        os.remove(os.path.join(outdir, 'samples.txt'))
        os.remove(os.path.join(outdir, 'wells.txt'))
        os.remove(os.path.join(outdir, 'fnames.txt'))
        os.remove(os.path.join(outdir, 'bpms.txt'))
        os.remove(os.path.join(outdir, 'bpis.txt'))
        os.remove(os.path.join(outdir, 'linenum.txt'))
        os.remove(os.path.join(outdir, 'massstart.txt'))
        os.remove(os.path.join(outdir, 'massend.txt'))
        df[['Sample', 'Plate', 'Well', 'Fname', 'intensitymax', 'massatmax',]].to_csv(os.path.join(outdir, 'samples', 'extracted_samples.csv'), index=False)

        samples_dir = os.path.join(outdir, 'samples')
        zip_base = os.path.join(outdir, 'samples')
        archive_path = shutil.make_archive(zip_base, 'zip', root_dir=samples_dir)
        print(f"\n\nExtracted and saved {len(df)} samples to the Samples folder;\nzipped folder is at {archive_path}.\n\n")
        return df[['Sample', 'Plate', 'Well', 'Fname', 'intensitymax', 'massatmax',]]
    return df



def trim_csv_columns_in_file(csv_path: str, ncols: int = 25) -> Tuple[bool, int]:
    """Cut every row of a CSV to its first ncols fields, in place.

    For EnVision exports where someone typed in columns past the plate
    (row label + 24 columns = 25). Blank lines and quoting are kept, and the
    file is only rewritten when a row was longer than ncols.

    Returns:
        (changed, rows_trimmed)
    """
    with open(csv_path, newline='') as f:
        text = f.read()
    newline = '\r\n' if '\r\n' in text else '\n'
    rows = list(csv.reader(text.splitlines()))
    trimmed = sum(len(r) > ncols for r in rows)
    if not trimmed:
        return (False, 0)
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f, lineterminator=newline)
        for r in rows:
            writer.writerow(r[:ncols])
    logger.info(f"Trimmed {trimmed} row(s) of {csv_path} to {ncols} columns")
    return (True, trimmed)


def _wells_per_read(csv_path, plate_type='384'):
    """Filled-well count of each read (measurement block) in one EnVision CSV."""
    import tempfile, shutil
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(csv_path, os.path.join(tmp, os.path.basename(csv_path)))
        try:
            parsed = parse_envision_files(plate_type=plate_type, load_dir=tmp)
        except Exception as e:  # unreadable file: say so rather than stop the tidy
            logger.warning(f"Could not parse {csv_path} to count wells per read: {e}")
            return []
    meas = [c for c in parsed.columns if c not in ('plate', 'row', 'col', 'well', 'Well')]
    return [int(parsed[c].notna().sum()) for c in meas]


def pad_single_background_row_in_file(csv_path: str) -> Tuple[bool, int]:
    """Where a block's Background information has one data row instead of
    two, add a blank line after it so row A of the results grid is on line 12
    of the block, as for two-row exports. Edits in place only when needed.

    Returns:
        (changed, blocks_padded)
    """
    with open(csv_path, newline='') as f:
        text = f.read()
    newline = '\r\n' if '\r\n' in text else '\n'
    rows = list(csv.reader(text.splitlines()))

    def blank(i):
        return i >= len(rows) or not any(str(v).strip() for v in rows[i])

    def first(i):
        return str(rows[i][0]).strip() if i < len(rows) and rows[i] else ''

    padded, out, i = 0, [], 0
    while i < len(rows):
        out.append(rows[i])
        # Background information, header, ONE row, blank, results title
        if (first(i).startswith('Background information') and not blank(i + 2) and blank(i + 3)
                and first(i + 4).startswith(('Calculated results', 'Results'))):
            out += rows[i + 1:i + 4]
            out.append([''] * max(1, len(rows[i + 3])))
            padded += 1
            i += 4
            continue
        i += 1
    if not padded:
        return (False, 0)
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f, lineterminator=newline).writerows(out)
    logger.info(f"Added a blank line after the single Background information row in {padded} block(s) of {csv_path}")
    return (True, padded)


def tidy_envision_files(paths, ncols: int = 25) -> pd.DataFrame:
    """Tidy the given EnVision CSVs in place and check what can't be tidied.

    Fixes: trims extra columns a user added after the plate (ncols = row
    label + 24), adds a blank line after a single Background information row
    (so row A is on line 12, as with two rows), and fills empty Z' blocks. Checks (reported, not fixed): the
    block layout (row A on line 12 counting each block's 'Plate information'
    line as 1, i.e. Plate and Background information first, one blank line
    between blocks), values pasted outside the grid, and reads covering
    different numbers of wells. See check_envision_csv.

    Meant for the raw files make_dr_platemaps_and_folders copied into the
    toload folders (plate_maps.raw_files(info)), so only the files that will
    be uploaded are touched and the originals on the share are left alone.
    Non-CSV paths (e.g. .rpt) are skipped.

    Returns:
        DataFrame with one row per CSV: rows trimmed and Z' blocks filled (0 =
        nothing to fix), wells per read, and problems (empty = compliant).
    """
    paths = [str(p) for p in ([paths] if isinstance(paths, (str, os.PathLike)) else paths)]
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(f"File(s) not found: {missing}")
    csvs = [p for p in paths if p.lower().endswith('.csv') and not os.path.basename(p).startswith('._')]
    skipped = len(paths) - len(csvs)
    rows = []
    for p in csvs:
        _, trimmed = trim_csv_columns_in_file(p, ncols)
        _, bg_padded = pad_single_background_row_in_file(p)
        _, zblocks = fill_zprime_blocks_in_file(p)
        counts = _wells_per_read(p)
        problems = check_envision_csv(p)
        rows.append({'file': p, 'rows trimmed': trimmed, 'background rows padded': bg_padded,
                     "Z' blocks filled": zblocks,
                     'wells per read': '/'.join(map(str, counts)), 'reads differ': len(set(counts)) > 1,
                     'problems': '; '.join(problems)})
    cols = ['file', 'rows trimmed', 'background rows padded', "Z' blocks filled", 'wells per read',
            'reads differ', 'problems']
    out = pd.DataFrame(rows, columns=cols)
    fixed = ['rows trimmed', 'background rows padded', "Z' blocks filled"]
    changed = int((out[fixed].sum(axis=1) > 0).sum()) if len(out) else 0
    for _, r in out[out['problems'] != ''].iterrows():
        logger.warning(f"{r['file']} doesn't look like an untouched EnVision export, and tidying can't fix it: "
                       f"{r['problems']}. Fix it by hand (leave missing wells blank, not 0).")
    note = f" (skipped {skipped} non-CSV file(s))" if skipped else ""
    if not csvs:
        logger.warning(f"No CSV files to tidy{note}.")
    elif changed:
        logger.info(f"Found {len(csvs)} CSV file(s){note}: tidied {changed}, {len(csvs) - changed} already tidy")
    else:
        logger.info(f"Found {len(csvs)} CSV file(s){note}: all already tidy, nothing changed")
    return out


_GRID_ROW = re.compile(r'^[A-Z]{1,2}$')


_PP_GRID_LINE = re.compile(r'^(\w[,\t])?[\d\.\-\t,]*$')


def _grid_row_labels(nrows):
    """Row letters for a plate grid: A..P for 384, A..Z then AA..AF for 1536."""
    labels = []
    for number in range(1, nrows + 1):
        label = ''
        while number:
            number, rem = divmod(number - 1, 26)
            label = chr(65 + rem) + label
        labels.append(label)
    return labels


def _block_layout_problems(lines, nrows=16):
    """Check each data block's layout, counting its 'Plate information' line as 1:

        1  Plate information          7-8  background rows
        2  header row                 9    blank
        3  plate row                  10   Calculated results... / Results for...
        4  blank                      11   column numbers
        5  Background information     12   row A
        6  header row                 ...  then a blank line after the last row

    i.e. Plate information and Background information are exported before
    the results, with one blank line between blocks. Returns a list of problems.
    """
    def blank(i):
        return i >= len(lines) or not any(str(v).strip() for v in lines[i])

    def first(i):
        return str(lines[i][0]).strip() if i < len(lines) and lines[i] else ''

    problems = []
    starts = [i for i in range(len(lines)) if first(i).startswith('Plate information')]
    stop = next((i for i in range(len(lines)) if first(i).startswith('Basic assay information')), len(lines))
    results = [i for i in range(stop) if first(i).startswith(('Calculated results', 'Results for'))]
    if results and (not starts or results[0] < starts[0]):
        problems.append(f"line {results[0] + 1}: a results block comes before any Plate information. Export with "
                        "Plate information and Background information placed before the results.")
    for b in starts:
        n = b + 1  # 1-based line number of this block's Plate information
        expect = [(1, lambda i: first(i).startswith('Plate'), "the Plate information header row"),
                  (3, blank, "a blank line"),
                  (4, lambda i: first(i).startswith('Background information'), "'Background information'"),
                  (8, blank, "a blank line"),
                  (9, lambda i: first(i).startswith(('Calculated results', 'Results')), "the results title"),
                  (11, lambda i: first(i) == 'A', "row A")]
        bad = [(off, what) for off, test, what in expect if not test(b + off)]
        if bad:
            off, what = bad[0]
            a_row = next((i for i in range(b + 1, min(stop, b + 40)) if first(i) == 'A'), None)
            where = f"row A is on line {a_row - b + 1} of the block" if a_row is not None else "no row A found"
            single_bg = (first(b + 4).startswith('Background information') and not blank(b + 6)
                         and blank(b + 7) and first(b + 8).startswith(('Calculated results', 'Results')))
            hint = (" The Background information has one row: the tidy cell adds the missing blank line."
                    if single_bg else " Plate information and Background information must come first, "
                    "each followed by one blank line.")
            problems.append(f"line {n} (block starting here): expected {what} on line {off + 1} of the block, "
                            f"so row A lands on line 12 ({where}).{hint}")
            continue
        # Every grid row keeps its letter even when its wells are blank, so a
        # blank row P still reads 'P'. Check the letters in order for each block
        # to tell a shifted block from a plate with an empty last row.
        wrong = next(((i, lab) for i, lab in enumerate(_grid_row_labels(nrows))
                      if first(b + 11 + i) != lab), None)
        if wrong is not None:
            i, lab = wrong
            got = first(b + 11 + i)
            problems.append(f"line {b + 12 + i}: expected row {lab} of the grid that starts on line {b + 12}, "
                            f"found {got!r}" + (" (blank line)." if blank(b + 11 + i) else "."))
            continue
        # Pipeline Pilot's EnVision Reader keeps only grid lines matching this
        # pattern and letters rows by their order, not by the label in the file.
        # A row with text (e.g. 1.2E+05, NaN, #DIV/0!) is dropped and every row
        # below it moves up one letter. Its 16-row limit means 1536 is skipped.
        if nrows <= 16:
            for i in range(nrows):
                cells = lines[b + 11 + i]
                if not any(str(v).strip() for v in cells[1:]):
                    continue  # an empty row (label only) is kept by PP as empty wells
                line = ','.join(str(v) for v in cells)
                if not _PP_GRID_LINE.match(line):
                    problems.append(f"line {b + 12 + i}: row {first(b + 11 + i)} has a value that is not a plain "
                                    "number (for example text, NaN or 1.2E+05). Pipeline Pilot drops this row and "
                                    "shifts the rows below it up one letter.")
                    break
        last = b + 11 + nrows - 1
        if not blank(last + 1):
            problems.append(f"line {last + 2}: expected a blank line after the results grid that starts on line "
                            f"{b + 12}.")
    return problems


def check_envision_csv(csv_path: str, plate_type: str = None, max_problems: int = 10) -> list:
    """Check that an EnVision CSV still looks like an instrument export.

    Catches the hand edits that make Pipeline Pilot's EnVision reader fail:
    values pasted outside the plate grid or over the Plate information /
    Background information lines, extra rows or columns in a plate grid, and
    reads (measurement blocks) that cover different numbers of wells, which is
    what makes Merge Plates fail with "Well property size mismatch". Only the
    part before "Basic assay information" is checked; the protocol details
    after it have no fixed layout. The plate format (96/384/1536) is read from
    the first plate grid's column numbers unless plate_type is given.

    Returns:
        list[str]: problems, each naming the line. Empty = looks compliant.
    """
    with open(csv_path, newline='') as f:
        lines = list(csv.reader(f))
    problems = []
    if plate_type is None:
        widths = [sum(1 for v in r[1:] if str(v).strip()) for r in lines
                  if r and not r[0].strip() and len(r) > 1 and str(r[1]).strip().lstrip('0') == '1']
        by_cols = {v[2]: k for k, v in plate_dict.items()}
        plate_type = next((by_cols[w] for w in widths if w in by_cols), '384')
    ncols, nrows = plate_dict[str(plate_type)][2], plate_dict[str(plate_type)][1]

    def filled(row, start=0):
        return [i for i, v in enumerate(row) if i >= start and str(v).strip()]

    section, width, grid_rows = None, 0, 0
    for n, row in enumerate(lines, start=1):
        first = row[0].strip() if row else ''
        if first.startswith('Basic assay information'):
            break
        cells = filled(row)
        if not cells:                                   # blank line ends a section
            section = None
            continue
        if first.startswith(('Plate information', 'Background information')):
            section, width = 'info_title', 0
            extra = filled(row, 1)
        elif first.startswith(('Calculated results', 'Results for')) or (first and section is None
                                                                           and not _GRID_ROW.match(first)):
            section = 'title'
            extra = filled(row, 1)
        elif section == 'info_title':                   # header row: its named columns are the limit
            width = next((i for i, v in enumerate(row) if not str(v).strip()), len(row))
            section = 'info'
            extra = filled(row, width)
        elif section == 'info':
            extra = filled(row, width)
        elif first == '' and [str(row[i]).strip().lstrip('0') for i in cells] == [str(k) for k in range(1, len(cells) + 1)]:
            section, grid_rows = 'grid', 0              # column-number header of a plate grid
            extra = filled(row, ncols + 1)
        elif section == 'grid' and _GRID_ROW.match(first):
            grid_rows += 1
            extra = filled(row, ncols + 1)
            if grid_rows > nrows:
                problems.append(f"line {n}: plate grid has more than {nrows} rows")
        else:
            extra = cells                               # values where no data belongs
        if extra:
            cols = ', '.join(str(i + 1) for i in extra[:6]) + (' ...' if len(extra) > 6 else '')
            problems.append(f"line {n} ({first[:30] or 'blank first cell'}): unexpected values in column(s) {cols}")
        if len(problems) >= max_problems:
            problems.append('... more problems not listed')
            return problems

    problems += _block_layout_problems(lines, nrows)
    counts = _wells_per_read(csv_path, plate_type)
    if len(set(counts)) > 1:
        problems.append(f"reads have data in different numbers of wells ({'/'.join(map(str, counts))}); "
                        "Pipeline Pilot's Merge Plates fails on this ('Well property size mismatch')")
    return problems
