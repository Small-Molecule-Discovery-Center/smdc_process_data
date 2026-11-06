# plate_maps.py reference

Tooling for building compound plate maps for the SMDC high-throughput screening
workflow, and for staging screen results so they can be loaded into HiTS.

The job this solves: a screen produces two kinds of information that have to be
reconciled before anything can be uploaded.

1. **What the scientist intended**, recorded by hand in an Excel workbook. One
   sheet per destination plate, drawn as grids that look like the physical plate.
2. **What the inventory database knows**: which compound lot physically sits in
   which well of which source plate.

This module parses (1), looks up (2), joins them, and writes out the per-plate
manifests plus the folder tree the HiTS loader expects.

---

## Where this lives now

This reference came over from the retired data-proc repo (Nov 2026). The code is
`database_interactions/plate_maps.py` in smdc_process_data. The notebook that drives it
is `database_interactions/prepare_hits_upload.ipynb`.

## Contents

- [Setup](#setup)
- [Quick start](#quick-start)
- [The input workbook](#the-input-workbook)
- [The `Source_plate` dispatch language](#the-source_plate-dispatch-language)
- [How lots get resolved](#how-lots-get-resolved)
- [Empty wells and zero concentrations](#empty-wells-and-zero-concentrations)
- [What gets written](#what-gets-written)
- [Function reference](#function-reference)
- [Tests](#tests)
- [Troubleshooting](#troubleshooting)
- [Known bugs and rough edges](#known-bugs-and-rough-edges)

---

## Setup

The utilities are notebook-driven and can also be installed as modules from a
local checkout. Database access expects an open SSH tunnel to the local MySQL
port.

### Python dependencies

From a checkout of this repo: `python -m pip install -e '.[upload]'`.

### Three things you must have before anything works

| Requirement | Detail |
|---|---|
| **Database credentials** | `mysql_user`, `mysql80_pwd` and `mysql80_port` in the environment (see `mysql_cnx.create_mysql_cnx`). Reads go to `QB3_HiTS`, or to `hits_db` if that is set. |
| **SSH tunnel to HiTS** | Open before you start. If HiTS can't be reached, the first lookup stops with `ConnectionError: Cannot reach HiTS ... Is the SSH tunnel up`. |
| **The Montreal share** | Mounted. The code finds it in this order: `smdc_share_root`, `/Volumes/Shared/` (Mac), `/mnt/z/` (WSL), `Z:/` (Windows). You can also pass `share_root=`. |

### Import

```python
import smdc_paths   # only needed outside database_interactions/
from plate_maps import (
    make_dr_platemaps_and_folders,
    make_platemaps_from_barcode,
    draw_plate_heatmap,
    clean_up_folders,
)
```

Importing configures root logging at INFO and sets a seaborn theme. That is
deliberate for notebook use, and it means importing this into anything larger
will change that program's logging level.

The pipeline itself is **non-interactive**: nothing prompts on stdin, and it runs
from a script as well as from a notebook.

Plates that share a raw data file, assay and screen date share one `toload_<i>`
folder and one `grid_num.txt`. Each plate still uses its own `Source_plate` and
`Lot_info` from the `Info` sheet.

---

## Quick start

The common case, a dose-response screen:

```python
info = make_dr_platemaps_and_folders(
    infile='/path/to/my_screen.xlsx',      # the workbook, see below
    users=[("Firstname", "Lastname")],     # whoever registered the compound lots
    basedir='/path/to/hits_data_uploads/', # where the folder tree gets built
    dose_response=True,
    replace_dict={},                       # extra control labels, see below
)
```

This will:

1. Read the `Info` sheet to find out which plates were run, when, and where the
   raw data lives.
2. For each plate, resolve every well to a compound, lot, source well and
   concentration.
3. Write `<grid_num>_DR_plate.xlsx` manifests.
4. Build `<basedir>/<assay>/toload_<i>/<date>/` and copy the raw data in.

Then, once HiTS ingest is confirmed:

```python
clean_up_folders(
    infile='/path/to/my_screen.xlsx',  # MUST be the same file as above
    basedir='/path/to/hits_data_uploads/',
    delete_orig=False,                 # leave False until you trust the copy
)
```

Before or after, to eyeball a plate:

```python
draw_plate_heatmap(barcode='OTPARP1234-Eco1', values='SMDC_ID', labels='SMDC_ID')
```

---

## The input workbook

> `TODO(template)`: link and describe the example workbooks in the
> `hits_data_uploads` Drive folder here, and say which one to copy for which kind
> of screen. The format below is what the code enforces, derived from the code
> rather than from the templates.

One `.xlsx` file with:

- **One `Info` sheet** (or `Sheet1`), the driver table. One row per destination plate.
- **One sheet per destination plate**, named after that plate's `Grid_num`.

### The `Info` sheet

| Column | Meaning |
|---|---|
| `Grid_num` | Destination plate barcode. **Also the name of that plate's sheet.** |
| `Source_plate` | Where the compounds came from. This is a small dispatch language, not just a barcode. See below. |
| `Raw_assay_ID` | The HiTS assay ID. Either a single number for the whole plate, or the literal word `map` if it varies per well. |
| `Lot_info` | Either blank, or `use map` to read lots from the grid. |
| `Date_of_screen` | Date the plate was run. Becomes a folder name. |
| `Raw_data_file` | Path to the instrument output. Typed in by hand; see `_normalize_share_path`. |
| `Assay_name` | Human-readable name. Becomes part of a folder name. |

Column names are matched **case- and punctuation-insensitively**, so `Grid num`,
`grid_num` and `GridNum` all work. See `_standardize_info_columns` for the exact
set of accepted spellings.

Rows with no `Date_of_screen` are dropped. A workbook with a title row above the
real header is handled.

### The plate sheets: the block convention

**This is the thing to understand about the format.** A single sheet stacks
several plate-shaped grids vertically. Each block has a header row of column
numbers, then one row per plate row with the row label in the first column.
Completely blank rows are dropped before reading, so blank lines between blocks
don't matter:

```
      1    2    3   ...  24     <- header row
 A  [                      ]
 B  [                      ]   <- BLOCK 0: SMDC IDs
...  [                      ]
 P  [                      ]
      1    2    3   ...  24     <- header row
 A  [                      ]
 B  [                      ]   <- BLOCK 1: concentrations in uM
...  [                      ]
 P  [                      ]
      1    2    3   ...  24     <- header row
 A  [                      ]
 B  [                      ]   <- BLOCK 2: context-dependent, see below
...  [                      ]
 P  [                      ]
```

Blocks are found by their row labels (`A` to `P`, or `A` to `AF` for 1536) in the
first column, so a sheet that leaves out the per-block header rows also reads
correctly. If a populated row is added inside a block, that block's labels don't
run A to P, it isn't found, and the run stops with an error that says so. (Before
Nov 2026 the code used a fixed stride of rows + 1. That was right for the real
maps, which have header rows, but would have shifted a sheet without them.)

**Block 1 is optional.** A cherry-pick map may leave it out, or use it for a
plate suffix instead of concentrations. In both cases `[compound] uM` is left
blank, and the suffix and a whole-plate concentration are passed in when the
plate is registered. A block 1 that mixes numbers and text stops with an error.

**Block 2 is overloaded.** Its meaning depends on `Source_plate` and `Lot_info`:

| Condition | Block 2 holds |
|---|---|
| `Lot_info` is `use map` | Lot numbers |
| `Source_plate` contains `map` | Source plate barcodes |
| `Raw_assay_ID` is `map` | Assay IDs |

A consequence worth knowing: **a plate can only have one per-well third field.**
You cannot have both per-well lots and per-well assay IDs on the same plate.

### Shorthand in block 0

- An SMDC ID cell may be written `285699-2` to pin an exact lot (spaces around
  the dash are fine). Only cells that are exactly `<digits>-<digits>` are split.
  Pinned lots always win; the other wells get their lot the usual way (block 2,
  source plate lookup, then newest lot by user). If every compound well is
  pinned, block 2 isn't needed. A cell like `12346-1-2` stops with an error.
- Non-compound wells: `DMSO`, `Water`, `Background` (any case) are recognised
  automatically and nulled out.
- Any other text label will raise a `ValueError` naming the offending string.
  That is the intended failure mode. Add it to `replace_dict`:

  ```python
  replace_dict={'MyControl': np.nan, 'PositiveCtrl': 285699}
  ```

  Values may be `np.nan` (not a compound) or a real SMDC ID.

### Sheet name truncation

Excel caps sheet names at 31 characters. The code slices `Grid_num[0:31]` when
looking up the sheet. If two of your barcodes agree in their first 31
characters, they will collide.

---

## The `Source_plate` dispatch language

The `Source_plate` cell tells the pipeline how to find where compounds came from.
Matched as a **case-insensitive substring**:

| Value | Behaviour |
|---|---|
| `n.a.` | No source plate exists (compounds added by hand). Source plate and well become the literal string `NA`. |
| `lot` (e.g. `use lots`) | Lots are drawn in block 2. Source plates are then looked up from (SMDC_ID, Lot). |
| blank | **An error** for dose-response uploads, raised before any folder is made. It used to mean `use lots`, which hid a forgotten Source_plate until the spot-check. |
| `map` | Source plate barcodes are drawn in block 2, i.e. wells came from different plates. |
| `biomek` | The sheet is already a long-format Biomek worklist, not a block grid. Source wells are read directly from it; block parsing is skipped entirely. This branch also **overwrites** the worklist's SMDC ID column with the database's answer. |
| anything else | Treated as one or more literal barcodes, separated by `;` and/or newlines. |

**Watch out with literal barcodes**: every listed plate is pulled and joined on
`SMDC_ID`. If a compound appears on two of the listed plates, the join fans out
into duplicate rows. This is exactly what the "Was a compound pulled from two
plates?" warning catches.

---

## How lots get resolved

The lot is what makes a result reproducible: same compound, different synthesis
batch. It is also the field most often missing from a hand-drawn map, hence a
fallback ladder. Each stage only fills wells the previous ones left empty, in
decreasing order of confidence:

1. **Already present.** Drawn in the map, or split from a `285699-2` cell.
2. **Read from the map or the obvious lookup**, depending on the dispatch value:
   `use map` reads block 2; `map` looks up by (SMDC_ID, source plate); `biomek`
   looks up by (source plate, source well).
3. **By (SMDC_ID, any known source plate).**
4. **By user**: the newest lot registered by whoever you passed in `users`.
5. **Give up**: leave empty, log a warning, write `NA`.

**Stage 4 is an assumption, not a fact.** It attaches the newest lot that person
registered, with no constraint that the lot was ever on the plate. It logs at
INFO and is not distinguishable in the output. If lot provenance matters for what
you are doing, check the log for `Getting lot info from user` and confirm those
wells by hand.

---

## Empty wells and zero concentrations

A zero in the concentration column is not loadable: dose-response curves are fit
on a log scale, so zero has no place on the axis. In practice a zero means the
well held no compound.

**Zero-concentration wells are dropped from the manifest**, which is how an empty
well is represented. This happens automatically, with a log line:

```
INFO: Dropping 16 well(s) with zero concentration; treating them as empty wells.
```

Nothing is asked of you. An earlier version prompted on stdin for a replacement
value; pressing Enter or typing `NA` produced exactly this behaviour and was the
recommended answer, so the prompt is gone and the recommended answer is now the
only one.

The drop happens **after** the plate row-count check, so the count still sees the
full plate. If you see the "Was a compound pulled from two plates?" warning, it is
about a genuine fan-out, not about your empty wells.

If you ever have zeros that genuinely mean "below the lowest dilution" rather than
"empty", fix them in the workbook. Do not reintroduce the prompt.

---

## What gets written

```
<basedir>/
  <raw_assay_id>_<assay_name>/
    toload_<i>/
      <YYYYMMDD>/
        <grid_num>_DR_plate.xlsx   <- the manifest, one row per well
        grid_num.txt               <- which plates this folder's data covers
        <raw_data_file>            <- copied in from /Volumes/Shared/
    loaded/                        <- clean_up_folders moves things here after ingest
```

### The manifest

One row per destination well, columns:

```
SMDC_ID | Lot | Source plate | Source well | Destination plate | Destination well | [compound] uM | Assay ID
```

Empty fields are the **literal string `NA`**, which is what the HiTS loader
expects. Not blank, not NaN.

### The `toload_<i>` numbering

Worth knowing because it is not obvious:

- If `Raw_data_file` ends in **`.rpt`**, `i` is pinned to `0`. One report covers
  the whole run, so all its plates share `toload_0`.
- Otherwise `i` is the row index within that file's group, giving `toload_0`,
  `toload_1`, and so on.

Both the pipeline and `clean_up_folders` get this from `_iter_load_groups`, and
both get the folder name from `_assay_dir_name`, so they cannot disagree about
where a folder lives. **If you change the layout, change those two helpers, not
the call sites.**

You must still pass `clean_up_folders` **the same workbook**, since it recomputes
the paths from it.

---

## Function reference

### Entry points

| Function | Purpose |
|---|---|
| `make_dr_platemaps_and_folders(infile, users, basedir, dose_response, replace_dict)` | The main pipeline. Everything above. |
| `make_platemaps_from_barcode(barcode, outdir, cnx)` | Given a physical plate barcode, pull its contents from HiTS and write an Excel workbook: one sheet per attribute in plate geometry, plus a long-form sheet. Pass `outdir=None` to get the DataFrame back without writing. |
| `draw_plate_heatmap(barcode, df, wells, values, labels, ...)` | Visual QC of one plate. Reads from the database if you give it a barcode, or plots a DataFrame you already have. |
| `clean_up_folders(infile, basedir, delete_orig)` | Post-ingest housekeeping. Destructive, no dry-run. |

### Useful pieces

| Function | Purpose |
|---|---|
| `plate_well_list(size)` | `['A01', 'A02', ... 'P24']` for a given plate format. |
| `assign_wells_by_idx(df, ...)` | Map a within-plate index column onto well names, optionally reserving control columns. |
| `read_platemap_block(plate_map, block, value_name, size)` | Read one grid block out of a sheet into long form. This is the parser for the block convention. |
| `lookup_id_and_lot_by_plate_and_well(plates, wells)` | Most trustworthy lookup: plate + well identifies exactly one physical well. |
| `lookup_lot_by_plate_and_ic_id(ids, plates, ...)` | Lot and source well for compounds known to sit on given plates. |
| `lookup_plate_by_lot_num(ids, lot_nums, ...)` | Which plates and wells hold a given (compound, lot) pair. |
| `lookup_lot_by_user(ids, users)` | Newest lot a person registered. Last resort only. |

### Shared plumbing

These exist so the pipeline and clean-up cannot drift apart. They are the single
source of truth for their piece; change them here, not at the call sites.

| Function | Purpose |
|---|---|
| `_read_info_sheet(infile)` | Read and normalise the Info sheet. Returns `(dp, info)`. |
| `_iter_load_groups(info)` | Yield `(grouped_df, rpt)` per raw data file. Owns the `toload_<i>` numbering. |
| `_assay_dir_name(raw_assay_id, assay_name)` | Build the `<assay_id>_<assay_name>` folder name. |
| `_normalize_share_path(raw)` | Normalise a hand-typed `Raw_data_file` to a path relative to the share root. |

The other `_`-prefixed functions handle column-name normalisation (`_norm_key`,
`_canonicalize_columns`, `_standardize_info_columns`, `_standardize_dids_columns`)
and control-label setup (`_build_replace_metadata`). The `get_*_info` functions
each own one column of the manifest and are called in dependency order by
`make_dr_platemaps`; that order matters, since each reads columns the previous
one added.

---

## Tests

`tests/test_plate_maps*.py` cover the parts that run without a database:

```bash
python -m pytest tests/
```

It stubs `mysql_cnx` (and, if absent, `humanize` / `seaborn`), then checks:

- `_normalize_share_path` reproduces the original inline logic on 13 real-world
  path shapes, so the extraction was faithful.
- The one deliberate divergence, `Montreal/Shared/`, now yields the corrected
  result. The old buggy output is asserted too, as documentation of what changed.
- **Every prefix branch strips exactly its own length.** This is the check that
  catches the whole class of off-by-one the `Montreal/Shared/` bug belonged to;
  add a case here whenever you add a branch.
- No result starts with `/`, which would make `os.path.join` discard the share
  root.
- `_assay_dir_name` gives the same answer for `123`, `123.0`, `'123'` and
  `'123.0'`, the spellings that used to break clean-up.
- `handle_zero_concentrations` drops zeros without prompting.
- `read_platemap_block` rejects an unknown plate size with a clear error.

Everything else needs a live HiTS connection and is untested. The database-facing
functions are not written to be testable; see the `TODO: decompose further` on
`make_platemaps_from_barcode`.

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'plate_maps'` (or `mysql_cnx`)**
Run `pip install -e '.[upload]'` from this repo, and add `import smdc_paths`
before the import when the notebook is outside `database_interactions/`.

**`ConnectionError: Cannot reach HiTS ...`**
The SSH tunnel is down, or the `mysql_*` variables are wrong. Reopen the tunnel
and rerun.

**`ValueError: could not convert string to float: 'SomeLabel'`**
An unrecognised text label in block 0. Add it to `replace_dict`.

**`KeyError` on a sheet name**
The `Grid_num` in your `Info` sheet does not match a sheet name, or the barcode
exceeds 31 characters and Excel truncated the sheet name differently than you
expect.

**`KeyError: make_platemaps_from_barcode returned unexpected columns...`**
The `SELECT` list in `make_platemaps_from_barcode` changed and no longer matches
what `get_source_plate_info` renames. The message names the missing column. This
used to be a silent mislabelling.

**`FileNotFoundError` on the raw data**
The `Raw_data_file` path is typed by hand from whichever machine acquired the
data, so it may arrive as a Windows drive letter (`Z:\...`), a UNC path, an
`smb://` URL, a macOS `/Volumes` mount, or a raw IP. `_normalize_share_path`
strips each known prefix. **The first time someone uses a new machine, you will
need to add a branch there**, and a matching case in `test_plate_maps.py`. The
normalised path is logged; compare it against what actually exists under
`/Volumes/Shared/`.

The copy now uses `shutil.copy2` and **raises** on failure rather than scrolling
past, so a missing file stops the run instead of producing a silently empty
folder.

**"Your plate has N wells. Was a compound pulled from two plates?"**
A merge fanned out, usually because a compound was found on more than one listed
source plate. This is no longer triggered by empty wells.

**`clean_up_folders` says "No such folder, skipping"**
It looked for a `toload_<i>` that does not exist. Confirm you passed the same
workbook you gave the pipeline. Both compute the path the same way now, so a
genuine mismatch means the workbook changed between the two runs.


---

## Known bugs and rough edges

### Fixed in this version

Recorded so you know what changed if you compare against an older copy.

| | What was wrong |
|---|---|
| **Source plates never uppercased** | In `lookup_lot_by_plate_and_ic_id`, `.upper()` ran *before* `pd.concat`, so it applied to the last loop iteration's frame and was then discarded. Barcodes were never normalised and case-sensitive merges could silently drop rows. Now concat first, then uppercase. An empty result returns an empty frame with a warning rather than raising `NameError`. |
| **1536-well plates unreadable** | `read_platemap_block` tested `size==1596`, a typo, so `size=1536` fell through every branch and raised `UnboundLocalError`. Fixed, and an unknown size now raises a `ValueError` that says what it wanted. |
| **1536 row slice** | `draw_plate_heatmap` used `plate_rows[0:48]` where the list has 32 entries. Now `[0:32]`, canonical for 1536. |
| **Pipeline and clean-up disagreed on folders** | They computed the assay folder name and the `toload_<i>` index by different rules, so clean-up could fail to find folders that existed. Both now go through `_read_info_sheet`, `_iter_load_groups` and `_assay_dir_name`. |
| **Shell-outs** | `os.system("rm {path}")` with an *unquoted* path, and `os.system("cp '{path}'")` whose return code was ignored. Now `os.remove` and `shutil.copy2`: no shell, and failures raise. |
| **Positional column rename** | `df.columns = [...]` in `get_source_plate_info` silently mislabelled everything if the `SELECT` order changed. Now a dict rename with an explicit check that raises and names the missing column. |
| **Interactive prompt** | `handle_zero_concentrations` called `input()`. Now drops zeros as empty wells, and runs after the row-count check so that check stays meaningful. |
| **Copy duration** | Logged `start - end`, a negative number. Now `end - start`. |
| **`basicConfig(force=True)`** | Tore out whatever logging the host had configured. Now `basicConfig(level=INFO)` plus an explicit `setLevel`, which makes this module's output visible without stomping anyone else's handlers. |
| **`Montreal/Shared/` path branch** | Sliced `[17:]` against a 16-character prefix, eating the next character (`SMDC/run.txt` → `MDC/run.txt`), so any path arriving in that form could never resolve. Now `[16:]`. `test_plate_maps.py` asserts every branch strips exactly its own prefix length. |
| **Dead twin** | `_make_platemaps_from_barcode` and its four helpers (`fetch_plate_data`, `filter_common_name`, `format_plate_data`, `save_plate_layouts`) duplicated the public function and were never called. Removed. `make_platemaps_from_barcode` carries a `TODO: decompose further`. |

### Still outstanding

**1. Lot stage 4 attaches an unverified lot.** Covered under
[How lots get resolved](#how-lots-get-resolved). The "newest lot this user
registered" fallback is a guess recorded as fact, distinguishable only in the log.
If provenance matters, make it visible in the manifest.

Items 2 (swallowed database failure) and 3 (f-string SQL) were fixed in Nov 2026:
connection failures now raise, and all queries use bound parameters. The Copilot
review fixes from the same move are listed in smdc_process_data PRs #3 and #5.

**4. Smaller things.**

- `read_platemap_block` rebinds `plate_rows` as a local int, shadowing the
  module-level list of row labels. Confusing, not a live bug.
- `assign_wells_by_idx` does chained assignment on a slice, which will break under
  copy-on-write.
- `draw_plate_heatmap` infers plate size from observed rows/columns, so a sparse
  384w plate using only A-H/1-12 reads as 96w. Pass `plate_size` explicitly for
  partial plates.
- `draw_plate_heatmap` uses `str.strip(col)` to remove the column suffix, which
  strips *characters*, not a suffix. Works only because row letters and digits do
  not overlap. `removesuffix` would say what is meant.
- Importing the module still reconfigures root logging and sets a seaborn theme.
  Fine for notebooks, rude anywhere else.
