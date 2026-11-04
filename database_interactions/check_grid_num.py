#!/usr/bin/env python3
"""
Audit grid_num.txt files against EnVision plate-reader CSV exports.

The parser supports both common EnVision CSV styles:

    old:  plate/background information blocks, but no data/result grid
    new:  one or more Results for / Calculated results data blocks

Style is controlled by the presence of a result grid (`has_grid`), rather than
by the position of the metadata. In new exports, Plate information and
Background information may occur before *or after* the associated data block;
both positions are parsed.

Usage:
    python check_grid_num.py ROOT [-o report.csv] [-d detail.csv]
                                     [--min-score 0.55] [--strict-dates] [-v]
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import re
import sys
from dataclasses import dataclass, field
from datetime import datetime
from difflib import SequenceMatcher
from typing import Iterable

try:
    from scipy.optimize import linear_sum_assignment
    import numpy as _np
    _HAVE_SCIPY = True
except ImportError:  # pragma: no cover
    _HAVE_SCIPY = False

try:
    from rapidfuzz.distance import JaroWinkler
    _HAVE_RAPIDFUZZ = True
except ImportError:
    _HAVE_RAPIDFUZZ = False

GRID_FILENAMES = {"grid_num.txt"}
DATA_EXTS = {".csv", ".txt", ".xlsx", ".xls", ".json", ".xml", ".tsv", ".rpt"}
SNIFF_BYTES = 64_000
BLANK_BARCODE = {"", "N/A", "NA", "NONE", "NULL", "-", "0", "-0-",}

PLATE_INFORMATION_RE = re.compile(r"^\s*Plate information\s*$", re.I | re.M)
BACKGROUND_INFORMATION_RE = re.compile(r"^\s*Background information\s*$", re.I | re.M)
GRID_RE = re.compile(r"^\s*(Results for|Calculated results)\b", re.I | re.M)

FORMATS = (
    "%m/%d/%y %H:%M",
    "%m/%d/%Y %I:%M:%S %p",
)


def read_text(path: str, limit: int | None = None) -> str:
    """Read text defensively; EnVision exports commonly contain a BOM."""
    with open(path, "rb") as fh:
        raw = fh.read(limit) if limit else fh.read()
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            pass
    return raw.decode("latin-1", errors="replace")


def norm(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", (s or "").upper())


def digit_runs(s: str) -> list[int]:
    runs = re.findall(r"\d+", s or "")
    runs.sort(key=len, reverse=True)
    return [int(r) for r in runs]


def strip_trailing_alpha(s: str) -> str:
    return re.sub(r"[A-Z]+$", "", norm(s))


def similarity(a: str, b: str) -> float:
    if _HAVE_RAPIDFUZZ:
        return JaroWinkler.similarity(a, b)
    return SequenceMatcher(None, a, b).ratio()


def is_datetime(value: str) -> bool:
    if not isinstance(value, str):
        return False
    for fmt in FORMATS:
        try:
            datetime.strptime(value.strip(), fmt)
            return True
        except ValueError:
            pass
    return False


def parse_grid_num(path: str) -> list[str]:
    lines = [ln.strip() for ln in read_text(path).splitlines()]
    entries = [ln for ln in lines if ln]
    if entries and norm(entries[0]) in {"GRIDNUM", "GRID", "BARCODE", "PLATE", "PLATEID"}:
        entries = entries[1:]
    return entries


@dataclass
class PlateRecord:
    index: int | None
    barcode: str
    measurement_date: str
    blank_barcode: bool


@dataclass(frozen=True)
class EnvisionMetadataBlock:
    kind: str
    row: int
    index: int | None
    barcode: str
    measurement_date: str


def _clean_row(row: list[str]) -> list[str]:
    return [c.strip() for c in row]


def _is_blank_row(row: list[str]) -> bool:
    return not any(c.strip() for c in row)


def _parse_plate_block(rows: list[list[str]], start: int) -> EnvisionMetadataBlock | None:
    """Parse a Plate information block wherever it appears in the file."""
    if start + 2 >= len(rows):
        return None
    header = _clean_row(rows[start + 1])
    values = _clean_row(rows[start + 2])
    if not header or not values:
        return None
    rec = dict(zip(header, values))
    if "Plate" not in rec or "Barcode" not in rec:
        return None

    raw_index = (rec.get("Plate") or "").strip()
    barcode = (rec.get("Barcode") or "").strip('"= ')
    try:
        idx = int(raw_index)
    except (TypeError, ValueError):
        idx = None

    return EnvisionMetadataBlock(
        kind="plate",
        row=start,
        index=idx,
        barcode=barcode,
        measurement_date=(rec.get("Measurement date") or "").strip(),
    )


def _find_markers(rows: list[list[str]]) -> tuple[list[int], list[int], list[int]]:
    plate_rows: list[int] = []
    background_rows: list[int] = []
    grid_rows: list[int] = []
    for i, row in enumerate(rows):
        cell = row[0].strip() if row else ""
        if PLATE_INFORMATION_RE.match(cell):
            plate_rows.append(i)
        elif BACKGROUND_INFORMATION_RE.match(cell):
            background_rows.append(i)
        elif GRID_RE.match(cell):
            grid_rows.append(i)
    return plate_rows, background_rows, grid_rows


def classify_csv(path: str) -> tuple[str, str, bool]:
    """
    Return (kind, reason, has_grid).

    kind is one of envision_old, envision_new, tabular_csv, not_csv_like.
    EnVision style is selected by has_grid, never by metadata placement.
    """
    try:
        text = read_text(path)
    except OSError as exc:
        return "not_csv_like", f"unreadable: {exc}", False
    if "\x00" in text:
        return "not_csv_like", "binary content", False

    rows = list(csv.reader(io.StringIO(text, newline="")))
    _, _, grid_rows = _find_markers(rows)
    plate_marker = bool(_find_markers(rows)[0])
    has_grid = bool(grid_rows)

    if plate_marker:
        if has_grid:
            return "envision_new", "Plate information + result grids", True
        return "envision_old", "Plate information blocks without result grids", False
    if has_grid:
        return "envision_new", "result grids without visible Plate information", True

    first = next((ln for ln in text.splitlines() if ln.strip()), "")
    if "," in first or "\t" in first:
        return "tabular_csv", "flat delimited file", False
    return "not_csv_like", "no EnVision markers or delimiter", False


def parse_envision_plates(path: str) -> list[PlateRecord]:
    """
    Extract one record per physical plate from Plate information blocks.

    Metadata is located independently of result grids, so a new-style export
    works whether Plate/Background information occurs before or after the grid.
    Duplicate blocks for the same plate are merged, retaining the first useful
    barcode and date encountered in file order.
    """
    text = read_text(path)
    rows = list(csv.reader(io.StringIO(text, newline="")))
    plate_rows, _, _ = _find_markers(rows)

    by_index: dict[int | None, PlateRecord] = {}
    order: list[int | None] = []

    for start in plate_rows:
        block = _parse_plate_block(rows, start)
        if block is None:
            continue
        idx = block.index
        barcode = block.barcode
        date = block.measurement_date
        blank = barcode.upper() in BLANK_BARCODE or is_datetime(barcode)

        if idx not in by_index:
            by_index[idx] = PlateRecord(idx, barcode, date, blank)
            order.append(idx)
            continue

        prev = by_index[idx]
        if prev.blank_barcode and not blank:
            prev.barcode = barcode
            prev.blank_barcode = False
        if not prev.measurement_date and date:
            prev.measurement_date = date

    if all(k is not None for k in order):
        order = sorted(order)
    return [by_index[k] for k in order]


def dates_monotonic(plates: list[PlateRecord]) -> bool | None:
    parsed: list[datetime] = []
    formats = (
        "%m/%d/%Y %I:%M:%S %p",
        "%d/%m/%Y %H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
        "%m/%d/%Y %H:%M:%S",
    )
    for plate in plates:
        if not plate.measurement_date:
            return None
        for fmt in formats:
            try:
                parsed.append(datetime.strptime(plate.measurement_date, fmt))
                break
            except ValueError:
                pass
        else:
            return None
    return all(a <= b for a, b in zip(parsed, parsed[1:]))


def pair_score(grid_entry: str, barcode: str) -> tuple[float, str]:
    if not barcode:
        return 0.0, "barcode_blank"
    g_norm, b_norm = norm(grid_entry), norm(barcode)
    if not g_norm or not b_norm:
        return 0.0, "empty_after_normalisation"
    if g_norm == b_norm:
        return 1.00, "exact"
    g_core, b_core = strip_trailing_alpha(grid_entry), strip_trailing_alpha(barcode)
    if g_core and g_core == b_core:
        return 0.98, "exact_ignoring_suffix"
    g_runs, b_runs = digit_runs(grid_entry), digit_runs(barcode)
    if g_runs and b_runs:
        if g_runs[0] == b_runs[0]:
            return 0.95, "longest_number_equal"
        if set(g_runs) & set(b_runs):
            return 0.85, "shared_number"
    if b_norm in g_norm:
        return 0.80, "barcode_is_substring"
    return min(similarity(g_norm, b_norm) * 0.75, 0.74), "fuzzy_string"


def align(grid_entries: list[str], plates: list[PlateRecord], min_score: float) -> list[dict]:
    ng, nb = len(grid_entries), len(plates)
    scores = [[pair_score(g, p.barcode) for p in plates] for g in grid_entries]
    if ng and nb:
        if _HAVE_SCIPY:
            cost = _np.array([[1.0 - s for s, _ in row] for row in scores])
            gi, bi = linear_sum_assignment(cost)
            assign = dict(zip(gi.tolist(), bi.tolist()))
        else:
            assign, taken = {}, set()
            ranked = sorted(
                ((scores[g][b][0], g, b) for g in range(ng) for b in range(nb)),
                reverse=True,
            )
            for _, g, b in ranked:
                if g not in assign and b not in taken:
                    assign[g] = b
                    taken.add(b)
    else:
        assign = {}

    out: list[dict] = []
    for g in range(ng):
        b = assign.get(g)
        if b is None:
            out.append(dict(grid_position=g + 1, grid_entry=grid_entries[g],
                            matched_position=None, barcode=None, score=0.0,
                            method="no_plate_available", accepted=False))
            continue
        score, method = scores[g][b]
        accepted = score >= min_score
        out.append(dict(
            grid_position=g + 1,
            grid_entry=grid_entries[g],
            matched_position=(b + 1) if accepted else None,
            barcode=plates[b].barcode,
            score=round(score, 3),
            method=method if accepted else f"below_threshold({method})",
            accepted=accepted,
        ))
    return out


def align_strings(grid_entries: list[str], candidates: list[str], min_score: float) -> list[dict]:
    """Fuzzy-match grid entries to arbitrary strings using the same pair scoring."""
    ng, nc = len(grid_entries), len(candidates)
    scores = [[pair_score(g, c) for c in candidates] for g in grid_entries]
    if ng and nc:
        if _HAVE_SCIPY:
            cost = _np.array([[1.0 - score for score, _ in row] for row in scores])
            gi, ci = linear_sum_assignment(cost)
            assign = dict(zip(gi.tolist(), ci.tolist()))
        else:
            assign, taken = {}, set()
            ranked = sorted(
                ((scores[g][c][0], g, c) for g in range(ng) for c in range(nc)),
                reverse=True,
            )
            for _, g, c in ranked:
                if g not in assign and c not in taken:
                    assign[g] = c
                    taken.add(c)
    else:
        assign = {}

    out: list[dict] = []
    for g in range(ng):
        c = assign.get(g)
        if c is None:
            out.append(dict(grid_position=g + 1, grid_entry=grid_entries[g],
                            matched_position=None, candidate=None, score=0.0,
                            method="no_candidate_available", accepted=False))
            continue
        score, method = scores[g][c]
        accepted = score >= min_score
        out.append(dict(
            grid_position=g + 1,
            grid_entry=grid_entries[g],
            matched_position=(c + 1) if accepted else None,
            candidate=candidates[c],
            score=round(score, 3),
            method=method if accepted else f"below_threshold({method})",
            accepted=accepted,
        ))
    return out


def describe_permutation(pairs: list[dict]) -> str:
    acc = [(p["grid_position"], p["matched_position"]) for p in pairs if p["accepted"]]
    if not acc:
        return "no confident matches"
    gp, mp = [a for a, _ in acc], [b for _, b in acc]
    if gp == mp:
        return "in order"
    if mp == sorted(mp, reverse=True):
        return "REVERSED"
    offsets = [b - a for a, b in acc]
    if len(set(offsets)) == 1:
        return f"SHIFTED by {offsets[0]:+d} throughout"
    if sorted(mp) == mp and len(set(offsets)) == 2:
        for k in range(1, len(offsets)):
            if len(set(offsets[:k])) == 1 and len(set(offsets[k:])) == 1:
                return f"SHIFTED by {offsets[k] - offsets[0]:+d} from grid position {acc[k][0]} onward"
    if sorted(mp) == mp:
        return "ascending but with gaps"
    return "REORDERED (see detail rows)"


@dataclass
class FolderResult:
    folder: str
    grid_file: str
    n_grid_entries: int
    data_files: list[str] = field(default_factory=list)
    csv_files: list[str] = field(default_factory=list)
    envision_files: list[str] = field(default_factory=list)
    styles: list[str] = field(default_factory=list)
    verdict: str = ""
    detail: str = ""
    n_plates: int = 0
    n_matched: int = 0
    order: str = ""
    filename_order: str = ""
    barcode_filename_match: str = ""
    blank_barcodes: int = 0
    dates_ok: bool | None = None
    rows: list[dict] = field(default_factory=list)


def audit_folder(folder: str, grid_path: str, min_score: float, strict_dates: bool) -> FolderResult:
    entries = parse_grid_num(grid_path)
    res = FolderResult(folder=folder, grid_file=os.path.basename(grid_path), n_grid_entries=len(entries))

    for name in sorted(os.listdir(folder)):
        full = os.path.join(folder, name)
        if not os.path.isfile(full) or name.lower() in GRID_FILENAMES or name.startswith("._"):
            continue
        if os.path.splitext(name)[1].lower() in DATA_EXTS:
            res.data_files.append(name)
        if name.lower().endswith((".csv", ".tsv")):
            res.csv_files.append(name)

    if not res.data_files:
        res.verdict, res.detail = "NO_DATA_FILES", "grid_num.txt with nothing beside it"
        return res
    if not res.csv_files:
        res.verdict, res.detail = "NO_CSV", "data present but not CSV: " + ", ".join(res.data_files[:6])
        return res

    kinds: dict[str, tuple[str, str, bool]] = {}
    for name in res.csv_files:
        kind, reason, has_grid = classify_csv(os.path.join(folder, name))
        kinds[name] = (kind, reason, has_grid)
        if kind in {"envision_old", "envision_new"}:
            res.envision_files.append(name)
            res.styles.append("new" if has_grid else "old")

    if not res.envision_files:
        res.verdict = "CSV_NOT_ENVISION"
        res.detail = "; ".join(f"{n}: {kinds[n][1]}" for n in res.csv_files[:4])
        return res

    per_file = {name: parse_envision_plates(os.path.join(folder, name)) for name in res.envision_files}

    # The downstream pipeline sorts EnVision filenames alphanumerically before
    # pairing them with grid_num. A filename/grid_num order check is meaningful
    # when each file contains exactly one plate, or when the number of grid_num
    # entries equals the number of EnVision files. For multi-plate files whose
    # filenames legitimately skip numbers, do not flag the filename sequence;
    # the plate-level count check below still catches a missing file/plate.
    if len(res.envision_files) > 1:
        filename_stems = [os.path.splitext(name)[0] for name in sorted(res.envision_files)]
        filename_pairs = align_strings(entries, filename_stems, min_score)
        filename_unmatched = [p["grid_entry"] for p in filename_pairs if not p["accepted"]]
        one_plate_per_file = all(len(per_file[name]) == 1 for name in res.envision_files)
        filename_order_check = one_plate_per_file or len(entries) == len(res.envision_files)
        if filename_order_check:
            res.filename_order = describe_permutation(filename_pairs)
        else:
            res.filename_order = "n/a (multi-plate files)"
    else:
        filename_pairs = []
        filename_unmatched = []
        filename_order_check = False

    # When an EnVision file contains exactly one plate, compare that plate's
    # barcode directly to the filename stem using the same fuzzy scoring logic.
    # This applies independently to each file, including folders with multiple
    # one-plate EnVision exports.
    barcode_filename_results = []
    for name in sorted(res.envision_files):
        file_plates = per_file[name]
        if len(file_plates) != 1 or not file_plates[0].barcode or file_plates[0].blank_barcode:
            continue
        plate = file_plates[0]
        score, method = pair_score(os.path.splitext(name)[0], plate.barcode)
        status = "PASS" if score >= min_score else "FAIL"
        barcode_filename_results.append(
            f"{name}: {status} {score:.2f} ({method})"
        )
    res.barcode_filename_match = "; ".join(barcode_filename_results) if barcode_filename_results else "n/a"

    # Multiple files can be complementary exports or repeated reads of the same plates.
    if len(per_file) > 1:
        sets = [{norm(p.barcode) for p in v if p.barcode} for v in per_file.values()]
        overlap = any(a & b for i, a in enumerate(sets) for b in sets[i + 1:])
        if overlap:
            best = max(per_file, key=lambda k: len(per_file[k]))
            plates = per_file[best]
            res.detail = f"{len(per_file)} EnVision CSVs share barcodes; compared against {best}. "
        else:
            plates = [p for name in res.envision_files for p in per_file[name]]
            res.detail = f"concatenated {len(per_file)} EnVision CSVs with disjoint barcodes. "
    else:
        plates = next(iter(per_file.values()))

    style_set = sorted(set(res.styles))
    style_label = style_set[0] if len(style_set) == 1 else "/".join(style_set)
    res.detail += f"style={style_label} (has_grid={'true' if 'new' in style_set else 'false'}). "

    res.n_plates = len(plates)
    res.blank_barcodes = sum(p.blank_barcode for p in plates)
    res.dates_ok = dates_monotonic(plates)

    if not plates:
        res.verdict, res.detail = "NO_PLATE_BLOCKS", res.detail + "no parseable Plate information blocks"
        return res
    all_barcodes_blank = res.blank_barcodes == len(plates)
    res.rows = align(entries, plates, min_score) if not all_barcodes_blank else [
        dict(grid_position=i + 1, grid_entry=e, matched_position=None,
             barcode="", score=0.0, method="barcode_blank", accepted=False)
        for i, e in enumerate(entries)
    ]

    # Add the filename associated with each grid position. For multiple EnVision
    # files this comes from the filename/grid fuzzy alignment above; for a single
    # file it is simply that filename. Keep rows in grid_num order for both
    # console and detail CSV output.
    if len(res.envision_files) == 1:
        row_filenames = [res.envision_files[0]] * len(entries)
    else:
        row_filenames = [p.get("candidate") or "" for p in filename_pairs]
    for row, filename in zip(res.rows, row_filenames):
        row["filename"] = filename

    res.n_matched = sum(r["accepted"] for r in res.rows)
    res.order = describe_permutation(res.rows)

    notes = []
    if len(entries) != len(plates):
        notes.append(f"count mismatch: {len(entries)} grid entries vs {len(plates)} plates")
    if res.blank_barcodes:
        notes.append(f"{res.blank_barcodes} blank barcode field(s)")
    unmatched = [r["grid_entry"] for r in res.rows if not r["accepted"]]
    if unmatched:
        notes.append(f"{len(unmatched)} unmatched: " + ", ".join(unmatched[:4]))
    weak = [r for r in res.rows if r["accepted"] and r["score"] < 0.9]
    if weak:
        notes.append(f"{len(weak)} fuzzy match(es)")
    if res.dates_ok is False:
        notes.append("measurement timestamps NOT monotonic with plate index")
    elif res.dates_ok is None:
        notes.append("timestamps unparseable, acquisition order unconfirmed")

    barcode_order_bad = (not all_barcodes_blank) and res.order != "in order"
    filename_order_bad = filename_order_check and res.filename_order != "in order"
    filename_count_bad = len(res.envision_files) > 1 and len(entries) != len(res.envision_files)
    barcode_filename_bad = "; FAIL " in res.barcode_filename_match or res.barcode_filename_match.startswith("FAIL ")

    # Barcode blanks mean the barcode-based checks have nothing to evaluate;
    # they must not suppress an independent filename/grid_num failure.
    if filename_order_bad or filename_count_bad or barcode_filename_bad or barcode_order_bad:
        res.verdict = "ORDER_MISMATCH"
    elif all_barcodes_blank:
        res.verdict = "BARCODES_BLANK"
    elif not unmatched and len(entries) == len(plates):
        res.verdict = "OK" if (res.dates_ok is not False or not strict_dates) else "REVIEW"
    elif res.n_matched == 0:
        res.verdict = "NO_MATCH"
    else:
        res.verdict = "REVIEW"

    if all_barcodes_blank:
        notes.insert(0, f"all {len(plates)} plate blocks have empty/unusable barcode fields")
    if barcode_order_bad and not all_barcodes_blank:
        notes.insert(0, f"grid_num order does not match barcode order: {res.order}")
    if filename_order_bad:
        notes.insert(0, f"grid_num order does not match filename order: {res.filename_order}")
    if barcode_filename_bad:
        notes.append(f"single-file barcode does not match filename: {res.barcode_filename_match}")
    if filename_count_bad:
        notes.append(f"count mismatch: {len(entries)} grid entries vs {len(res.envision_files)} EnVision CSV filenames")
    if filename_order_check and filename_unmatched:
        notes.append(f"{len(filename_unmatched)} filename-order unmatched: " + ", ".join(filename_unmatched[:4]))
    res.detail += "; ".join(notes) if notes else "all entries matched in order"
    return res


def find_grid_files(root: str, follow_links: bool = False) -> Iterable[tuple[str, str]]:
    for dirpath, dirnames, filenames in os.walk(root, followlinks=follow_links):
        dirnames[:] = [d for d in dirnames if not d.startswith((".", "__MACOSX"))]
        for name in filenames:
            if name.lower() in GRID_FILENAMES:
                yield dirpath, os.path.join(dirpath, name)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", help="directory to walk")
    ap.add_argument("-o", "--report", help="write folder-level summary CSV here")
    ap.add_argument("-d", "--detail", help="write per-plate match rows CSV here")
    ap.add_argument("--min-score", type=float, default=0.55)
    ap.add_argument("--strict-dates", action="store_true")
    ap.add_argument("--follow-links", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    if not _HAVE_SCIPY:
        print("note: SciPy not found, falling back to greedy matching", file=sys.stderr)

    results = []
    for folder, grid_path in find_grid_files(args.root, args.follow_links):
        try:
            results.append(audit_folder(folder, grid_path, args.min_score, args.strict_dates))
        except Exception as exc:
            r = FolderResult(folder=folder, grid_file=os.path.basename(grid_path), n_grid_entries=0)
            r.verdict, r.detail = "ERROR", f"{type(exc).__name__}: {exc}"
            results.append(r)

    if not results:
        print(f"No grid_num.txt found under {args.root}")
        return 0

    flagged = {"ORDER_MISMATCH", "NO_MATCH", "BARCODES_BLANK", "CSV_NOT_ENVISION",
               "NO_CSV", "NO_DATA_FILES", "NO_PLATE_BLOCKS", "REVIEW", "ERROR"}

    for r in results:
        mark = "ok  " if r.verdict == "OK" else "FLAG"
        print(f"[{mark}] {r.verdict:<18} {r.folder}")
        print(f"         grid entries: {r.n_grid_entries}  plates: {r.n_plates}  matched: {r.n_matched}  order: {r.order or 'n/a'}")
        if len(r.envision_files) > 1:
            print(f"         filename order: {r.filename_order or 'n/a'}")
        print(f"         barcode vs filename: {r.barcode_filename_match or 'n/a'}")
        print(f"         envision style: {('/'.join(sorted(set(r.styles)))) if r.styles else 'n/a'}")
        if r.envision_files:
            print(f"         envision CSVs: {', '.join(r.envision_files)}")
        elif r.csv_files:
            print(f"         CSVs: {', '.join(r.csv_files[:5])}")
        if r.detail:
            print(f"         {r.detail}")
        if r.rows and (args.verbose or r.verdict in flagged):
            for row in r.rows:
                flag = " " if row["accepted"] and row["grid_position"] == row["matched_position"] else "!"
                print(f"          {flag} #{row['grid_position']:<3} {row['grid_entry']:<28} -> file {str(row.get('filename') or '?'):<28} plate #{str(row['matched_position'] or '?'):<3} barcode {str(row['barcode'] or ''):<14} [{row['score']:.2f} {row['method']}]")
        print()

    counts: dict[str, int] = {}
    for r in results:
        counts[r.verdict] = counts.get(r.verdict, 0) + 1
    print("summary: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))

    if args.report:
        with open(args.report, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["folder", "grid_file", "verdict", "style", "order", "filename_order", "barcode_filename_match", "n_grid_entries",
                        "n_plates", "n_matched", "n_blank_barcodes", "timestamps_monotonic",
                        "envision_csvs", "all_csvs", "other_data_files", "detail"])
            for r in results:
                w.writerow([r.folder, r.grid_file, r.verdict, "/".join(sorted(set(r.styles))), r.order, r.filename_order, r.barcode_filename_match,
                            r.n_grid_entries, r.n_plates, r.n_matched, r.blank_barcodes,
                            "" if r.dates_ok is None else r.dates_ok,
                            "|".join(r.envision_files), "|".join(r.csv_files),
                            "|".join(f for f in r.data_files if f not in r.csv_files), r.detail])
        print(f"wrote {args.report}")

    if args.detail:
        with open(args.detail, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["folder", "grid_position", "grid_entry", "filename", "matched_plate_position",
                        "barcode", "score", "method", "accepted", "position_preserved"])
            for r in results:
                for row in r.rows:
                    w.writerow([r.folder, row["grid_position"], row["grid_entry"], row.get("filename", ""), row["matched_position"],
                                row["barcode"], row["score"], row["method"], row["accepted"],
                                row["accepted"] and row["grid_position"] == row["matched_position"]])
        print(f"wrote {args.detail}")

    return 1 if any(r.verdict in flagged for r in results) else 0


if __name__ == "__main__":
    sys.exit(main())