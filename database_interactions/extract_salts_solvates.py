#!/usr/bin/env python
"""
extract_salts_solvates.py

Split combinatorial registration SMILES (parent compound + salts/solvates
joined by ".") into a clean parent SMILES plus salt/solvate name+count
columns, matching fragments against reference SMARTS tables. Also handles
the residual multi-fragment "co-formulation" rows left after salt/solvate
stripping: reduces repeated-identical-fragment stoichiometry, collapses
genuine stereoisomer pairs into a single CXSMILES string, and classifies +
annotates true mixtures of distinct compounds.

Requires rdkit + pandas (use the `cheminf` conda env or any env with rdkit).

Reference tables
----------------
The salt and solvate reference SMARTS come from the HiTS tables INV_LOT_SALT
(ILSA_NAME, ILSA_SMARTS) and INV_LOT_SOLVATE (ILSO_NAME, ILSO_SMARTS). Read
them straight from the database by passing a SQLAlchemy engine, or from CSV
exports of those tables. If you pass both, the database is tried first and the
CSVs are used only if the database read fails or comes back empty.

The database read only ever runs SELECTs on those two tables, inside a
read-only transaction on MySQL, and rolls back afterwards. It never writes.
This module does not open connections itself: build the engine however you
normally do (e.g. ``safe_register.connect("QB3_HiTS")``) and pass it in.

Usage
-----
From a notebook, reading the reference tables from the database::

    import extract_salts_solvates as ess
    out = ess.run_pipeline("Registration Template.csv", engine=engine,
                           out_csv="extracted.csv", review_csv="review.csv",
                           header_row=0, smiles_col="SMILES")

From the command line with CSV exports (the old way, still supported)::

    python extract_salts_solvates.py \\
        --template "Registration Template.csv" \\
        --salts INV_LOT_SALT.csv \\
        --solvates INV_LOT_SOLVATE.csv \\
        --out extracted.csv \\
        --review review.csv \\
        --header-row 0 --smiles-col SMILES

From the command line reading from the database, give a SQLAlchemy URL in an
environment variable (never on the command line, where it lands in shell
history)::

    export HITS_DB_URL="mysql+pymysql://user:password@localhost:3306/QB3_HiTS"
    python extract_salts_solvates.py --template ... --db-url-env HITS_DB_URL

Output columns added to every row (original columns are preserved):
    original_smiles     -- the input SMILES, untouched
    extracted_smiles     -- parent compound only, salts/solvates stripped;
                             for the special cases below this holds either a
                             single reduced fragment or a single CXSMILES
                             string, never a "."-joined multi-fragment string
                             (except genuine true_mixture_diff_compounds rows)
    salt                  -- reference salt name(s), "; "-joined if several
    salt_num              -- matching count(s) for `salt` (may be fractional,
                             e.g. "0.5" for a hemi-salt)
    solvate               -- reference solvate name(s), "; "-joined
    solvate_num           -- matching count(s) for `solvate`
    extraction_flag       -- '' | 'parse_fallback' | 'all_components_salt_review'
                             | 'empty' | 'dative bonds processed to single bonds'
                             (semicolon-joined if more than one applies)
    mixture_class         -- '' | 'reduced_stoichiometry' | 'stereo_mixture_cxsmiles'
                             | 'true_mixture_diff_compounds'
                             (only set for rows where >1 parent fragment
                             remained after salt/solvate stripping)
    IC_STEREO             -- free-text stereochemistry/regiochemistry note,
                             populated only for stereo_mixture_cxsmiles and
                             regiochemical true_mixture rows (e.g.
                             "Racemic mixtures", "Diastereomeric mixture",
                             "Regioisomers", "Isomeric mixture")

True-mixture rows (mixture_class == 'true_mixture_diff_compounds') are, by
default, exploded into one row per distinct component (set
expand_mixtures=False / --no-expand-mixtures to keep them as a single
"."-joined row instead). Exploded rows get:
    <id_col>              -- suffixed '_part1', '_part2', ... (default
                             id_col is "CompoundIndex"; pass id_col= /
                             --id-col to use e.g. "lot_alias" instead)
    mixture_part_of       -- the original <id_col> value before suffixing
    mixture_part_count    -- how many times that exact fragment occurred in
                             the original row (e.g. 3 for Chromium
                             picolinate's picolinate ligand)
All other original columns are duplicated verbatim onto every part
(including extraction_flag, unchanged, on every part). Any salt/solvate
found in the original row is assigned to exactly one part -- the one whose
formal charge is complementary to the counter-ion/co-crystal partner -- and
cleared ('' / 0) on the others; if charge complementarity can't decide, it's
assigned to _part1 only.
"""
import argparse
import os
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, FindMolChiralCenters, StereoGroupType

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# Repairs for malformed SMARTS in older reference salt files (ring-closure
# digits inside atom brackets, a stray trailing digit/slash). Newer reference
# files may already ship these fixed at the source -- this dict is a no-op
# for any name/SMARTS pair that's already valid.
# ---------------------------------------------------------------------------
SMARTS_FIXES = {
    "formate":              "[O;D1]C=[O;D1]",
    "salicylate":           "c1([O;D1])c(C(=[O;D1])[O;D1])cccc1",
    "sulfosalicylate":      "c1([O;D1])c(C(=[O;D1])[O;D1])cc(S([O;D1])(=[O;D1])=[O;D1])cc1",
    "naphthalenesulfonate": "c1(S(=[O;D1])(=[O;D1])[O;D1])cccc2c(S(=[O;D1])(=[O;D1])[O;D1])cccc12",
    "triethylamine":        "[C;H3][C;H2][N;D3]([C;H2][C;H3])[C;H2][C;H3]",
    "maleate":              "[O;D1]C(=[O;D1])[C;H1]=[C;H1]C(=[O;D1])[O;D1]",
}


# ---------------------------------------------------------------------------
# Reference table building + per-fragment classification
# ---------------------------------------------------------------------------
def ss_has_dative_bonds(mol):
    """True if mol contains any RDKit DATIVE bond (see ss_strip_dative_bonds)."""
    return any(b.GetBondType() == Chem.BondType.DATIVE for b in mol.GetBonds())


def ss_strip_dative_bonds(mol):
    """Convert RDKit DATIVE bonds to ordinary SINGLE bonds on a copy, so
    downstream MolToSmiles/MolToCXSmiles never emit RDKit's non-standard
    '->'/'<-' dative-bond arrow notation (RDKit auto-perceives certain
    metal-donor bonds, e.g. macrocycle N/S -> transition metal, as dative
    on parse/sanitize, and always re-emits them with arrows -- most other
    SMILES parsers, and downstream registration/mass-balancing tools,
    reject that token). Preserves formula, formal charge, and connectivity;
    only changes how the metal-ligand bond is notated.
    """
    if not ss_has_dative_bonds(mol):
        return mol
    rw = Chem.RWMol(mol)
    for b in rw.GetBonds():
        if b.GetBondType() == Chem.BondType.DATIVE:
            b.SetBondType(Chem.BondType.SINGLE)
    return rw.GetMol()


DATIVE_FLAG_TEXT = "dative bonds processed to single bonds"


def ss_is_bare_h_ion(mol):
    """True if mol is a single hydrogen atom carrying a nonzero formal
    charge -- a bare [H+]/[H-] fragment. These show up in vendor SMILES as
    protonation-state bookkeeping (re-neutralizing an ionized macrocycle
    arm) rather than a real salt-former, so ss_process_row drops them
    entirely instead of matching them against the salt/solvate reference
    tables or carrying them through as a mixture component; the parent
    fragment keeps its native formal charge."""
    return (mol.GetNumAtoms() == 1
            and mol.GetAtomWithIdx(0).GetAtomicNum() == 1
            and mol.GetAtomWithIdx(0).GetFormalCharge() != 0)


BARE_H_ION_FLAG_TEXT = "stray H+/H- fragment removed (charge retained on parent)"


def ss_build_refs(df, name_col, smarts_col, fixes):
    refs = []
    for _, r in df.iterrows():
        val = r[name_col]
        nm = "" if val is None or (isinstance(val, float) and val != val) else str(val).strip()
        if not nm or nm == "nan":
            continue
        sval = r[smarts_col]
        sma = "" if sval is None or (isinstance(sval, float) and sval != sval) else str(sval).strip()
        if nm in fixes:
            sma = fixes[nm]
        q = Chem.MolFromSmarts(sma) if sma else None
        refs.append(dict(name=nm, query=q, multi=("." in sma), qok=(q is not None)))
    return refs


def ss_frag_full_matches(frag, refs):
    """A reference matches a fragment only if it covers every heavy atom
    (equal atom count + full substructure match) -- partial matches (e.g. a
    hydroxyl SMARTS hitting one -OH of a larger parent) are excluded."""
    n = frag.GetNumAtoms()
    hits = []
    for r in refs:
        q = r["query"]
        if q is None or r["multi"] or q.GetNumAtoms() != n:
            continue
        m = frag.GetSubstructMatch(q)
        if m and len(m) == n:
            hits.append(r)
    return hits


def ss_classify_frag(frag, salt_refs, solv_refs):
    """Classify one fragment as parent / salt / solvate.

    Ambiguity rule: a fragment matching BOTH a salt and a solvate reference
    (e.g. ammonia/ammonium, ethylenediamine, water) is resolved by charge --
    charged -> salt (counter-ion), neutral -> solvate (co-crystal solvent).
    """
    s = ss_frag_full_matches(frag, salt_refs)
    v = ss_frag_full_matches(frag, solv_refs)
    if not s and not v:
        return ("parent", None)
    charged = any(a.GetFormalCharge() != 0 for a in frag.GetAtoms())
    if s and v:
        pool, role = (s, "salt") if charged else (v, "solvate")
    elif s:
        pool, role = s, "salt"
    else:
        pool, role = v, "solvate"
    return (role, pool[0]["name"])


def ss_process_row(smi, salt_refs, solv_refs):
    """Split one SMILES string into parent fragment(s) + salt/solvate hits.
    Also tracks salt_charge / solvate_charge -- the formal charge of the
    single matched salt/solvate fragment, when there's exactly one distinct
    match (None when there are zero or multiple distinct names) -- used by
    ss_split_mixture_parts to decide which true-mixture component a salt or
    solvate counter-ion/co-crystal partner belongs with.

    Bare [H+]/[H-] fragments (a single charged hydrogen atom) are dropped
    entirely before classification -- these are vendor protonation-state
    bookkeeping (re-neutralizing an ionized macrocycle arm), not a real
    salt-former or solvate, so they're never matched against the reference
    tables and never appear as a mixture component; the parent fragment
    keeps its native formal charge. Dropping one sets
    BARE_H_ION_FLAG_TEXT on the row."""
    if not isinstance(smi, str) or smi.strip() == "":
        return dict(extracted_smiles=smi, salt="", salt_num=0, solvate="", solvate_num=0, flag="empty",
                    salt_charge=None, solvate_charge=None)
    mol = Chem.MolFromSmiles(smi)
    used_fallback = False
    if mol is not None:
        pieces = [("mol", f) for f in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)]
    else:
        used_fallback = True
        pieces = []
        for part in smi.split("."):
            fm = Chem.MolFromSmiles(part)
            pieces.append(("mol", fm) if fm is not None else ("txt", part))

    dropped_h_ion = False
    filtered_pieces = []
    for kind, obj in pieces:
        if kind == "mol" and obj is not None and ss_is_bare_h_ion(obj):
            dropped_h_ion = True
            continue
        filtered_pieces.append((kind, obj))
    pieces = filtered_pieces

    roles = []
    for kind, obj in pieces:
        if kind == "mol" and obj is not None:
            role, name = ss_classify_frag(obj, salt_refs, solv_refs)
            roles.append((obj, role, name, "mol"))
        else:
            roles.append((obj, "parent", None, "txt"))
    flag = "parse_fallback" if used_fallback else ""
    all_removed = False
    if roles and not any(t[1] == "parent" for t in roles):
        all_removed = True
        big = max(roles, key=lambda t: t[0].GetNumAtoms() if t[3] == "mol" else 0)
        roles = [(t[0], "parent", None, t[3]) if t is big else t for t in roles]

    dative_found = any(t[3] == "mol" and ss_has_dative_bonds(t[0]) for t in roles)

    def rep(t):
        return Chem.MolToSmiles(ss_strip_dative_bonds(t[0])) if t[3] == "mol" else t[0]

    extracted = ".".join(rep(t) for t in roles if t[1] == "parent")
    salt_order, salt_ct = [], Counter()
    solv_order, solv_ct = [], Counter()
    salt_charge_by_name, solv_charge_by_name = {}, {}
    for obj, role, name, _ in roles:
        if role == "salt":
            if name not in salt_ct:
                salt_order.append(name)
                salt_charge_by_name[name] = Chem.GetFormalCharge(obj)
            salt_ct[name] += 1
        elif role == "solvate":
            if name not in solv_ct:
                solv_order.append(name)
                solv_charge_by_name[name] = Chem.GetFormalCharge(obj)
            solv_ct[name] += 1
    if all_removed and (salt_order or solv_order):
        flag = (flag + ";" if flag else "") + "all_components_salt_review"
    if dative_found:
        flag = (flag + ";" if flag else "") + DATIVE_FLAG_TEXT
    if dropped_h_ion:
        flag = (flag + ";" if flag else "") + BARE_H_ION_FLAG_TEXT
    salt_charge = salt_charge_by_name[salt_order[0]] if len(salt_order) == 1 else None
    solv_charge = solv_charge_by_name[solv_order[0]] if len(solv_order) == 1 else None
    return dict(
        extracted_smiles=extracted,
        salt="; ".join(salt_order),
        salt_num="; ".join(str(salt_ct[n]) for n in salt_order) if salt_order else 0,
        solvate="; ".join(solv_order),
        solvate_num="; ".join(str(solv_ct[n]) for n in solv_order) if solv_order else 0,
        flag=flag,
        salt_charge=salt_charge,
        solvate_charge=solv_charge,
    )


def ss_varying_centers(molA, molB):
    """Atom-index-aligned diff of defined stereocenters between two mols
    that share the same (achiral) connectivity. Returns the list of atom
    indices in molA whose CIP code differs from the matched atom in molB,
    or None if the achiral skeletons don't align 1:1."""
    a_noster = Chem.Mol(molA)
    Chem.RemoveStereochemistry(a_noster)
    match = molB.GetSubstructMatch(a_noster, useChirality=False)
    if not match or len(match) != molA.GetNumAtoms():
        return None
    Chem.AssignStereochemistry(molA, cleanIt=True, force=True)
    Chem.AssignStereochemistry(molB, cleanIt=True, force=True)
    varying = []
    for i, atom in enumerate(molA.GetAtoms()):
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        atomB = molB.GetAtomWithIdx(match[i])
        if atom.GetPropsAsDict().get("_CIPCode") != atomB.GetPropsAsDict().get("_CIPCode"):
            varying.append(i)
    return varying


def ss_reduce_and_classify_fragments(extracted_smiles, salt_num_str, solv_num_str):
    """Post-process one row's `extracted_smiles` (already salt/solvate-
    stripped). Returns dict: extracted_smiles, salt_num, solvate_num,
    mixture_class, IC_STEREO. No-ops (mixture_class == '') for single-
    fragment rows or anything that doesn't parse."""
    out = dict(extracted_smiles=extracted_smiles, salt_num=salt_num_str, solvate_num=solv_num_str,
               mixture_class="", IC_STEREO="")
    if not isinstance(extracted_smiles, str) or "." not in extracted_smiles:
        return out
    frags = extracted_smiles.split(".")
    mols = [Chem.MolFromSmiles(f) for f in frags]
    if any(m is None for m in mols):
        return out  # leave unparseable fragments alone

    canos = [Chem.MolToSmiles(m) for m in mols]
    canos_clean = [Chem.MolToSmiles(ss_strip_dative_bonds(m)) for m in mols]
    n = len(frags)

    # Case 1: every fragment identical, incl. stereo -> repeated stoichiometry
    if len(set(canos)) == 1:
        def reduce_counts(num_str):
            if num_str in (0, "0", ""):
                return num_str
            parts = [float(x) for x in str(num_str).split(";")]
            new_parts = [p / n for p in parts]
            return "; ".join(str(int(p)) if float(p).is_integer() else str(p) for p in new_parts)
        out["extracted_smiles"] = canos_clean[0]
        out["salt_num"] = reduce_counts(salt_num_str)
        out["solvate_num"] = reduce_counts(solv_num_str)
        out["mixture_class"] = "reduced_stoichiometry"
        return out

    formulas = [rdMolDescriptors.CalcMolFormula(m) for m in mols]
    if len(set(formulas)) != 1 or n != 2:
        # different compounds entirely, or a >2-fragment co-formulation --
        # leave untouched, just classify as a true mixture.
        out["mixture_class"] = "true_mixture_diff_compounds"
        return out

    achiral = [Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(m, isomericSmiles=False))) for m in mols]
    if achiral[0] != achiral[1]:
        # same formula, different connectivity -> regiochemical/constitutional
        # isomer mixture, not a stereo difference.
        out["mixture_class"] = "true_mixture_diff_compounds"
        out["IC_STEREO"] = "Regioisomers"
        return out

    varying = ss_varying_centers(mols[0], mols[1])
    if not varying:
        out["mixture_class"] = "true_mixture_diff_compounds"
        out["IC_STEREO"] = "Isomeric mixture"
        return out

    total_centers = len(FindMolChiralCenters(mols[0], includeUnassigned=True, useLegacyImplementation=False))
    rwm = Chem.RWMol(mols[0])
    sg = Chem.CreateStereoGroup(StereoGroupType.STEREO_AND, rwm, varying)
    rwm.SetStereoGroups([sg])
    m2 = rwm.GetMol()
    Chem.AssignStereochemistry(m2, cleanIt=True, force=True)
    m2_clean = ss_strip_dative_bonds(m2)
    if m2_clean is not m2:
        m2_clean.SetStereoGroups(m2.GetStereoGroups())
    out["extracted_smiles"] = Chem.MolToCXSmiles(m2_clean)
    out["mixture_class"] = "stereo_mixture_cxsmiles"
    out["IC_STEREO"] = "Racemic mixtures" if len(varying) == total_centers else "Diastereomeric mixture"
    return out


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------
def ss_split_mixture_parts(extracted_smiles, salt_charge, solvate_charge):
    """Given a '.'-joined extracted_smiles for a true_mixture_diff_compounds
    row, group fragments by canonical SMILES (preserving first-seen order,
    with dative bonds re-stripped after the parse round-trip, since RDKit
    can re-perceive them on re-parsing certain metal-donor ring fragments)
    and decide which group the salt / solvate should attach to, using a
    charge-complementarity heuristic: if exactly one group carries the
    opposite sign of charge from the salt/solvate's own (single,
    unambiguous) charge, assign there; otherwise default to the first group
    (part1). Returns None if the string doesn't split into >=2 distinct
    groups or any fragment fails to parse (e.g. residual parse_fallback
    text)."""
    if not isinstance(extracted_smiles, str) or "." not in extracted_smiles:
        return None
    frags = extracted_smiles.split(".")
    mols = [Chem.MolFromSmiles(f) for f in frags]
    if any(m is None for m in mols):
        return None
    mols = [ss_strip_dative_bonds(m) for m in mols]
    order, counts, charge_by_group = [], Counter(), {}
    for m in mols:
        c = Chem.MolToSmiles(m)
        if c not in charge_by_group:
            order.append(c)
            charge_by_group[c] = Chem.GetFormalCharge(m)
        counts[c] += 1
    if len(order) < 2:
        return None

    def pick_target(counter_charge):
        if counter_charge is None or counter_charge == 0:
            return 0
        opposite = [i for i, g in enumerate(order)
                    if charge_by_group[g] != 0 and (charge_by_group[g] > 0) != (counter_charge > 0)]
        return opposite[0] if len(opposite) == 1 else 0

    return dict(order=order, counts=counts, charge_by_group=charge_by_group,
                salt_target=pick_target(salt_charge), solvate_target=pick_target(solvate_charge))


def ss_expand_true_mixtures(df, id_col, id_fallback_col=None):
    """Explode every row with mixture_class == 'true_mixture_diff_compounds'
    into one row per distinct component, appending '_part1', '_part2', ... to
    `id_col` (or, when the `id_col` value for that row is missing/blank and
    `id_fallback_col` is given, appending the suffix to `id_fallback_col`
    instead -- e.g. id_col="lot_alias", id_fallback_col="chemist_ref" for
    rows that only have a chemist reference and no lot alias yet). All other
    columns are duplicated verbatim (including extraction_flag, unchanged,
    on every part) except extracted_smiles (set to that part's own
    re-stripped SMILES) and salt/salt_num/solvate/solvate_num, which are kept
    only on the part chosen by ss_split_mixture_parts (via charge
    complementarity) and cleared ('' / 0) on the others. Adds
    `mixture_part_count` (occurrence count of that fragment in the original
    row) and `mixture_part_of` (original id before suffixing, resolved
    through the same fallback). Rows that aren't true_mixture_diff_compounds,
    or that can't be split (parse failure, <2 distinct groups), pass through
    unchanged."""
    import pandas as pd

    def _blank(v):
        return v is None or (isinstance(v, float) and pd.isna(v)) or (isinstance(v, str) and v.strip() == "")

    out_rows = []
    for _, row in df.iterrows():
        rowd = row.to_dict()
        if rowd.get("mixture_class") != "true_mixture_diff_compounds":
            out_rows.append(rowd)
            continue
        info = ss_split_mixture_parts(rowd.get("extracted_smiles"), rowd.get("salt_charge"), rowd.get("solvate_charge"))
        if info is None:
            out_rows.append(rowd)
            continue
        orig_id = rowd.get(id_col)
        used_col = id_col
        if _blank(orig_id) and id_fallback_col and id_fallback_col in rowd:
            fb = rowd.get(id_fallback_col)
            if not _blank(fb):
                orig_id, used_col = fb, id_fallback_col
        for i, cano in enumerate(info["order"]):
            new_row = dict(rowd)
            new_row[used_col] = f"{orig_id}_part{i + 1}"
            new_row["mixture_part_of"] = orig_id
            new_row["extracted_smiles"] = cano
            new_row["mixture_part_count"] = info["counts"][cano]
            if i != info["salt_target"]:
                new_row["salt"] = ""
                new_row["salt_num"] = 0
            if i != info["solvate_target"]:
                new_row["solvate"] = ""
                new_row["solvate_num"] = 0
            out_rows.append(new_row)
    return pd.DataFrame(out_rows)


# ---------------------------------------------------------------------------
# Reference tables: database (read-only) or CSV exports
# ---------------------------------------------------------------------------
SALT_TABLE = "INV_LOT_SALT"
SOLVATE_TABLE = "INV_LOT_SOLVATE"


def read_reference_table(engine, table_name, name_col, smarts_col):
    """Read the name + SMARTS columns of one reference table, read-only.

    Builds a SELECT of just those two columns (identifiers are quoted by
    SQLAlchemy, so odd names can't inject SQL). On MySQL it runs inside
    ``START TRANSACTION READ ONLY``, so the server itself refuses any write;
    on every backend the transaction is rolled back afterwards.
    `engine` may be an Engine or an already-open Connection.
    """
    import pandas as pd
    from sqlalchemy import column, select, table
    from sqlalchemy.engine import Connection

    query = select(column(name_col), column(smarts_col)).select_from(table(table_name))

    def _read(conn):
        if conn.dialect.name == "mysql":
            conn = conn.execution_options(isolation_level="AUTOCOMMIT")
            conn.exec_driver_sql("START TRANSACTION READ ONLY")
            try:
                rows = conn.execute(query).fetchall()
            finally:
                conn.exec_driver_sql("ROLLBACK")
        else:
            with conn.begin() as tx:
                rows = conn.execute(query).fetchall()
                tx.rollback()
        return pd.DataFrame(rows, columns=[name_col, smarts_col])

    if isinstance(engine, Connection):
        return _read(engine)
    with engine.connect() as conn:
        return _read(conn)


def load_reference_tables(engine=None, salts_csv=None, solvates_csv=None,
                          salt_name_col="ILSA_NAME", salt_smarts_col="ILSA_SMARTS",
                          solv_name_col="ILSO_NAME", solv_smarts_col="ILSO_SMARTS"):
    """Return (salt_df, solvate_df) reference tables.

    With `engine`, reads INV_LOT_SALT and INV_LOT_SOLVATE from the database
    (see `read_reference_table`). If that read fails or returns an empty
    table and both CSV paths are given, falls back to the CSVs with a
    warning; without CSVs the error is raised. Without `engine`, both CSV
    paths are required. Prints which source was used.
    """
    import pandas as pd

    have_csv = salts_csv is not None and solvates_csv is not None
    if engine is None and not have_csv:
        raise ValueError("Pass engine=... to read the reference tables from the database, "
                         "or both salts_csv and solvates_csv.")

    if engine is not None:
        try:
            salt_df = read_reference_table(engine, SALT_TABLE, salt_name_col, salt_smarts_col)
            solv_df = read_reference_table(engine, SOLVATE_TABLE, solv_name_col, solv_smarts_col)
            empty = [t for t, d in ((SALT_TABLE, salt_df), (SOLVATE_TABLE, solv_df)) if d.empty]
            if empty:
                raise ValueError(f"reference table(s) came back empty: {', '.join(empty)}")
        except Exception as exc:
            if not have_csv:
                raise
            print(f"WARNING: could not read reference tables from the database ({exc}); "
                  f"using CSVs {salts_csv} and {solvates_csv} instead.")
        else:
            print(f"Reference tables read from the database: {len(salt_df)} salts, "
                  f"{len(solv_df)} solvates.")
            return salt_df, solv_df

    salt_df = pd.read_csv(salts_csv)
    solv_df = pd.read_csv(solvates_csv)
    print(f"Reference tables read from CSV: {len(salt_df)} salts, {len(solv_df)} solvates.")
    return salt_df, solv_df


def run_pipeline(template_csv, salts_csv=None, solvates_csv=None, out_csv=None, review_csv=None,
                  header_row=1, smiles_col="smiles",
                  salt_name_col="ILSA_NAME", salt_smarts_col="ILSA_SMARTS",
                  solv_name_col="ILSO_NAME", solv_smarts_col="ILSO_SMARTS",
                  id_col="CompoundIndex", id_fallback_col=None, expand_mixtures=True,
                  engine=None):
    """Extract salts/solvates from combinatorial SMILES, then reduce repeated
    stoichiometry, collapse stereoisomer pairs into CXSMILES, and (if
    `expand_mixtures` and `id_col` is present in the data) explode true
    mixtures of distinct compounds into one row per component with `id_col`
    (falling back to `id_fallback_col` for rows where `id_col` is
    missing/blank, e.g. id_col="lot_alias", id_fallback_col="chemist_ref")
    suffixed '_part1'/'_part2'/... .

    Requires rdkit + pandas. Writes out_csv (full result) and review_csv
    (rows that still contain "." or carry a non-empty extraction_flag or
    mixture_class) if given. Returns the full result DataFrame with original
    columns preserved plus original_smiles, extracted_smiles, salt,
    salt_num, solvate, solvate_num, extraction_flag, mixture_class,
    IC_STEREO, and (for exploded rows) mixture_part_of / mixture_part_count.

    Reference tables: pass `engine` (a SQLAlchemy engine or connection) to read
    INV_LOT_SALT / INV_LOT_SOLVATE from the database read-only, and/or
    `salts_csv` / `solvates_csv` as CSV exports. See `load_reference_tables`.
    """
    import pandas as pd
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    raw = pd.read_csv(template_csv, header=header_row, dtype=str, keep_default_na=False)
    salt_df, solv_df = load_reference_tables(
        engine=engine, salts_csv=salts_csv, solvates_csv=solvates_csv,
        salt_name_col=salt_name_col, salt_smarts_col=salt_smarts_col,
        solv_name_col=solv_name_col, solv_smarts_col=solv_smarts_col,
    )
    salt_refs = ss_build_refs(salt_df, salt_name_col, salt_smarts_col, SMARTS_FIXES)
    solv_refs = ss_build_refs(solv_df, solv_name_col, solv_smarts_col, {})
    bad = [r["name"] for r in salt_refs + solv_refs if not r["qok"]]
    if bad:
        print("WARNING: reference SMARTS still not compiling:", bad)

    res = raw[smiles_col].apply(lambda s: ss_process_row(s, salt_refs, solv_refs)).apply(pd.Series)
    final = raw.copy()
    final["original_smiles"] = raw[smiles_col]
    for c in ["extracted_smiles", "salt", "salt_num", "solvate", "solvate_num", "salt_charge", "solvate_charge"]:
        final[c] = res[c]
    final["extraction_flag"] = res["flag"]

    # post-process residual multi-fragment parent SMILES
    post = final.apply(
        lambda row: ss_reduce_and_classify_fragments(row["extracted_smiles"], row["salt_num"], row["solvate_num"]),
        axis=1,
    ).apply(pd.Series)
    final["extracted_smiles"] = post["extracted_smiles"]
    final["salt_num"] = post["salt_num"]
    final["solvate_num"] = post["solvate_num"]
    final["mixture_class"] = post["mixture_class"]
    final["IC_STEREO"] = post["IC_STEREO"]

    if expand_mixtures and id_col in final.columns:
        final = ss_expand_true_mixtures(final, id_col, id_fallback_col=id_fallback_col)

    final = final.drop(columns=["salt_charge", "solvate_charge"])

    if out_csv:
        final.to_csv(out_csv, index=False)
    if review_csv:
        has_dot = final["extracted_smiles"].astype(str).str.contains(".", regex=False)
        review = final[has_dot | (final["extraction_flag"] != "") | (final["mixture_class"] != "")]
        review.to_csv(review_csv, index=False)

    print(f"rows: {len(final)} | with salt: {(final['salt'] != '').sum()} | "
          f"with solvate: {(final['solvate'] != '').sum()} | "
          f"flags: {dict(final['extraction_flag'].value_counts())} | "
          f"mixture_class: {dict(final['mixture_class'].value_counts())}")
    return final



def _cli():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--template", required=True, help="registration template CSV")
    p.add_argument("--salts", default=None, help="salt reference CSV (ILSA_* columns); "
                    "with --db-url-env, used only as a fallback")
    p.add_argument("--solvates", default=None, help="solvate reference CSV (ILSO_* columns); "
                    "with --db-url-env, used only as a fallback")
    p.add_argument("--db-url-env", default=None,
                    help="name of an environment variable holding a SQLAlchemy URL; if given, the "
                         "reference tables are read from that database (read-only)")
    p.add_argument("--out", default=None, help="output CSV (full result)")
    p.add_argument("--review", default=None, help="output CSV (flagged/dotted rows only)")
    p.add_argument("--header-row", type=int, default=1, help="header row index for --template (0-based)")
    p.add_argument("--smiles-col", default="smiles", help="SMILES column name in --template")
    p.add_argument("--salt-name-col", default="ILSA_NAME")
    p.add_argument("--salt-smarts-col", default="ILSA_SMARTS")
    p.add_argument("--solv-name-col", default="ILSO_NAME")
    p.add_argument("--solv-smarts-col", default="ILSO_SMARTS")
    p.add_argument("--id-col", default="CompoundIndex",
                    help="row-identifier column to suffix _part1/_part2/... when exploding true mixtures")
    p.add_argument("--id-fallback-col", default=None,
                    help="fallback row-identifier column used (instead of --id-col) for any row where "
                         "--id-col is missing/blank, e.g. --id-col lot_alias --id-fallback-col chemist_ref")
    p.add_argument("--no-expand-mixtures", action="store_true",
                    help="keep true-mixture rows as a single '.'-joined row instead of exploding into parts")
    args = p.parse_args()
    engine = None
    if args.db_url_env:
        url = os.environ.get(args.db_url_env)
        if not url:
            p.error(f"environment variable {args.db_url_env} is not set")
        from sqlalchemy import create_engine
        engine = create_engine(url)
    elif not (args.salts and args.solvates):
        p.error("give --salts and --solvates, or --db-url-env")
    run_pipeline(
        template_csv=args.template, salts_csv=args.salts, solvates_csv=args.solvates,
        out_csv=args.out, review_csv=args.review, header_row=args.header_row,
        smiles_col=args.smiles_col, salt_name_col=args.salt_name_col,
        salt_smarts_col=args.salt_smarts_col, solv_name_col=args.solv_name_col,
        solv_smarts_col=args.solv_smarts_col, id_col=args.id_col,
        id_fallback_col=args.id_fallback_col,
        expand_mixtures=not args.no_expand_mixtures, engine=engine,
    )


if __name__ == "__main__":
    _cli()
