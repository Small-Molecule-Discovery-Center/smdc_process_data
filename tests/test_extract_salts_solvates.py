"""Small end-to-end tests for extract_salts_solvates.run_pipeline.

Uses tiny made-up reference tables, so no real HiTS data is needed.
Run with:  pytest   (from the repo root)
"""

import pandas as pd
import pytest
from sqlalchemy import create_engine, event

pytest.importorskip("rdkit")  # needs an env with rdkit, e.g. `cheminf`

import extract_salts_solvates as ess

SALTS = pd.DataFrame({
    "ILSA_NAME": ["hydrochloride", "sodium"],
    "ILSA_SMARTS": ["[Cl&D0,Cl&H1]", "[Na&D0]"],
})
SOLVATES = pd.DataFrame({
    "ILSO_NAME": ["water"],
    "ILSO_SMARTS": ["[O&H2&D0,O&D0]"],
})
ROWS = [
    ("plain",   "CCO"),                                  # nothing to strip
    ("salt",    "CCN.Cl"),                               # amine HCl
    ("solvate", "CC(=O)O.O"),                            # monohydrate
    ("both",    "CCN.Cl.O.O"),                           # HCl dihydrate
    ("hemi",    "CCN.CCN.Cl"),                           # repeated parent -> 0.5 salt
    ("stereo",  "C[C@H](N)C(=O)O.C[C@@H](N)C(=O)O"),     # enantiomer pair
    ("mixture", "c1ccccc1.CCCCCC"),                      # two different compounds
    ("ionmix",  "CC(=O)[O-].[Na+].c1ccncc1"),            # mixture with a counter-ion
]


@pytest.fixture
def result(tmp_path):
    template = tmp_path / "template.csv"
    salts = tmp_path / "salts.csv"
    solvates = tmp_path / "solvates.csv"
    pd.DataFrame(ROWS, columns=["CompoundIndex", "smiles"]).to_csv(template, index=False)
    SALTS.to_csv(salts, index=False)
    SOLVATES.to_csv(solvates, index=False)
    out = ess.run_pipeline(str(template), str(salts), str(solvates), header_row=0)
    return out.set_index("CompoundIndex")


def test_plain_row_untouched(result):
    r = result.loc["plain"]
    assert r.extracted_smiles == "CCO"
    assert r.salt == "" and r.solvate == "" and r.mixture_class == ""


def test_salt_extracted(result):
    r = result.loc["salt"]
    assert r.extracted_smiles == "CCN"
    assert r.salt == "hydrochloride" and float(r.salt_num) == 1


def test_solvate_extracted(result):
    r = result.loc["solvate"]
    assert r.extracted_smiles == "CC(=O)O"
    assert r.solvate == "water" and float(r.solvate_num) == 1


def test_salt_and_solvate_counts(result):
    r = result.loc["both"]
    assert r.extracted_smiles == "CCN"
    assert float(r.salt_num) == 1 and float(r.solvate_num) == 2


def test_repeated_parent_gives_hemi_salt(result):
    r = result.loc["hemi"]
    assert r.mixture_class == "reduced_stoichiometry"
    assert r.extracted_smiles == "CCN" and float(r.salt_num) == 0.5


def test_enantiomer_pair_becomes_one_cxsmiles(result):
    r = result.loc["stereo"]
    assert r.mixture_class == "stereo_mixture_cxsmiles"
    assert "." not in r.extracted_smiles
    assert r.IC_STEREO == "Racemic mixtures"


def test_true_mixture_is_split_into_parts(result):
    parts = result[result["mixture_part_of"] == "mixture"]
    assert sorted(parts.index) == ["mixture_part1", "mixture_part2"]
    assert set(parts.extracted_smiles) == {"c1ccccc1", "CCCCCC"}


def test_counter_ion_goes_to_one_part_only(result):
    parts = result[result["mixture_part_of"] == "ionmix"]
    assert len(parts) == 2
    assert (parts.salt == "sodium").sum() == 1
    assert parts.loc[parts.salt == "sodium", "extracted_smiles"].item() == "CC(=O)[O-]"


# ---------------------------------------------------------------------------
# Reading the reference tables from a database (SQLite stands in for HiTS)
# ---------------------------------------------------------------------------
@pytest.fixture
def ref_engine():
    engine = create_engine("sqlite://")
    SALTS.assign(ILSA_ID=range(len(SALTS))).to_sql("INV_LOT_SALT", engine, index=False)
    SOLVATES.assign(ILSO_ID=range(len(SOLVATES))).to_sql("INV_LOT_SOLVATE", engine, index=False)
    return engine


@pytest.fixture
def template_csv(tmp_path):
    path = tmp_path / "template.csv"
    pd.DataFrame(ROWS, columns=["CompoundIndex", "smiles"]).to_csv(path, index=False)
    return str(path)


def test_db_gives_same_result_as_csv(result, ref_engine, template_csv):
    from_db = ess.run_pipeline(template_csv, engine=ref_engine, header_row=0).set_index("CompoundIndex")
    pd.testing.assert_frame_equal(from_db, result)


def test_db_read_only_runs_selects(ref_engine, template_csv):
    seen = []
    event.listen(ref_engine, "before_cursor_execute",
                 lambda conn, cur, stmt, *a: seen.append(stmt.strip().split()[0].upper()))
    ess.run_pipeline(template_csv, engine=ref_engine, header_row=0)
    assert seen and set(seen) == {"SELECT"}
    # the reference tables are unchanged afterwards
    assert pd.read_sql("SELECT COUNT(*) AS n FROM INV_LOT_SALT", ref_engine).n[0] == len(SALTS)


def test_db_selects_only_name_and_smarts(ref_engine):
    df = ess.read_reference_table(ref_engine, "INV_LOT_SALT", "ILSA_NAME", "ILSA_SMARTS")
    assert list(df.columns) == ["ILSA_NAME", "ILSA_SMARTS"]
    assert sorted(df.ILSA_NAME) == sorted(SALTS.ILSA_NAME)


def test_falls_back_to_csv_when_db_read_fails(result, tmp_path, template_csv, capsys):
    salts, solvates = tmp_path / "salts.csv", tmp_path / "solvates.csv"
    SALTS.to_csv(salts, index=False)
    SOLVATES.to_csv(solvates, index=False)
    empty_db = create_engine("sqlite://")          # tables don't exist
    out = ess.run_pipeline(template_csv, str(salts), str(solvates),
                           engine=empty_db, header_row=0).set_index("CompoundIndex")
    assert "WARNING: could not read reference tables from the database" in capsys.readouterr().out
    pd.testing.assert_frame_equal(out, result)


def test_falls_back_to_csv_when_db_table_empty(tmp_path, template_csv, capsys):
    engine = create_engine("sqlite://")
    SALTS.iloc[0:0].to_sql("INV_LOT_SALT", engine, index=False)
    SOLVATES.to_sql("INV_LOT_SOLVATE", engine, index=False)
    salts, solvates = tmp_path / "salts.csv", tmp_path / "solvates.csv"
    SALTS.to_csv(salts, index=False)
    SOLVATES.to_csv(solvates, index=False)
    out = ess.run_pipeline(template_csv, str(salts), str(solvates), engine=engine, header_row=0)
    assert "came back empty: INV_LOT_SALT" in capsys.readouterr().out
    assert (out.salt == "hydrochloride").any()


def test_db_failure_without_csv_raises(template_csv):
    with pytest.raises(Exception):
        ess.run_pipeline(template_csv, engine=create_engine("sqlite://"), header_row=0)


def test_no_source_raises(template_csv):
    with pytest.raises(ValueError, match="engine="):
        ess.run_pipeline(template_csv, header_row=0)


# Optional: point at a throwaway MySQL/MariaDB that has INV_LOT_SALT and
# INV_LOT_SOLVATE (never production) to check the MySQL read-only path:
#   ESS_TEST_MYSQL_URL=mysql+pymysql://user:pw@localhost/testdb pytest ...
@pytest.mark.skipif(not __import__("os").environ.get("ESS_TEST_MYSQL_URL"),
                    reason="set ESS_TEST_MYSQL_URL to run against a test MySQL")
def test_mysql_read_only_transaction():
    import os
    from sqlalchemy import text
    engine = create_engine(os.environ["ESS_TEST_MYSQL_URL"])
    df = ess.read_reference_table(engine, "INV_LOT_SALT", "ILSA_NAME", "ILSA_SMARTS")
    assert not df.empty
    with engine.connect() as conn:  # connection goes back to the pool in normal mode
        assert conn.execute(text("SELECT @@autocommit")).scalar() == 0
