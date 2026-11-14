"""scripts/build_manifest.py must credit a script with what it writes, not with what it reads (rule 15).

A frame assigned from pd.read_csv(OUTD / "a.csv") and then written with
frame.to_csv(OUTD / "b.csv") was credited with a.csv, because the receiver of
a method call was tried as the path before the arguments. Reads are
deliberately not counted; crediting one is the error that makes a manifest lie.
"""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_manifest", ROOT / "scripts" / "build_manifest.py")
bm = importlib.util.module_from_spec(spec); spec.loader.exec_module(bm)


def test_dataframe_writer_is_credited_with_its_argument_not_its_source(tmp_path):
    script = tmp_path / "s.py"
    script.write_text(
        'from pathlib import Path\nimport pandas as pd\nROOT = Path(".")\nOUTD = ROOT / "output" / "x"\n'
        'sc = pd.read_csv(OUTD / "a.csv")\nsc.to_csv(OUTD / "b.csv", index=False)\n'
        'OUT = ROOT / "output" / "findings" / "c.md"\nOUT.write_text("text")\n')
    outs, _ = bm.python_outputs(script)
    assert outs == {"output/x/b.csv", "output/findings/c.md"}


def test_r_write_json_is_counted(tmp_path):
    # Six probe R scripts write their only artifact with jsonlite::write_json(..., here(...)); the
    # R pattern matched writeLines, saveRDS, ggsave and write.csv only, so they were absent from the manifest.
    script = tmp_path / "s.R"
    script.write_text('library(here)\nx <- list(a = 1)\njsonlite::write_json(x, here("output", "probe", "s.json"), auto_unbox = TRUE)\n')
    assert bm.r_outputs(script) == {"output/probe/s.json"}


def test_declared_template_outputs_are_listed(tmp_path):
    # Rule 15 gap found by the whole-branch review (2026-10-19): a script whose output names come from command-line
    # arguments writes nothing the builder can resolve. It declares its templates in MANIFEST_OUTPUTS; the builder
    # lists them as written, braces and all, so the manifest names every artifact family and its producer.
    script = tmp_path / "s.py"
    script.write_text('from pathlib import Path\nROOT = Path(".")\nMANIFEST_OUTPUTS = ["output/x/table_{generator}.csv", "output/findings/t_{point}.md"]\n'
                      'def main(a):\n    (ROOT / "output" / "x" / f"table_{a}.csv").write_text("")\n')
    outs, _ = bm.python_outputs(script)
    assert outs == {"output/x/table_{generator}.csv", "output/findings/t_{point}.md"}
    bad = tmp_path / "t.py"; bad.write_text('MANIFEST_OUTPUTS = ["notes/x.md"]\n')
    import pytest
    with pytest.raises(SystemExit):
        bm.python_outputs(bad)                                   # a declared path outside the tracked directories is an error, not silence


def test_r_file_path_with_a_here_variable_is_counted(tmp_path):
    # validate_a1.R and a1_prior_check.R write their JSON summaries as write_json(x, file.path(OUTD, "name.json"))
    # with OUTD <- here("output", "probe_latent_groups"); the pattern resolved only here(...) chains in the call.
    script = tmp_path / "s.R"
    script.write_text('library(here)\nOUTD <- here("output", "probe")\nx <- list(a = 1)\n'
                      'write_json(x, file.path(OUTD, "summary.json"), auto_unbox = TRUE)\nwriteLines("a", here("output", "findings", "f.md"))\n')
    assert bm.r_outputs(script) == {"output/probe/summary.json", "output/findings/f.md"}
