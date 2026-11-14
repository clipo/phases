"""scripts/check_signs.py must flag a flipped sign and must not invent negatives.

The checker's first version passed three deliberately flipped signs because it
matched magnitudes anywhere in the outputs. These tests hold it to the cases
that matter: a value whose source line prints the other sign is flagged even
when an unrelated output prints the same magnitude with the manuscript's sign,
and hyphens that are not minus signs are not read as negatives.
"""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("check_signs", ROOT / "scripts" / "check_signs.py")
cs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cs)

# One output that is the passage's source, and one unrelated output that
# happens to print 7.3 as a negative number.
OUTPUTS = [
    ("partition.md", [[7.3, 4.6, 9.8], [6.3, 0.58, 0.65]]),
    ("elsewhere.md", [[-7.3, 112.0, 0.004]]),
]


def _verdicts(text):
    return [(tok, verdict) for verdict, _line, tok, _ctx, _why in cs.check_text(text, OUTPUTS)]


def test_correct_signs_pass():
    text = "The excess is +7.3 (interval +4.6 to +9.8), against +6.3; probabilities 0.58 and 0.65."
    assert all(v == "ok" for _t, v in _verdicts(text))


def test_a_flipped_sign_is_flagged_even_when_another_output_prints_it():
    # elsewhere.md prints -7.3, which is what cleared this in the first version.
    text = "The excess is -7.3 (interval +4.6 to +9.8), against +6.3; probabilities 0.58 and 0.65."
    assert ("-7.3", "FLIPPED?") in _verdicts(text)


def test_the_typographic_minus_is_a_minus():
    text = "The excess is \N{MINUS SIGN}7.3 (interval +4.6 to +9.8)."
    assert [v for t, v in _verdicts(text) if t.endswith("7.3")] == ["FLIPPED?"]


def test_hyphens_that_are_not_signs_are_not_negatives():
    text = ("Ranges such as 10-15 km, k-means, the form 1e-5, a citation [-@lipo_2001, 60], "
            "and a year span 1996-97 carry no signed value.")
    assert _verdicts(text) == []


def test_half_way_rounding_matches_either_way():
    outputs = [("corr.md", [[-0.140, -0.435, 0.204]])]
    text = "posterior median -0.14, interval -0.44 to +0.20"
    assert [v for _t, v in [(t, v) for v, _l, t, _c, _w in cs.check_text(text, outputs)]] == ["ok"] * 3
