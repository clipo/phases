"""phasecheck: do archaeological phases behave as bounded groups?

A general version of the phase tests of Lipo, DiNapoli and Madsen (the
`mls-emergence` / `phases` repository). It takes one table of assemblages
(class counts, latitude, longitude, phase) and asks, in order:

  1. Does location alone recover the phases?
  2. Do the phase boundaries separate the assemblages better than other
     lines drawn across the same map?
  3. Does each phase stand apart from the rest?
  4. Does each assemblage fit its own phase?
  5. How different are the phases, in sherds?

Questions 1 to 5 need no model of copying (`questions.py`). They compare the
phase lines with other lines on the same map, so on their own they cannot say
whether the phases differ by more than copying across distance explains.
With the copying model (`copying.py`, `--model` on the command line) four
more questions answer that:

  6. Can copying with no groups reproduce this record's diversity?
  7. Does the record differ between phases by more than copying across
     distance produces?
  8. Would assemblages made with no groups sort into these phases?
  9. Could this comparison detect a restriction at the phase boundaries?

Module map: `data` reads and checks the table; `geo` turns coordinates into
distances; `measures` holds the measures and their sign conventions;
`divisions` builds the alternative divisions of the map; `questions` answers
1 to 5; `model` and `copying` answer 6 to 9; `report` runs everything and
writes the report; `cli` is the command line.

Conventions used throughout: distances are in km; cultural F_ST is computed
on Gini-Simpson diversity; the boundary excess is within-phase minus
between-phase similarity at matched distance, so a positive value means the
phase lines separate assemblages more than distance alone predicts. Every
probability is a share of posterior draws or of simulated runs, never a
significance level.
"""
from .copying import calibrate, compare_with_model, power, sorting_under_model
from .data import Dataset, read_table
from .questions import (
    boundary_comparison,
    each_phase_against_rest,
    own_phase_fit,
    phase_profiles,
    recovery_by_location,
)
from .report import run

__all__ = [
    "Dataset",
    "boundary_comparison",
    "calibrate",
    "compare_with_model",
    "each_phase_against_rest",
    "own_phase_fit",
    "phase_profiles",
    "power",
    "read_table",
    "recovery_by_location",
    "run",
    "sorting_under_model",
]
__version__ = "0.1.0"
