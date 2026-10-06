"""phasecheck - do archaeological phases behave as bounded groups?

A general version of the phase tests of Lipo, DiNapoli and Madsen (the
`mls-emergence` / `phases` repository). It takes one table of assemblages
(class counts, latitude, longitude, phase) and asks, in order:

  1. Does location alone recover the phases?
  2. Do the phase boundaries separate the assemblages better than other
     lines drawn across the same map?
  3. Does each phase stand apart from the rest?
  4. Does each assemblage fit its own phase?
  5. How different are the phases, in sherds?

This is stage 1: the tests that need no copying model. It does not say what
copying across distance alone would produce on this map, so it cannot say
whether the phases differ by more than distance explains. Read every result
with that limit in mind.
"""
from .data import Dataset, read_table
from .questions import (boundary_comparison, each_phase_against_rest, own_phase_fit,
                        phase_profiles, recovery_by_location)
from .copying import calibrate, compare_with_model, power, sorting_under_model
from .report import run

__all__ = ["Dataset", "read_table", "boundary_comparison", "each_phase_against_rest",
           "own_phase_fit", "phase_profiles", "recovery_by_location", "run",
           "calibrate", "compare_with_model", "sorting_under_model", "power"]
__version__ = "0.1.0"
