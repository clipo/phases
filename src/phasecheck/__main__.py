"""Entry point for `python -m phasecheck`; the same as the installed `phasecheck` command."""
from .cli import main

raise SystemExit(main())
