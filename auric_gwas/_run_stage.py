"""Launcher: run one validated CHOP genotyping stage against a request-scoped common_paths.

Mirrors scripts/ad_cohort/run_stage.py, but the cohort path module is supplied as a file path instead
of being hard-coded, so a per-request module (written by genotype.py) can be used. Registers that
module in sys.modules under the name `common_paths` BEFORE runpy so the stage's own
`import common_paths` resolves to it regardless of sys.path, then asserts every write target resolves
under the request root and refuses to run if any resolves under the frozen panel.

Usage:  python -m auric_gwas._run_stage <cohort_common_paths.py> <stage_script.py> [stage args...]
"""
from __future__ import annotations

import importlib.util
import runpy
import sys
from pathlib import Path

_PANEL_ROOTS = (
    Path("/fastpool/sausnp/db"),
    Path("/fastpool/sausnp/profiles"),
    Path("/fastpool/sausnp/profiles_masked"),
)


def _under(child: Path, parent: Path) -> bool:
    try:
        Path(child).resolve().relative_to(Path(parent).resolve())
        return True
    except ValueError:
        return False


def main() -> int:
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__ + "\n")
        return 2
    cp_path = Path(sys.argv[1]).resolve()
    stage = Path(sys.argv[2]).resolve()

    spec = importlib.util.spec_from_file_location("common_paths", cp_path)
    cp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cp)
    sys.modules["common_paths"] = cp
    cp.ensure_dirs()

    # SAFETY: every write target must be under the request root, none under the frozen panel.
    write_targets = [cp.COH, cp.INPUT, cp.PREDICT, cp.BUCKETS, cp.CALLS, cp.BUCKETS_NC,
                     cp.CALLS_NC, cp.SWITCHED_NC, cp.MATRIX, cp.QC, cp.LOGS, cp.CDS_DB_UNMASKED]
    for t in write_targets:
        for root in _PANEL_ROOTS:
            if _under(t, root):
                raise SystemExit(f"REFUSING: write target {t} is under the frozen panel {root}")
        if not _under(t, cp.COH):
            raise SystemExit(f"REFUSING: write target {t} is not under the request root {cp.COH}")

    sys.argv = [str(stage)] + sys.argv[3:]
    print(f"[auric-gwas] stage {stage.name} | common_paths COH={cp.COH}", flush=True)
    try:
        runpy.run_path(str(stage), run_name="__main__")
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
        return code
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
