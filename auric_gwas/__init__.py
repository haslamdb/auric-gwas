"""AURIC-associate — genotype new S. aureus genomes against the frozen AURIC coordinate system and run
the manuscript's lineage-blind + lineage-aware association, returning annotated ranked variants.

Read-only against the frozen panel by construction (see paths.assert_not_panel). Prototype: the `scan`
and `annotate` subcommands are functional; `genotype` (project new assemblies) and `qc` are stubs that
point at the validated cohort-genotyping fork. See docs/AURIC_ASSOCIATE_SCOPE.md.
"""
from . import paths  # noqa: F401

paths.set_thread_caps(4)

__version__ = "0.0.1"
