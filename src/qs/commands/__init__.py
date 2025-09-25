# src/qs/commands/__init__.py
"""
Command package.

Submodules are imported explicitly by qs.cli to avoid circular imports.
Do NOT import submodules here.
"""
__all__ = [
    "init",
    "import_reads",
    "manifest",
    "metadata_validate",
    "metadata_tabulate",
    "trim_primers",
    "denoise_runs",
    "classify_sweep",
    "auto_run",
]
