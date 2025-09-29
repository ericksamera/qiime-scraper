# src/qs/config/load.py
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from qs.config.schema import Params


def _read_yaml_or_json(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    # Prefer YAML if available, fall back to JSON
    try:
        import yaml  # type: ignore
        data = yaml.safe_load(text)
    except Exception:
        data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError("Params file must contain a mapping/object at the top level.")
    return data


def load_params_file(path: Optional[Path]) -> Optional[Params]:
    if not path:
        return None
    data = _read_yaml_or_json(path)
    return Params.from_mapping(data)


def apply_params_defaults(args, params: Optional[Params]) -> None:
    """
    Merge params into argparse args **only where the user left CLI at defaults**.
    CLI must always win; params just save typing. We compare against _DEFAULTS stored in auto_run.
    """
    if not params:
        return
    from qs.commands.auto_run import _DEFAULTS  # avoids a hard cycle at import time

    p = params.to_dict()
    # map dataclass into args fields that exist
    for name, value in p.items():
        if not hasattr(args, name):
            continue
        if name in ("groups", "beta_group_cols") and isinstance(value, list):
            # CLI expects comma string (existing codepath); keep that style here
            value = ",".join(value)
        current = getattr(args, name)
        default = _DEFAULTS.get(name, None)
        if current == default:
            setattr(args, name, value)
