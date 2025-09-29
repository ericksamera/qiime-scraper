# src/qs/config/io.py
from __future__ import annotations
import json
from pathlib import Path
from typing import Any, Dict, Optional
from .schema import Params, Plan

def _read_mapping(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    try:
        import yaml  # type: ignore
        data = yaml.safe_load(text)
    except Exception:
        data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError("Params file must contain a mapping at the top level.")
    return data

def load_params_raw(path: Path) -> Dict[str, Any]:
    return _read_mapping(path)

def load_params_typed(path: Path) -> tuple[Params, Optional[Plan]]:
    data = _read_mapping(path)
    params_data = data.get("params", data)
    plan_data = data.get("plan")
    return Params(**params_data), (Plan(**plan_data) if isinstance(plan_data, dict) else None)

def write_params_and_plan(path: Path, params: Dict[str, Any], plan: Dict[str, Any]) -> None:
    payload = {"params": params, "plan": plan}
    try:
        import yaml  # type: ignore
        text = yaml.safe_dump(payload, sort_keys=False)
    except Exception:
        text = json.dumps(payload, indent=2)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")

def apply_params_defaults(args, params: Dict[str, Any], defaults: Dict[str, Any]) -> None:
    """Set argparse args from params only where current value equals our known defaults."""
    def _canon(v: Any) -> Any:
        # turn lists to comma-strings for argparse fields we expect as str
        return ",".join(v) if isinstance(v, (list, tuple)) else v
    for k, v in params.items():
        if not hasattr(args, k):
            continue
        cur = getattr(args, k)
        if k in defaults and cur == defaults[k]:
            setattr(args, k, _canon(v))
