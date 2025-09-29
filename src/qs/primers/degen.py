# src/qs/primers/degen.py
from __future__ import annotations
from itertools import product
from typing import Dict, Iterable, List, Sequence, Set, Tuple

IUPAC: Dict[str, str] = {
    "A":"A","C":"C","G":"G","T":"T",
    "R":"AG","Y":"CT","S":"GC","W":"AT","K":"GT","M":"AC",
    "B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT",
}
REV: Dict[frozenset, str] = {frozenset(v): k for k, v in IUPAC.items()}
DNA: Set[str] = set("ACGT")
ALL: Set[str] = set(IUPAC.keys())

def clean(seq: str) -> str:
    return (seq or "").strip().upper().replace("U", "T")

def has_degenerate(seq: str) -> bool:
    return any(ch not in DNA for ch in clean(seq) if ch.isalpha())

def degenerate_bases(seq: str) -> List[str]:
    s = clean(seq)
    return sorted({ch for ch in s if ch in ALL and ch not in DNA})

def expand_iupac(seq: str, limit: int | None = 64) -> Tuple[List[str], int]:
    s = clean(seq)
    per: List[str] = []
    total = 1
    for ch in s:
        opts = IUPAC.get(ch, "ACGT")
        per.append(opts)
        total *= len(opts)
        if limit is not None and total > limit and limit > 0:
            # too many — bail without constructing large cartesian product
            return [], total
    variants = ["".join(p) for p in product(*per)] if (limit is None or total <= (limit or 0)) else []
    return variants, total

def consensus_iupac(variants: Sequence[str]) -> str:
    if not variants:
        return ""
    vs = [clean(v) for v in variants if clean(v)]
    if not vs:
        return ""
    L = min(len(v) for v in vs)
    out: List[str] = []
    for i in range(L):
        bases = {v[i] for v in vs if len(v) > i and v[i] in DNA}
        code = REV.get(frozenset(bases), "N") if bases else "N"
        out.append(code)
    return "".join(out)

def summarize_variants(variants: Iterable[str], expand_limit: int = 64) -> Dict[str, object]:
    uniq = sorted({clean(v) for v in variants if clean(v)})
    cons = consensus_iupac(uniq) if uniq else ""
    expanded, total = expand_iupac(cons, limit=expand_limit) if cons else ([], 0)
    return {
        "consensus_iupac": cons,
        "degenerate_bases": degenerate_bases(cons) if cons else [],
        "expanded_variants": expanded,
        "expanded_count": total,
        "observed_variants_count": len(uniq),
        "observed_variants_example": uniq[:8],
    }
