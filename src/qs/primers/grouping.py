# src/qs/primers/grouping.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Set

from qs.primers import degen
from qs.utils.text import slugify


def _set_for(ch: str) -> Set[str]:
    return set(degen.IUPAC.get((ch or "N").upper().replace("U", "T"), "ACGT"))


def iupac_edit_distance(a: str, b: str, max_edits: int = 1) -> int:
    """
    Levenshtein distance with IUPAC-aware substitutions:
      - substitution cost = 0 if the IUPAC sets intersect, else 1
      - insertion/deletion cost = 1
    Early-exits when distance will exceed max_edits.
    """
    a = (a or "").upper().replace("U", "T")
    b = (b or "").upper().replace("U", "T")
    n, m = len(a), len(b)
    if abs(n - m) > max_edits:
        return max_edits + 1

    prev = list(range(m + 1))
    cur = [0] * (m + 1)
    for i in range(1, n + 1):
        cur[0] = i
        A = _set_for(a[i - 1])
        for j in range(1, m + 1):
            B = _set_for(b[j - 1])
            sub_cost = 0 if (A & B) else 1
            cur[j] = min(
                prev[j] + 1,       # deletion
                cur[j - 1] + 1,    # insertion
                prev[j - 1] + sub_cost,   # substitution
            )
        if min(cur) > max_edits:
            return max_edits + 1
        prev, cur = cur, prev
    return prev[m]


@dataclass(frozen=True)
class Cluster:
    key: str                      # canonical "Fcons|Rcons"
    slug: str                     # filesystem-safe
    samples: List[str]            # sample-ids in this cluster
    fwd_variants: List[str]       # concrete observed F primers
    rev_variants: List[str]       # concrete observed R primers


def canonicalize_primer_pairs(
    per_sample: Dict[str, Tuple[str, str]],
    *,
    max_edits: int = 1,
) -> Tuple[Dict[str, List[str]], Dict[str, str], Dict[str, Tuple[List[str], List[str]]], List[Cluster]]:
    """
    Cluster primer pairs so that pairs within ≤ max_edits (IUPAC-aware) on BOTH F and R
    collapse together. Returns:
      groups_map: canonical_key -> [sample-ids]
      sample_to_key: sample-id -> canonical_key
      per_group_variants: canonical_key -> (list_fwd, list_rev)
      clusters: rich list with slug + variant lists
    """
    # Deduplicate identical pairs first
    pair_to_samples: Dict[Tuple[str, str], List[str]] = {}
    for sid, (f, r) in per_sample.items():
        pair_to_samples.setdefault((f or "", r or ""), []).append(sid)
    pairs = list(pair_to_samples.keys())

    # Union-Find
    parent = list(range(len(pairs)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    for i in range(len(pairs)):
        f1, r1 = pairs[i]
        for j in range(i + 1, len(pairs)):
            f2, r2 = pairs[j]
            if iupac_edit_distance(f1, f2, max_edits) <= max_edits and \
               iupac_edit_distance(r1, r2, max_edits) <= max_edits:
                union(i, j)

    # Gather clusters
    root_to_idx: Dict[int, List[int]] = {}
    for i in range(len(pairs)):
        root_to_idx.setdefault(find(i), []).append(i)

    groups_map: Dict[str, List[str]] = {}
    sample_to_key: Dict[str, str] = {}
    per_group_variants: Dict[str, Tuple[List[str], List[str]]] = {}
    clusters: List[Cluster] = []

    for idxs in root_to_idx.values():
        fset: Set[str] = set()
        rset: Set[str] = set()
        sids: List[str] = []
        for k in idxs:
            f, r = pairs[k]
            fset.add(f)
            rset.add(r)
            sids.extend(pair_to_samples[(f, r)])

        fsum = degen.summarize_variants(sorted(fset), expand_limit=64)
        rsum = degen.summarize_variants(sorted(rset), expand_limit=64)
        key = f"{fsum['consensus_iupac']}|{rsum['consensus_iupac']}"
        slug = slugify(key if key else "all")

        groups_map[key] = sorted(sids)
        for sid in sids:
            sample_to_key[sid] = key
        per_group_variants[key] = (sorted(fset), sorted(rset))
        clusters.append(Cluster(
            key=key, slug=slug, samples=sorted(sids),
            fwd_variants=sorted(fset), rev_variants=sorted(rset),
        ))

    return groups_map, sample_to_key, per_group_variants, clusters
