# src/qs/config/schema.py
from __future__ import annotations
from typing import Dict, List, Optional
from pydantic import BaseModel, Field, field_validator

class PrimerSideSummary(BaseModel):
    consensus_iupac: str = ""
    degenerate_bases: List[str] = Field(default_factory=list)
    expanded_variants: List[str] = Field(default_factory=list)
    expanded_count: int = 0
    observed_variants_count: int = 0
    observed_variants_example: List[str] = Field(default_factory=list)

class PrimerGroupSummary(BaseModel):
    group_key: str = ""
    group_slug: str = ""
    n_samples: int = 0
    forward: PrimerSideSummary = Field(default_factory=PrimerSideSummary)
    reverse: PrimerSideSummary = Field(default_factory=PrimerSideSummary)

class Params(BaseModel):
    # selection / resume
    groups: Optional[List[str]] = None
    resume: bool = True

    # denoise
    split_by_primer_group: bool = True
    cores: int = 0
    discard_untrimmed: bool = False
    no_indels: bool = False

    # auto truncation
    auto_trunc: bool = False
    trunc_len_f: int = 0
    trunc_len_r: int = 0
    trunc_lower_frac: float = 0.80
    trunc_min_len: int = 50
    trunc_step: int = 10
    trunc_refine_step: int = 5
    trunc_quick_learn: int = 250000

    # dada2 misc
    trim_left_f: int = 0
    trim_left_r: int = 0
    max_ee_f: int = 2
    max_ee_r: int = 2
    trunc_q: int = 2
    min_overlap: int = 12
    pooling_method: str = "independent"
    chimera_method: str = "consensus"

    # filtering
    do_filter: bool = True
    rare_freq_frac: float = 0.001
    rare_min_samples: int = 1
    contam_exclude: str = "mitochondria,chloroplast"
    filter_unclassified_phylum: bool = False
    min_sample_depth: int = 0

    # core metrics
    retain_fraction: float = 0.90
    min_depth: int = 1000
    alpha_tsv: bool = False
    taxa_barplot: bool = True

    # stats
    run_stats: bool = False
    beta_group_cols: Optional[List[str]] = None
    beta_method: str = "permanova"
    pairwise: bool = False
    permutations: int = 999
    adonis_formula: Optional[str] = None

    # display
    show_qiime: bool = True

    # classification
    classifiers_map: Optional[Dict[str, str | None]] = None

    # detected primer pairs (lightweight)
    primer_pairs: Dict[str, Dict[str, str]] = Field(default_factory=dict)

    # NEW: grouping tolerance for indel/mismatch-tolerant grouping
    group_edit_max: int = 1

    @field_validator("beta_method")
    @classmethod
    def _check_beta(cls, v: str) -> str:
        if v not in ("permanova", "anosim", "permdisp"):
            raise ValueError("beta_method must be one of: permanova, anosim, permdisp")
        return v

class Plan(BaseModel):
    generated: str
    detected_groups_n: int = 0
    primer_pairs: List[PrimerGroupSummary] = Field(default_factory=list)
    detected_primers_per_sample: Dict[str, List[str] | str | tuple] = Field(default_factory=dict)
    suggested_command: str = ""
    notes: str = ""
