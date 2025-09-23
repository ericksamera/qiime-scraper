# qiime-scraper · End-to-end QIIME 2 pipeline with smart sweeps

**qiime-scraper** streamlines a QIIME 2 workflow from FASTQs → trimming/classifier selection → phylogeny → diversity analyses and visualizations. It adds a **diversity sweep** to iterate over nested metadata subsets (e.g., by `treatment_group` within `colony_source`), with **per‑subset sampling depths** computed from each subset’s `table.qzv`.

## Features

- 🔧 **Global `--dry-run`** everywhere (prints commands, no execution).
- 📣 **`--show-qiime` on by default** (live logs); use `--no-show-qiime` to capture.
- 🏆 **Automatic winner staging**: links rep‑seqs, feature table, taxonomy into `analysis/`.
- 🌲 **Phylogeny** (MAFFT → mask → FastTree → root).
- 📊 **Core diversity** + alpha/beta tests + Emperor.
- 🔁 **Two‑level diversity sweep** over metadata columns with automatic skipping of tiny subsets.
- 🎯 **Smart sampling depth per subset** via `--retain-fraction` (e.g., `0.95` for small cohorts).
- 🧾 **Provenance**: `analysis/RUN.json` and per‑subset `subset.json`.

## Quick start

> Requires QIIME 2 in your environment.

```bash
# 0) Generate manifest + import
python main.py import \
  --fastq-dir ./fastqs \
  --project-dir ./my_project

# 1) Trimming search → best pair (writes work/optimal_trimming)
python main.py trim \
  --fastq-dir ./fastqs \
  --project-dir ./my_project

# 2) Classifier sweep (uses rep_seqs from best trim)
python main.py classify \
  --project-dir ./my_project \
  --classifiers-dir /path/to/classifiers

# 3) Stage winners + build tree + core diversity + rarefaction + taxa
python main.py downstream \
  --project-dir ./my_project \
  --metadata-file ./meta.tsv \
  --beta-cols "treatment_group,colony_source" \
  --time-column date_str
````

## Diversity sweep (nested subsets)

Run all categorical metadata columns in a 2‑layer loop (outer × inner), skipping subsets with fewer than `--min-samples` rows, and choosing depth per subset from its `table.qzv` so that about `--retain-fraction` of samples are retained.

```bash
# Iterate ALL categorical columns (outer) × ALL categorical columns (inner)
python main.py diversity-sweep \
  --project-dir ./my_project \
  --metadata-file ./meta.tsv \
  --min-samples 5 \
  --retain-fraction 0.90 \
  --beta-cols "treatment_group,colony_source" \
  --time-column date_str
```

Limit to specific columns:

```bash
# Outer loop by treatment_group; inner loop by colony_source
python main.py diversity-sweep \
  --project-dir ./my_project \
  --metadata-file ./meta.tsv \
  --by treatment_group \
  --within colony_source \
  --min-samples 5 \
  --retain-fraction 0.95
```

Outputs land under:

```
analysis/
  rep-seqs.qza, table.qza, taxonomy.qza, rooted-tree.qza, ...
  RUN.json
  subsets/
    treatment_group/Control/colony_source/BC/
      table.qza, table.qzv, core-metrics-phylo/, subset.json, ...
    ...
    index.json
```

## Sampling depth logic

For each subset, we summarize the filtered table (`table.qzv`) and select a depth `d` so that roughly `retain_fraction` of samples have frequency ≥ `d`. Example: with `--retain-fraction 0.90`, about 90% of samples survive rarefaction at that depth.

## Metadata checks

Before long runs, the pipeline validates that requested columns (in `--by`, `--within`, `--beta-cols`, `--time-column`) exist in the metadata. If not, it exits with a friendly message listing available columns.

## Tips

* Use `--no-show-qiime` to keep the console quieter; errors still print captured stdout/stderr.
* Use `--if-exists overwrite|new` on sweep reruns to avoid manual cleanup.
* Keep your metadata’s `#q2:types` row accurate so the sweep auto-detects categorical columns.

## License

MIT

```