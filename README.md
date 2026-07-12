# Perturbation PSF

**A source-aware calibration audit for sparse spatial perturbation screens.**

Spatial Perturb-seq measures many recipient cells around a smaller number of genetically perturbed sources. This repository tests a critical inference failure mode: the same 20 biological contrasts produce 16 nominal hits when 24,453-31,189 source-recipient pairs are treated as independent, but 0/20 when uncertainty is computed at the 78 or 101 perturbation and control sources.

The project was built with Claude Science by an electrical-engineering researcher entering an unfamiliar biomedical field. Human-directed statistical review repeatedly changed the analysis; Claude maintained the executable research state, implemented the reruns, and helped turn the final calibration result into reusable software.

## Read the paper

The compiled article is at [`paper/PerturbationPSF_submission.pdf`](paper/PerturbationPSF_submission.pdf). LaTeX source is intentionally excluded from this submission repository.

## Repository contents

- `src/perturbation_psf/` - reusable source-aware inference library.
- `scripts/run_audit.py` - command-line audit for long-form source-recipient tables.
- `scripts/make_synthetic_data.py` - deterministic smoke-test dataset generator.
- `scripts/00_download_data.sh` - download instructions for the public GEO data.
- `tests/` - unit tests for the central inference code.
- `notebooks/` - analysis notebook.
- `analysis/` - archived full Claude Science analysis pipeline.
- `data/results/` - review-validated result tables used in the paper.
- `data/synthetic/` - small deterministic example data.
- `docs/` - reproducibility notes, project evolution, decision log, and review response.
- `SUBMISSION_SUMMARY.md` - 200-word project description.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
python scripts/make_synthetic_data.py
python scripts/run_audit.py \
  --input data/synthetic/source_recipient_frame.csv \
  --outdir data/generated \
  --outcomes program_a program_b \
  --covariates log_depth \
  --wild-reps 199
pytest -q
```

The command writes:

- `data/generated/replication_unit_audit.csv`
- `data/generated/wild_cluster_bootstrap.csv`
- `data/generated/equal_source.csv`
- `data/generated/summary.json`

## Full public dataset

The approximately 2 GB raw archive is not vendored. Download it with:

```bash
bash scripts/00_download_data.sh
```

- GEO accession: [GSE274447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274447)
- Associated study: [10.1038/s41467-026-69677-6](https://doi.org/10.1038/s41467-026-69677-6)
- Prior PSF work: [arXiv:2602.10495](https://arxiv.org/abs/2602.10495)

## License

MIT. See [`LICENSE`](LICENSE).

