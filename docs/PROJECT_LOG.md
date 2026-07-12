# Perturbation PSF — Project Decision Log

**Project:** Spatial transfer functions of genetic perturbations in intact tissue.
**Hackathon deadline:** July 13, 2026.
**Dataset:** Spatial Perturb-seq, mouse hippocampus/brain — *Nat Commun* 2026, DOI `10.1038/s41467-026-69677-6`.

---

## Decision 1 — Novelty check: PASSED, with one reframing (see novelty_comparison.csv)

Ran an OpenAlex + arXiv scan across the four named comparators and the dataset paper's citers.

- **CONCERT** (bioRxiv preprint, **8 Nov 2025**, DOI `10.1101/2025.11.08.686890`) is the closest
  prior work — a niche-aware GP-VAE that learns spatial kernels to **predict** perturbation
  responses (patch/border/niche tasks), evaluated on Perturb-map **lung** data (extended clonal
  sources). It predates the dataset paper by ~2 months and MUST be cited as a co-equal comparator.
- **Perturb-FISH** (*Cell*, Feb 2025, `10.1016/j.cell.2025.02.012`): tests directly neighbouring /
  density-dependent effects — no distance-resolved kernel.
- **PerturbView** (*Nat Biotechnol* 2024): an assay, not an effect estimator.
- **Spatial CCC tools**: impose/learn distance weighting for ligand–receptor flux, not a
  perturbation impulse response.
- Only 7 works cite the dataset paper (all 2026); none estimate a distance-resolved kernel.

**Defensible novelty claim (tightened):** "To our knowledge, a null-calibrated empirical
perturbation response kernel — reporting effective radius, gain, and anisotropy with source-level
uncertainty — has not been reported for sparse single-cell spatial perturbation data, nor for this
dataset." Position PSF as **inference of a calibrated physical quantity on sparse point sources**,
contrasted with CONCERT's **prediction on extended sources**. Do NOT claim "first spatial kernel."

## Decision 2 — Data access: GO (see data_access_spike.md)

- Spatial data = 3 Stereo-seq `.cellbin.gef` files (GEO **GSE274447**), Deepcell cell-segmented,
  ~2.1 GB total. Verified Mouse 3: 65,992 cells, x/y/area/cellBorder polygons, 24,296 features
  including **all 18 sgRNA barcodes** + full transcriptome. Barcode→cell join is intrinsic.
- Dissociated scRNA-seq reference = GEO **GSE274058** (for cell-type transfer).
- GitHub `kimberle9/spatialperturbseq` ships analysis functions + a 3,757-cell DEMO object.
- **Gaps:** (1) cell-type labels NOT in .gef (cellTypeID all 0) — must re-derive (~½ day);
  (2) coordinate unit ≈ 0.5 µm (multiply coords ×0.5 for microns);
  (3) NO Xenium data deposited anywhere found — cross-platform validation demoted to stretch.
- **Power constraint:** sources are sparse PER SECTION (Lrrk2 = 6, Cfap410 = 14, msafe control = 13
  in Mouse 3). Paper's n=29 Lrrk2 is pooled across all 3 mice → MUST pool all sections. This is why
  source-level bootstrap + permutation nulls + pathway-program-first outcomes are mandatory.

## Decision 3 — MVP scope

Pool 3 sections. Perturbations: **Lrrk2** (broad neighbour effect, +control), **Cfap410** (most
neighbourhood DEGs), one weak comparator. Outcomes: pre-specified pathway programs first; per-gene
only for positive controls (Sparc, Vps35l, Dock10, Gpr37). All CPU-scale, local — no remote compute.

---

## Key files (this project)
- `novelty_comparison.csv` — 5-method comparison table
- `data_access_spike.md` — full data feasibility report
- GEO: GSE274447 (spatial), GSE274058 (scRNA-seq)
- Code/demo: github.com/kimberle9/spatialperturbseq ; Zenodo 17959756 (repo snapshot only)

---
## Session 2 — Pipeline execution (Steps 3–10)

### Decision 4 — Cell typing
Relative-enrichment marker scoring (z-score each cell-type score across 27 fine Leiden
clusters, argmax) recovered 10 mouse-brain cell types + 15 BANKSY niches. Absolute
winner-take-all failed (collapsed to 3 types). **Endothelial/Mural merged to "Vascular"**
(co-segment at Stereo-seq resolution; both express vascular markers); low-depth cluster 24
relabeled "Low-quality" (auditor-flagged, corrected + propagated downstream).

### Decision 5 — Kernel model
OLS spline bs(dist, knots=[30,60,100,160]) + treat×spline interaction, adjust
C(niche)+C(section)+C(recip_type)+log_depth, **cluster-robust SE by source_id**. K(r) =
contrast of treat + treat:spline over r∈[2,250].

### Key finding — honest, quantified null
- 170 kernels (17 guides × 10 programs), source-clustered: **2 nominal peak-sig, 0 survive**
  (min FDR 0.60).
- Four null routes agree: permutation, hierarchical bootstrap, angular-rotation,
  distant-field placebo — all consistent with no distance-resolved cell-extrinsic signal.
- Synthetic recovery: unbiased (slope 0.55), **2σ detection floor ≈0.34 near-source gain**.
  Real effects (~0.1–0.34) at/below floor at n=28–51 sources → **power limit, not absence**.
- Lrrk2 direction vs paper: 3/5 single-gene directions match (Vps35l↑, Dock10↑, Lrp1↓),
  2/5 invert (Sparc, Gpr37), none significant.

### Deliverables
unified_cells / typed_cells / neighborhoods / program_scores / kernels / all_kernels
(parquet); kernel_params / all_kernel_params / significance_table / synthetic_recovery /
leakage_qc (csv); psf_pipeline.py (runs all 18); PerturbationPSF_pipeline.ipynb;
6 result figures + figure_deck.png; SUBMISSION_SUMMARY.md.

### Xenium stretch — NOT met
247-gene Xenium panel not deposited in any locatable GEO accession. Cross-platform
validation infeasible with available data.


---

## Revision round (peer-review response) — 2026-07-09

**Trigger:** detailed peer review (verdict "do not submit unchanged"), 10 major methodological problems + biology concerns. User approved an 8-step revision reframing the work as a **calibration audit of spatial perturbation inference**.

**Decision R1 — reframe.** Headline is now the naïve-vs-source-aware ablation, not the "null kernel." Title: "How many neighbours count as one replicate? A calibration audit of spatial perturbation inference."

**Key new results:**
- Overlap audit: 50.6% of exposed recipients neighbour ≥2 perturbations; 29.5% of primary-frame recipients appear in ≥2 arms; isolated sources rare (3/guide).
- **Headline ablation:** naïve per-pair OLS 11/20 significant → source-clustered 0/20 (median SE inflation 5.3×).
- Matched pseudo-source anatomical null (replaces mislabeled angular null): anatomy-null mean is nonzero & same-signed as observed; 1/20 exceeds it.
- Global Wald test (replaces pointwise peak-sig): 5/30 nominal, 0 survive FDR (min 0.20). Source covariates now in model. T_signed & T_energy with w(r)∝2πr.
- Corrected power surface: estimator null-centred (FPR 0.04-0.06), near-unbiased on clean signal (slope ~1.0); MDAE at 80% power ≈0.5-1.0 for n≈25-50 sources; real effects (0.1-0.26) below floor.
- Reproducibility appendix: all 10 program gene sets, params, deposited-vs-target guide-name mapping (dpp5→Dpp6 etc., verified verbatim in raw .gef).

**Wording fixes applied:** regression (not penalised) spline; complementary diagnostics (not four independent nulls); null-centred/attenuated (not unbiased); pre-adjustment balance (no post-matching claim); cell-type-resolved only where stratified; PSF as analogy; "primary replication unit conditional on mouse/tissue dependence."

**Deliverables (v3):** manuscript.md/pdf v3, RESPONSE_TO_REVIEW.md, figure_deck.png v2, reproducibility_appendix.md, + audit/ablation/pseudo-source/global-test/power-surface CSVs and figures.

**Scoped-not-done:** like-for-like re-implementation of the original 15-NN Seurat test through naïve vs source-aware inference (cleanest benchmark; estimand differs from our program/near-field contrast).


---

## Round-2 revision (2026-07-10) — sealing the calibration audit

Second peer review moved the verdict to "genuinely compelling but not fully sealed" with 7 priority fixes. All completed:

1. **Bookkeeping reconciled.** 40,087 = 3-arm pool; per-contrast: Lrrk2+msafe 24,453 pairs/78 sources, Cfap410+msafe 31,189/101. Removed all "honest/true N≈90"; Fig 6c now paired per-contrast bars.
2. **Small-cluster inference primary.** Restricted wild cluster bootstrap + source-label permutation; sandwich relabelled "one-way approximation". Ablation 11→0→0→1 /20; global 5/30 sandwich → **0/30 wild** (min FDR 0.73). New Fig 6d.
3. **15-NN benchmark completed.** Original estimand, identical cells, 3 schemes: 3/5 directions reproduced, 0/5 significant; SE inflation 0.9–1.7× at 15-NN vs 5.3× at 250 µm. New Fig 8b (§3.5).
4. **Pseudo-source spec exact.** Both arms pseudo → zero-centred null; 2/20 exceed (Lrrk2 Endo/Astro, p≈0.004, FDR 0.04). Non-overlap nearest-source sensitivity concordantly flags the same two (FDR 0.28); neither survives global wild → power-limited candidates. `pseudo_source_spec.md`.
5. **Estimand/weighting corrected.** 2πr *increases* far-field weight (prior claim was wrong); empirical support ≈ area weight; adopted near-field T_near (≤40 µm) as primary. Overparam diagnostics: 55 cols full rank, cond 1768; minimal-vs-full adjustment identical (0/10); omitted Vascular composition. New Fig 14.
6. **Power surface hardened.** 1000 null reps, Wilson CIs, logistic-fit MDE (0.29–0.33 @n25 … 0.10–0.12 @n200). Null FPR 0.052–0.074, mildly anti-conservative at n=25 (4/12 CIs exclude 0.05) → reinforces wild bootstrap. Softened "resolvable at 100-200".
7. **170-screen descriptive only.** No peak significance, no stars; reach/peak not interpreted. Fig 13 rebuilt.

Also: testing families split (20 programs / 10 pos-control genes, BH separate); title → "From neighbour counts to perturbation replicates: a calibration audit of spatial perturbation inference"; figure numbering de-duplicated (1–14). Manuscript md/pdf v7; RESPONSE_TO_REVIEW_R2.md; figure_deck v3.


---

## Round-3 revision (2026-07-11) — from "close" to submission-ready

Third review: "strong and credible … very close to submission quality," no fatal flaws. Five priority fixes + additional inconsistencies, all completed:

1. **Family-mixing bug fixed.** The prior "20-program" ablation had tested 7 programs + 3 genes. Reran on the correct primary family (10 programs) and separate secondary family (5 genes): naïve **16/20 → 0/20 → 0/20 → 0/20** (all source-aware methods zero); the Cfap410→Vps35l permutation hit moved to the gene family (FDR 0.33). Figure 6d rebuilt.
2. **Source-class covariate + neuronal sensitivity.** C(source class) in the ablation and a neuronal-only run: 3/20 and 1/20 nominal, **0/20 survive FDR** under both (min q 0.32, 0.13). New Figure 6e.
3. **Better-conditioned permutation.** Section×niche×source-class strata + Freedman–Lane residual permutation, both 0/20; relabelled "primary nonparametric sensitivity analysis."
4. **Pseudo-source contradiction resolved.** Reran at 5,000 draws: two *Lrrk2* candidates now q=0.052/0.086 — 0/20 at q≤0.05, 2/20 at q≤0.10; reframed as method-dependent leads requiring validation. Figure 10 redesigned.
5. **Power per-arm + observed allocations.** Axis relabelled "sources per arm"; added observed-allocation sims (28/50, 51/50): MDE 0.14–0.19, inside the real-effect range. New Figure 11d.
6. **Weighting arithmetic corrected.** Equal-per-µm = 16%/36%; area/empirical = 3%/59%. Figure 14a shows all three.
7. **Terminology/inconsistencies.** "Interaction test" → "joint global treatment-profile test"; 15-NN "closely aligned" + "SE ratio"; cell-type-resolved → recipient-type adjustment; animal language softened; conditioning de-oversold; low-quality sensitivity rerun across all 30 contrasts (0/30); raw pointwise scan deleted, profile-max removed from Figure 13.
8. **Abstract** simplified (3,500→2,200 chars). RESPONSE_TO_REVIEW_R3.md; manuscript md/pdf v8 (19 figures); figure_deck v4.


## Round 4 revision (2026-07-11) — mouse-aware inference + estimand robustness

Fourth peer review ("scientifically defensible, suitable for submission after targeted corrections"). Eight-step plan, all complete:

1. **Mouse/chip metadata.** Adopted deposited 3-chip=3-mouse mapping (C02943C3=M1, B03018A2=M2, A03599E2=M3); noted article's 4-mouse discrepancy. 82% of both perturbations' sources in Mouse2. Added §2.1, Figure 6f, per_mouse_breakdown.csv, per_mouse_source_counts.csv. Removed "unresolved mapping" language.
2. **Leave-one-mouse-out.** 0/20 survive FDR under every mouse drop (Figure 6g, leave_one_mouse.csv). No 3-point animal bootstrap; inference conditional on 3 tissues.
3. **Equal-source estimand.** One adjusted T_s per source, Mann-Whitney + inverse-cluster-size OLS; 0/20 both (Figure 6h, equal_source_estimand.csv).
4. **Non-overlap all 20 contrasts.** Only 2 Lrrk2 candidates nominal, 0/20 FDR (Figure 10c). Overlap reframed as interference/SUTVA; exposure-basis model flagged as next step.
5. **Cell-diameter exclusion.** rmin=10/20/30µm + neuronal-only + large-mask; candidates persist/strengthen once sub-soma-diameter (<19µm) pairs removed (Figure 10d, exclusion_tests.csv). One-construct-per-gene caveat.
6. **Estimand clarity.** Total vs within-cell-type vs no-niche adjustment: 0/20 all (Figure 6i, estimand_sensitivity.csv, within_celltype_state.csv). Post-treatment-mediator flag.
7. **Power downgrade.** MDE bootstrap CIs, "idealised lower bounds ~10² sources/arm" (Figure 11d v2, mde_intervals_obs.csv). GEF resolution=500nm read from metadata. Module robustness: rank-based score weakens candidates (module_robustness.csv, module_gene_detection.csv). Title→"source-aware inference…screens". Pseudo-source "extreme tail" not "FDR support".
8. **Render.** Abstract updated; PDF v10 (25 figs, 5.02MB); RESPONSE_TO_REVIEW_R4.md; figure_deck v5.

Auditor bounce fixed: soma-diameter error (10.1µm ≈ half a 19µm diameter, not one); argument now rests on rmin≥20µm survival.

Manuscript v10: 68d436aa / PDF 85f0ff08 / response 70ef6e68 / deck 7dcf56b4 v5.
