# Response to Review — Round 4

We thank the reviewer for a careful, constructive read. The manuscript was already framed as a calibration/design audit; this round strengthens that framing and closes the seven highest-priority items and every additional point. Section references are to the revised manuscript; new artifacts are named inline.

## The seven highest-priority items

**1. Mouse/chip metadata is no longer "unresolved."**
The reviewer is right: the deposited GEO samples map the three chips one-to-one to three mice (C02943C3 = Mouse 1, B03018A2 = Mouse 2, A03599E2 = Mouse 3). We adopt this deposited mapping, note the article-text discrepancy (four mice), and treat mouse/chip as the top-level blocking variable (new §2.1). The design consequence is stark and now stated everywhere it matters: **82% of both perturbations' sources lie in Mouse 2** (23/28 *Lrrk2*, 42/51 *Cfap410*; Mouse 1 has 1–2 per guide), so effective animal replication is ≈1. We add per-mouse source counts and per-mouse effect estimates (Figure 6f), **leave-one-mouse-out** stability (0/20 survive FDR under every drop; Figure 6g), and we explicitly do **not** report a three-point animal bootstrap. Inference is stated as conditional on the three realized tissues, not a population sample. `per_mouse_breakdown.csv`, `per_mouse_source_counts.csv`, `leave_one_mouse.csv`.

**2. Equal-source-weighted primary estimand.**
We add an explicit equal-source estimand (§2.7, §3.4): each source is reduced to one adjusted near-field summary T_s (mean residual after pre-specified baseline adjustment), and arms are compared with one weight per source (Mann–Whitney/Welch), plus an inverse-cluster-size-weighted OLS. Both agree with the pair-weighted analysis — 0/20 surviving FDR (min p=0.10) — so the conclusion does not depend on recipient-count weighting (Figure 6h, `equal_source_estimand.csv`).

**3. Non-overlap analysis for all 20 contrasts; overlap framed as interference.**
The nearest-source non-overlap design (one recipient → one source, no cell shared between arms) is now run over **all 20 primary contrasts**, not just the candidates: only the two *Lrrk2* contrasts are nominal, 0/20 survive FDR (Figure 10c, `nearest_source_all20.csv`). We reframe recipient overlap as an interference/SUTVA problem, note that our delivered pairwise model does not fit the originally-proposed additive sum-over-sources term despite extensive overlap, and flag a per-recipient exposure-basis model plus two-way (source × recipient) clustering as the priority next step (Discussion §4).

**4. Cell-diameter exclusion for the two *Lrrk2* candidates.**
We add technical-exclusion tests (§3.6, Figure 10d, `exclusion_tests.csv`). At the ~19 µm median soma diameter, the nearest pair (10.1 µm) is about half a diameter, and 219 near-field pairs lie within 20 µm — these touching/overlapping pairs are exactly the leakage risk. Removing them (rmin=20 µm) does not eliminate either candidate (astrocyte-reactivity strengthens to p=0.001; endolysosomal holds at p=0.038); neuronal-only and large/fragmented-mask exclusions agree. Because there is one barcoded construct per gene, we relabel any surviving signal a "*Lrrk2*-targeting-construct-associated contrast." We additionally show the candidates are partly scoring-method-dependent (weaken under a rank-based module score; `module_robustness.csv`).

**5. 15-NN benchmark downgraded.**
Relabelled throughout as a "closely aligned independent reanalysis using our own guide calls/QC/expression," explicitly "not a bit-exact reproduction of the published pipeline" (§3.5).

**6. Power surface downgraded to illustrative lower bounds.**
§3.7 now states the DGP is deliberately optimistic (lists the omitted realism factors and notes it uses a source-clustered Wald fit, not our primary wild bootstrap/permutation), reports MDEs with bootstrap 95% CIs, and frames them as lower bounds "on the order of 10² sources per arm" (Figure 11d, `mde_intervals_obs.csv`).

**7. Technical exclusion tests.** Delivered as item 4 above (inner-shell, neuronal-only, large-mask exclusions), plus the depth-vs-distance leakage QC already in §2.9.

## Additional points

- **Coordinate scale read from metadata, not inferred.** The GEF root `resolution` attribute reads 500 (nm) in all three files; §2.1 now cites this directly (the ~19 µm soma cross-check is retained as corroboration).
- **Neuronal sources as primary / source class.** Source class enters as an explicit `C(source class)` covariate and a neuronal-only sensitivity; verdict unchanged (§2.2, §3.4, Figure 6e).
- **Module-score robustness.** Added a rank-based recomputation (per-gene quantile-normalised across cells, averaged over module members — a variance-insensitive alternative to `score_genes`, not the exact within-cell UCell statistic) and a per-gene detection table (`module_robustness.csv`, `module_gene_detection.csv`); the candidate signals are shown to be scoring-dependent, reinforcing "leads, not findings."
- **Estimand clarity.** New §3.4 subsection distinguishes the total local effect from the composition-adjusted within-cell-state contrast, runs all three adjustment sets (0/20 each; Figure 6i), and flags composition/niche as possible post-treatment mediators (§2.2).
- **Pseudo-source terminology.** "FDR support" replaced by "extreme tail of the matched null … not supported by the primary global test."
- **Title.** Changed to "From neighbour counts to source-aware inference: a calibration audit of spatial perturbation screens"; "Perturbation PSF" retained as the software name only.
- **"pre-specified"** used throughout (never "pre-registered").

## Claim-by-claim table (reviewer's calibration)
We align with the reviewer's assessment. Strongly supported: pair-level inference understates uncertainty. Supported conditional on estimand/model: no radial profile is resolved. Not supported / not claimed: identified non-cell-autonomous effects, a validated causal PSF, or refutation of the original DEGs — we now state these explicitly as out of scope. Exploratory only: the two *Lrrk2* candidates. Directionally plausible but numerically optimistic: the source-count targets, now presented as lower bounds.
