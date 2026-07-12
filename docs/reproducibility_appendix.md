# Reproducibility appendix — Perturbation PSF

## A. Pathway program gene sets

Programs were defined from canonical marker/pathway literature for mouse
neurodegeneration and glial biology **before** examining distance-dependence in
these data (they were fixed at the outcome-definition step, prior to any kernel
fit). Module scores use scanpy `score_genes` with `ctrl_size=50` random
background genes, computed on log-normalised (10⁴ counts/cell, log1p) whole-transcriptome
expression. Genes absent from the panel were dropped; counts below are genes **present**.

| Program | n genes | Genes |
|---|---|---|
| Synaptic_signaling | 15 | Snap25, Syt1, Syn1, Dlg4, Nrgn, Camk2a, Grin1, Gria1, Nlgn1, Nrxn1, Syp, Bsn, Homer1, Shank3, Slc17a7 |
| Calcium_homeostasis | 11 | Calb1, Calb2, Calm1, Camk2a, Atp2b1, Ryr2, Cacna1a, Cacna1b, Itpr1, Pvalb, Ncald |
| Neuroinflammation | 11 | Il1b, Tnf, Nfkb1, Ccl2, Cxcl10, Il6, Tlr4, Nlrp3, Ccl3, Ccl4, Il18 |
| Interferon_response | 11 | Irf7, Isg15, Ifit1, Ifit3, Stat1, Oasl2, Usp18, Rsad2, Mx1, Ifitm3, Bst2 |
| Microglial_activation | 11 | Cx3cr1, Tmem119, P2ry12, Csf1r, Aif1, Itgam, Cd68, Trem2, Tyrobp, C1qa, Ctss |
| Astrocyte_reactivity | 10 | Gfap, Vim, Serpina3n, Lcn2, C3, S100b, Aqp4, Cd44, Nes, Slc1a3 |
| ECM | 12 | Sparc, Col1a1, Fn1, Bgn, Dcn, Lama2, Postn, Mmp2, Timp1, Vcan, Hapln1, Ncan |
| Endolysosomal | 11 | Vps35, Vps35l, Lamp1, Lamp2, Ctsd, Ctsb, Rab7, Sqstm1, Gpr37, Atp6v0c, Cd63 |
| Mitochondrial_stress | 9 | Ndufa1, Ndufaf2, Sod1, Sod2, Hspa9, Ppargc1a, Cox4i1, Nnt, Bnip3 |
| LRP1_signaling | 8 | Lrp1, App, Apoe, A2m, Lrpap1, Ctsd, Calr, Hspg2 |

### Positive-control genes (from the original study's Lrrk2 neighbour DEGs)
- **Lrrk2_neighbor_up** (expected ↑): Vps35l, Dock10, Gpr37
- **Lrrk2_neighbor_down** (expected ↓): Sparc, Lrp1

Scored as **individual genes** (`expr_<Gene>`), separately from the programs.

### Testing-family clarification
`Vps35l`, `Dock10`, `Gpr37`, `Sparc`, `Lrp1` are **individual positive-control
genes**, *not* among the ten programs. When the text reports a per-gene contrast
(e.g. "Cfap410 → Vps35l"), the statistic belongs to the positive-control gene
family, and multiple-testing correction is applied within the panel actually
tested. `Vps35l` and `Gpr37` additionally appear inside the Endolysosomal
program, and `Sparc` inside ECM; the program score is a multi-gene module and is
not interchangeable with the single gene (this is why the Lrrk2 *Lrp1* gene can
trend down while the multi-gene *LRP1_signaling* program trends up).

## B. Guide calling and the deposited-vs-target name mapping

Guide-positive calling used a strict threshold of **≥2 guide UMIs** per cell.
The 18 guide barcodes are stored in the deposited `.gef` files under short
labels. These labels are the authors' deposited names, reproduced verbatim by
our pipeline (verified by reading `cellBin/gene` directly from the raw HDF5);
several are looser than the canonical gene symbol they most plausibly target:

| Deposited guide label | Most plausible target gene |
|---|---|
| sgrna_dpp5 | Dpp6 (or Dpp5/Dpp10 family) |
| sgrna_ndufaf | Ndufaf2 |
| sgrna_rbfox | Rbfox3 |
| sgrna_clu | Clu |
| sgrna_lrrk2 | Lrrk2 |
| sgrna_Cfap410 | Cfap410 |
| sgrna_msafe | safe-harbour control (non-targeting) |

The remaining labels (C9orf72, fasn, flcn, gfap, oligo2, rraga, sh3gl2, srf,
stk39, tbk1, trem2) match their gene symbols directly. **No truncation or
remapping occurred in our code** — the labels are as deposited. Analyses key on
the guide label, so the target-name ambiguity does not affect any count or
contrast; it only affects how a label is described in prose.

## C. Key parameters

| Parameter | Value |
|---|---|
| Coordinate scale | 0.5 µm / raw Stereo-seq unit |
| QC filter | ≥100 counts, ≥50 genes → 222,140 cells |
| Guide threshold | ≥2 UMI |
| Neighbourhood radius R | 250 µm |
| Near-field contrast | r ≤ 40 µm |
| Spline basis | cubic B-spline, interior knots 30/60/100/160 µm (**regression** spline, unpenalised OLS) |
| Kernel model | OLS, cluster-robust SE by source cell |
| Adjustments | recipient niche, section, cell type, log depth; + 14 source covariates (depth, local composition, density, boundary) |
| Global test | Wald test of joint {treat, treat×spline}=0 |
| Test statistics | T_signed = ∫w(r)K(r)dr and T_energy = ∫w(r)K(r)²dr, w(r) ∝ 2πr |
| Multiple testing | Benjamini–Hochberg across the perturbation×program family |
| Module score | scanpy score_genes, ctrl_size=50, on log-norm expression |

## D. Software
Python 3.11; scanpy 1.11.5 (typing/scoring; run with `NUMBA_CACHE_DIR` set
because the conda site-packages is read-only), statsmodels 0.14.6 (kernels,
Wald tests, cluster-robust SE), scipy cKDTree (neighbourhoods). The full
pipeline is `psf_pipeline.py`; the notebook is `PerturbationPSF_pipeline.ipynb`.
