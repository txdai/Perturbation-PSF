# Project evolution and Claude Science use

## Initial question

The project began with a cross-disciplinary hypothesis: if a genetically perturbed cell acts like a localized input to tissue, might the surrounding molecular response have a PSF-like spatial profile? This framing grew from prior electrical-engineering work on the spatial kernel of multi-resolution hash encodings for 3D reconstruction (arXiv:2602.10495). In that setting, a PSF translated a complex learned representation into resolution, anisotropy, and noise. I am not trained as a biomedical researcher; my leverage was statistical—asking what the estimand is, what unit is independent, and what result would falsify the apparent discovery. The hackathon asked whether that mathematical instinct, paired with Claude Science's domain and implementation fluency, could produce a useful biomedical measurement.

## First implementation

Claude Science located the public Spatial Perturb-seq dataset, identified the three Stereo-seq GEF files, parsed guide pseudo-genes and coordinates, harmonized genes across sections, re-derived cell classes and spatial niches, constructed source-centered recipient neighbourhoods, and fitted spline-based radial contrasts. The first analysis emphasized smooth kernels and a mostly null biological screen.

## Review-driven changes

The review was not a prewritten list of four prompts. It emerged from inspecting each result: every answer changed what I distrusted next, and Claude Science preserved the executable state as the investigation moved:

1. The apparent 16/20 result felt too precise, so I stopped the biological interpretation and challenged the replication unit. Claude implemented the naive-versus-source-aware ablation; the collapse reframed the work as a calibration audit.
2. A 0/20 result could also be fragile, so I tried to break the correction with small-cluster methods, a like-for-like 15-neighbour benchmark, corrected radial weighting, better power simulations, and strict testing-family separation. Claude implemented and reran them.
3. While inspecting the revised tables, I found that one corrected number did not align with its testing family. Claude traced and repaired the pipeline, reran it, and added source-class adjustment, a better-conditioned permutation, and observed-allocation power calculations.
4. The corrected statistics sent me back to the biological design: where did the sources actually come from? Claude resolved the mouse/chip mapping, quantified that 82% of both primary perturbations' sources lie in one mouse, and implemented equal-source, leave-one-mouse-out, non-overlap, alternative-estimand, module-score, and cell-diameter exclusion analyses.

The final finding is more modest biologically and stronger scientifically: the available data do not resolve a population-level radial response, and pair-level inference can be severely overconfident. Two Lrrk2 contrasts are retained only as method-dependent leads.

## Why the Claude use is substantive

Claude Science was not used as a prose wrapper or as the origin of the statistical critique. I directed the review; Claude maintained executable state across a multi-gigabyte dataset, generated and revised analysis modules, inspected metadata, created quantitative diagnostics, and tracked artifacts and decisions in response to my prompts. The visible evolution-from my initial PSF analogy and repeated challenges, through Claude's implementations and reruns, to a calibrated reusable audit-is the core demonstration of the collaboration.
