# thyroidBRS

<!-- badges: start -->
[![R-CMD-check](https://github.com/camilasauria/thyroidBRS/actions/workflows/check.yaml/badge.svg)](https://github.com/camilasauria/thyroidBRS/actions/workflows/check.yaml)
[![Bioconductor build](https://bioconductor.org/shields/build/devel/bioc/thyroidBRS.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/thyroidBRS/)
[![Bioconductor years](https://bioconductor.org/shields/years-in-bioc/thyroidBRS.svg)](https://bioconductor.org/packages/thyroidBRS)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23271086.svg)](https://doi.org/10.5281/zenodo.23271086)
<!-- badges: end -->

An R package that computes the **BRAF-RAS Score (BRS)** for papillary
thyroid carcinoma (PTC) — a 71-gene expression signature that quantifies how
much a tumor resembles a BRAF<sup>V600E</sup>-mutant or a RAS-mutant transcriptional
profile — following the description in:

> Agrawal N, Akbani R, Aksoy BA, et al. Integrated Genomic Characterization
> of Papillary Thyroid Carcinoma. *Cell*. 2014;159(3):676-690.
> [doi:10.1016/j.cell.2014.09.050](https://doi.org/10.1016/j.cell.2014.09.050)

## Why this package exists

The original TCGA-THCA study published a BRS for 391 of its samples, because
computing it required **both** exome sequencing (to define the BRAF-mutant and
RAS-mutant reference groups) **and** RNA sequencing (to compute the score).
But once the two reference centroids are fixed, scoring an *additional* tumor
only needs its gene expression — no mutation calling required. At the time of
writing no public package implements this classifier, and the gene list is
published only as a heatmap figure (Figure S7A), not as a table.

`thyroidBRS`:

- ships the 71-gene signature, read off Figure S7A and reconciled against
  current HGNC symbols (`brs_genes`);
- fits the two reference centroids (`brs_fit()`) from any cohort with known
  driver-mutation labels;
- scores *any* additional sample against those centroids (`predict()`),
  including samples with no mutation data at all.

## The score

Following Extended Experimental Procedures 14.1, the score of a tumor `t` is
the difference between its normalized Euclidean distance to the two centroids:

```
BRS(t) = ||v(t) - c(B)||₂ - ||v(t) - c(R)||₂
```

Negative is BRAF<sup>V600E</sup>-like, positive is RAS-like (`brs_class` spells
the two classes `"Braf-like"` and `"Ras-like"`). `predict()` returns both
this value (`brs_score`) and the version rescaled to [-1, +1] across the
scored samples (`brs_scaled`), which is the axis the published figures use.
The rescaling is a property of the *set* being scored, not of a sample on its
own — use `brs_score` when you need a per-sample value that does not move when
the cohort changes.

## Installation

`thyroidBRS` was accepted into Bioconductor in October 2026 and is in the
development branch (devel), which requires the current R:

```r
# install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::install("thyroidBRS")
```

From the next Bioconductor release on, the release version installs with
`BiocManager::install("thyroidBRS")` alone. To track the GitHub repository
directly:

```r
# install.packages("remotes")
remotes::install_github("camilasauria/thyroidBRS")
```

## Quick start

```r
library(thyroidBRS)

# `expr` and `labels` stand for your own data; the vignette builds a small
# simulated pair that runs as is.
# expr:   log2(TPM + 1) expression matrix, gene symbols x samples
# labels: named character vector of DRIVER MUTATION status, for the subset of
#         samples that have it — "BRAF_V600E" / "RAS"
fit <- brs_fit(expr, labels)
predict(fit, expr)   # scores every sample, including unlabeled ones
```

**Fit on mutation status, not on published BRS classes.** The two reference
sets in the paper are the BRAF-V600E-mutant and RAS-mutant tumors. The
TCGA-THCA column `BRAF_RAF_class` (also exposed by
`TCGAbiolinks::TCGAquery_subtype()`) is *not* that — it is exactly
`sign(BRAF_RAF_score)`, the output of the classifier. Fitting on it is
circular. Use `BRAFV600E_RAS`, restricted to `BRAF_V600E` and `RAS`, and keep
the published score and class for validation. See `vignette("thyroidBRS")`.

## What this reproduces, and what it does not

This is a re-implementation from the published description. The original code
was never released, and neither was the gene list except as a heatmap. So the
honest claim is not that this returns the published numbers, but that it
reproduces the published *classification*:

- against driver-mutation status, which is what the reference groups are
  defined by, 10-fold cross-validation gives **99.7%**;
- against the published class labels, **96.4%** on samples held out of the
  fit;
- against the published continuous score, **Spearman 0.941** — close, but not
  identical.

That last number is the point. Expression pipelines have moved (RSEM on hg19
then, STAR on GENCODE v36 now), and that difference cannot be reconstructed.
Use this to classify samples and to rank them along the BRAF–RAS axis. Do not
expect it to return the exact values printed in the 2014 supplementary table.

## How well does it work

Refitting on TCGA-THCA — centroids from the 234 BRAF-V600E-mutant and 52
RAS-mutant tumors, then scoring all 505 primary tumors available today:

| Comparison | Result |
| --- | --- |
| Spearman vs. published `BRAF_RAF_score` (n = 391) | 0.941 |
| Pearson, rescaled to [-1, 1] (n = 391) | 0.986 [^1] |
| Class concordance, samples **held out** of the fit (n = 111) | 96.4% |
| Class concordance, all 391 published samples | 99.0% |
| 10-fold CV against mutation status (n = 286) | 99.7% |
| Patients scored | 505, of which 114 have no published BRS |

**Read these carefully, because none of them is a clean generalization
estimate.**

The 10-fold CV re-fits the centroids in each fold, but the *signature* does
not move: those 71 genes were selected on these same TCGA BRAF/RAS tumors in
2014. Gene selection sits outside the cross-validation loop, so 99.7% is
optimistic in exactly the way this package criticizes resubstitution for
being. Treat it as an upper bound.

The 96.4% is the cleanest number here. Those 111 tumors are neither
BRAF-V600E- nor RAS-mutant, so they entered neither the signature derivation
nor the centroid fit — 58 carry a fusion driver, 10 another driver mutation,
and 43 no identified driver at all. All 111 were exome sequenced. (Of the
114 scored tumors without a published BRS, 94 were never exome sequenced, 9
were but received no score, and 11 are not in the clinical table at all;
none of them is in this set.) It still measures agreement with the original
classifier rather than with a biological truth.

Concordance on the samples that *defined* the centroids is resubstitution and
carries a large optimistic bias: at the dimensions of the published
comparison (70 genes, 391 samples split 272/119), pure noise with no signal
at all still resubstitutes at about 68%. `validate_brs()` reports how much of
any comparison is resubstitution, and warns when all of it is.

A genuinely independent estimate needs a cohort that played no part in
deriving the signature. See
[`data-raw/microarray_notes.md`](https://github.com/camilasauria/thyroidBRS/blob/main/data-raw/microarray_notes.md).

[^1]: `brs_scaled` is rescaled across whatever set is being scored, and the
    published values were rescaled across the 391 tumors that had one. The
    rescaling is linear within each sign but not across zero, so Pearson is
    not invariant to it and the two are only comparable when both extremes
    fall inside the shared subset. Here they do — the most negative and most
    positive of all 505 tumors are both among the 391, so scoring 505 or 391
    gives 0.986 either way — but that is a property of this cohort, not a
    guarantee. Rescale over the shared subset before comparing elsewhere.

Reproduce with
[`data-raw/validate_tcga.R`](https://github.com/camilasauria/thyroidBRS/blob/main/data-raw/validate_tcga.R).

## Gene list

The 71 genes were transcribed from the published heatmap, which is the only
form in which they were released, and then verified row by row against it.
Four checks back the transcription:

- Both row-clusters of Figure S7A are printed alphabetically — 13 genes and
  58 genes — and neither has a gene out of order, so an omission or a misread
  would show up. `13 + 58 = 71` matches the figure legend.
- Refitting the centroids on TCGA-THCA recovers the block structure exactly:
  all 13 block-1 genes are higher in the RAS group and all 57 usable block-2
  genes higher in the BRAF group, 70 of 70 with no exceptions.
- Symbols were resolved against `org.Hs.eg.db`. Four are stale in the 2014
  figure (`ARNTL` → `BMAL1`, `PVRL4` → `NECTIN4`, `FAM176A` → `EVA1A`,
  `TM7SF4` → `DCSTAMP`) and `FLJ23867` maps to nothing in current annotation,
  so the classifier uses 70 genes.
- Re-running the derivation itself
  ([`data-raw/derive_signature.R`](https://github.com/camilasauria/thyroidBRS/blob/main/data-raw/derive_signature.R))
  puts the
  published genes at the top of the ranking: over 200 limma-voom iterations
  only ~260 genes ever reach a top 100, and 68 of the 70 usable published
  genes are among them, at a median selection frequency of 93-97%. On a
  2014-like gene universe the paper's strict criterion keeps 27 genes, 26 of
  them published.

`brs_fit()` resolves each gene against whichever spelling your matrix actually
uses, so GENCODE v36 matrices (`ARNTL`, `PVRL4`) and current ones (`BMAL1`,
`NECTIN4`) both give the full 70.

What the re-derivation does *not* reproduce is the count of exactly 71, which
depends on choices the paper leaves unspecified. See
[`data-raw/signature_rederivation.md`](https://github.com/camilasauria/thyroidBRS/blob/main/data-raw/signature_rederivation.md)
for the full comparison, and
[`data-raw/brs_genes.md`](https://github.com/camilasauria/thyroidBRS/blob/main/data-raw/brs_genes.md)
for the provenance of the list.

## Provenance

Substantial parts of this package were written by Claude models (Anthropic),
working from the published description of the method and under the direction
of the authors, who reviewed the work and are responsible for its correctness
and maintenance: the first version was drafted with **Claude Sonnet 5**, the
audit and the preparation for Bioconductor with **Claude Opus 5**, and the
changes from the Bioconductor review with **Claude Fable 5.1**.

Concretely, of the code in this repository: `R/containers.R` and everything
under `data-raw/` were written by the model; `R/brs.R` and the test suite
were rewritten from a much smaller first version, and so were the vignette
and this README. Every commit with model-written code carries a
`Co-Authored-By: Claude ...` trailer, so the attribution is visible per commit
and can be checked with

```sh
git log --format=%b | grep Co-Authored-By | sort | uniq -c
```

The work is the audit described in `NEWS.md`: the classifier's reference
groups, score formula and eight silent failure modes were corrected against
Agrawal et al. (2014), and the results were validated against TCGA-THCA. The
findings of an independent review of that work, run on a second machine, are
recorded in the commit `Fix what a macOS audit of this package found`.

No code was copied from another software project. Under Anthropic's terms the
authors own the output, so the contributed code is redistributable under this
package's MIT licence.

Assisted-by: Claude Sonnet 5, Claude Opus 5 and Claude Fable 5.1 (Anthropic)

## Citation

To cite the software, use the Zenodo record of the package (concept DOI, which
always points to the latest version):
[doi:10.5281/zenodo.23271086](https://doi.org/10.5281/zenodo.23271086). The
study that defined the score is Agrawal et al. (2014), linked above. Run
`citation("thyroidBRS")` in R for both entries.

## License

MIT © Camila Barrios-Cervantes, Daniel Pérez-Calixto and Hugo Tovar
