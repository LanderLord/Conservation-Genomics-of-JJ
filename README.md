# Genomic vulnerability of brown bears (*Ursus arctos*) on the Tibetan Plateau and in the Himalaya

Code walkthroughs for the analyses in our study: variant functional annotation,
genotype–environment association, genome-wide selection scans, gradient-forest-based
genetic offset, and landscape connectivity (isolation by resistance).

Each `*_methods.md` file is a **self-contained, commented walkthrough**: the R/bash code
blocks can be pasted into an interactive session or split into scripts, and the comments
follow the order of the Methods section of the paper. `psmc_plot_avg.pl` is the only
standalone executable (it post-processes PSMC demographic results used alongside the
landscape-genomic analyses).

## Workflow

```
raw VCF / genotype matrix
  │
  ├─ 01  snpEff annotation + per-sample genotype counts     snpEff_annotation_methods.md
  │
  ├─ 02  genotype–climate association (LFMM2, RDA)
  │      + iHS selection scan                               GEA_LFMM_RDA_iHS_methods.md
  │
  ├─ 03  gradient forest + local/forward/backward
  │      genetic offset                                     GF_genetic_offset_methods.md
  │
  ├─ 04  resistance surfaces (ResistanceGA) + MLPE          landscape_connectivity_methods.md
  │
  └─ 05  PSMC trajectories + cross-sample mean curve        psmc_plot_avg.pl
```

## Files

| File                                | Language   | What it does                                                 | Main inputs                                                  | Main outputs                                                 |
| ----------------------------------- | ---------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `snpEff_annotation_methods.md`      | bash / awk | Annotates the merged VCF with snpEff (`Ursus_arctos` database), subsets four functional classes (synonymous, missense, loss-of-function, HIGH impact) and counts homozygous / heterozygous genotypes and the proportion of homozygous alleles per sample | annotated VCF                                                | `{Syn,NONSYN,LOF,HIGH}.vcf.gz`, `*_analysis.csv`, `snpEff_summary.csv` |
| `GEA_LFMM_RDA_iHS_methods.md`       | R / bash   | (i) LFMM2 (`LEA`, K = 2, q-value FDR) run separately for five bioclimatic variables, then a shell step that merges results and keeps SNPs significant in ≥ 2 variables; (ii) RDA (`vegan`) with a 3-SD loading cutoff; (iii) per-chromosome iHS (`rehh`), standardised genome-wide with 0.5 % / 99.5 % cutoffs | MAF matrix (SNP × sample), sample climate table, imputed VCF | `lfmm_*.txt`, `snp_env_count.tsv`, `multi_env_candidate.snp`, `rda_candidate_SNPs.csv`, `ihs_standardized_results.csv` |
| `GF_genetic_offset_methods.md`      | R          | Trains a gradient forest model on allele frequencies, predicts the transformed genetic composition on ~ 50 000 grid cells for current and future climate, and derives local, forward and backward genetic offset plus the corresponding geographic distance | allele-frequency matrix, sampling coordinates, current/future climate rasters | `gf_model.rds`, `GF_{Local,Forward,Backward}_Offset.csv`, `GF_*_Offset.tif` |
| `landscape_connectivity_methods.md` | R          | Optimises resistance surfaces with `ResistanceGA` (CIRCUITSCAPE v4 via Julia), then fits MLPE mixed models (geographic, single-variable, temperature, humidity, annual variability, elevation, climate, climate + elevation) and compares them by AICc, ΔAICc, Akaike weight and marginal/conditional R² | sampling coordinates, genetic distance matrix, environmental ASCII layers | optimised resistance matrices, `MLPE_AICc_AkaikeWeight_R2_results.tsv` |
| `psmc_plot_avg.pl`                  | Perl       | Single-file PSMC post-processing: parses `.psmc` output, rescales time (years) and *N*e (× 10⁴) with a given mutation rate and generation time, rebuilds every trajectory as a step function, evaluates all trajectories on the shared breakpoint grid, computes the equally weighted mean, and plots individuals (thin lines) plus the mean (thick red line) with gnuplot | `.psmc` file(s), `-M` sample names                           | `*.eps`, `*_all_coords.txt`, `*_avg_coords.txt`, `*.gp`      |

## Requirements

- **R ≥ 4.0** with `LEA` (3.4.0), `qvalue`, `vegan` (2.5-7), `rehh`, `gradientForest` (0.1-24),
  `raster`, `dplyr`, `fields`, `geosphere`, `foreach`, `doParallel`, `ResistanceGA` (4.2-2),
  `lme4`, `Matrix`
- **Julia 1.7** (required by the CIRCUITSCAPE optimiser inside `ResistanceGA`)
- **Java 21** + `snpEff` with the *Ursus arctos* database; `bcftools`, `bgzip`, `tabix`; Beagle 5.4
- **Perl 5** (core modules only); `gnuplot` for `psmc_plot_avg.pl`

## Usage

```bash
# 1–4: open each .md file and run its code blocks in order (R console / shell),
#      after editing the placeholder paths to your own working directory

# 5: PSMC — one command instead of three scripts
perl psmc_plot_avg.pl -g 10 -u 1.0e-08 -M 'QH9,SC1,SC4' -w 1 -W 3 -n 0 \
     -T "Ursus arctos" Hima_combine_plot Hima_combine.psmc
```

`psmc_plot_avg.pl` without arguments prints the full option list; `-n 0` keeps the
best-RI iteration instead of a fixed EM round (default `-n 20`).

## Notes and caveats

- **Placeholder paths**: file paths such as `D:/...` and `/public/home/...` in the
  walkthroughs refer to our working directories and must be edited before running.
- **Gradient forest**: `maxLevel = 1` was used because the number of populations is small;
  forward offset is computed over all grid cells, i.e. without restricting target cells to
  the suitable niche.
- **ResistanceGA**: `run_test = TRUE`, `run = 5`, `maxiter = 10` in the walkthrough are test
  settings; for production runs use `run_test = FALSE` and increase `run` / `maxiter` until
  the optimiser converges.
- **PSMC averaging**: the mean curve is an equally weighted arithmetic mean of *N*e across
  individuals, evaluated on the union of all breakpoints with a piecewise-constant (step)
  rule — no interpolation, no smoothing, no confidence interval.
- Genotypes are coded 0/1/2; the LFMM/RDA walkthroughs transpose the SNP × sample matrix to
  the sample × SNP layout required by `lfmm2` and `rda`.



## Contact

Questions and bug reports: please open an issue.
