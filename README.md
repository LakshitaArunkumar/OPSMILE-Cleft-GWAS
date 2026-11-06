# Orofacial cleft GWAS: QC, population structure, and subtype associations

**Research portfolio / public code reconstruction.** Work conducted with Dr. John Shaffer at the University of Pittsburgh. The public repository documents the *analytical methods and coding skills* used to investigate genetic associations across cleft subtypes. **No individual-level OPSMILE data or real association outputs are included.** The original analyses were performed in a controlled research environment; this code is a portable reconstruction and has **not** been verified against the original study dataset.

## Biological question
Do orofacial cleft subtypes have overlapping or distinct genetic associations?

## What I worked on
- Genotype/sample QC and investigation of missingness, sex discordance, relatedness and duplicates (manual exclusions were project-specific and are *not automated here*).
- Autosomal LD pruning and PCA for population structure (`--indep-pairwise 200 50 0.2` in project documentation).
- Logistic additive GWAS adjusted for sex and either PC1–PC3 or PC1–PC5, comparing five cleft phenotype definitions.
- GWAS Manhattan/QQ plots, genomic inflation, and variant-level follow-up.
- Comparison to previously reported cleft-associated SNPs, proxies from LDlinkR/1000 Genomes, and exploratory binomial enrichment testing.
- Regional visualization (LocusZoom) and functional mapping (FUMA) performed outside this public script.

## What is *actually runnable* in this public repository
`scripts/run_gwas.sh` takes **approved**, preprocessed private PLINK files and a de-identified, curated phenotype manifest as local inputs. It performs autosomal marker QC, LD pruning, PCA, five subtype association tests, and PC3/PC5 sensitivity comparisons. `scripts/gwas_analysis.R` prepares PLINK phenotype/covariate files and renders summary plots. This does **not** redo the original clinical adjudication, sample removal decisions, external LDlinkR queries or FUMA.

| Code | Phenotype | Cases (curated subtype codes) |
|---|---|---|
| ALL | Any cleft | 1–7, 99, or missing for cases |
| CL | Cleft lip only | 1, 2, 6 |
| CLP | Cleft lip + palate | 3, 4, 7 |
| CLPM | Cleft lip with/without palate | 1, 2, 3, 4, 6, 7 |
| CP | Cleft palate only | 5 |

All comparisons use controls with `case=0` and no recorded cleft subtype. **The input manifest must already have undergone authorized review and resolution of ambiguous coding.**

## Software & example
Requirements: PLINK 1.9, Bash, R with `data.table` and `ggplot2`. Run on your institutionally approved compute system; do not run real samples in a public GitHub Actions workflow.

```bash
bash scripts/run_gwas.sh \
  --bfile /secure/approved/autosomal_genotypes \
  --clinical /secure/approved/curated_manifest.csv \
  --outdir /secure/private_results \
  --maf 0.05 --threads 4
```

Private manifest schema (illustrative only; **not real people**): `FID,IID,case,subtype`. A sample template is in `examples/clinical_manifest_SYNTHETIC.csv`. `case` must be 0 or 1 and subtype must follow the authorized source dictionary. Phenotypes are converted to PLINK `1=control`, `2=case`. To investigate the documented sensitivity to MAF cutoff, run separate, labeled analyses using `--maf 0.05` and `--maf 0.01` rather than combining results.

## Historical analyses not independently reproducible here
The project documentation describes SNP QC, sex discrepancies, parent–offspring/duplicate identification, phenotype-data reconciliation, PCA by country, suggestive association results, candidate SNPs/LD proxies, LocusZoom, and FUMA. It also describes an exploratory comparison with published cleft-associated loci. These steps are acknowledged as work undertaken in the supervised project, **not claimed as fully implemented in the public scripts**. Exact sample counts and historical results are deliberately not republished without study authorization.

## Interpretation
The original research did not establish genome-wide significant associations in the recorded analyses. Nominal and suggestive associations and candidate SNP enrichment do not independently establish causal variants. The cohort size, phenotype classification, ancestry structure and multiple testing affect interpretation. Models with 3 versus 5 PCs should be evaluated in light of study design rather than selected solely by genomic inflation.

## Privacy, ethics & ownership
- **Never upload:** real `.bed/.bim/.fam`, person-level traits/IDs, sample or family maps, phenotype spreadsheets, curated removal lists, imputation files, PLINK logs, result `.assoc.*` files, original presentations with participant tables, or original `work.zip`.
- **No real data are needed to read the code.** Synthetic placeholders are explicitly labeled as synthetic and **are not results**.
- Obtain PI/institutional approval before sharing any original study-derived data, figures, statistics, manuscripts, or collaborators' materials.
- No license is assigned until code/data ownership and sharing permissions are confirmed.

**Status:** educational/portfolio code reconstruction; not an independently executable reproduction of the protected original project.
