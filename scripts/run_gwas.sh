#!/usr/bin/env bash
# Public, code-only reconstruction of OPSMILE GWAS association workflow.
# All genotype data, sample-level covariates, phenotype files and outputs MUST
# stay on a secure, approved system. This file contains no study data.
set -euo pipefail
usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_gwas.sh --bfile /secure/approved/geno \
    --clinical /secure/approved/analysis_manifest.csv \
    --outdir /secure/results [--maf 0.05] [--threads 4] [--remove /secure/approved/remove.txt]

Clinical manifest requires FID,IID,case,subtype; case 0=control, 1=case;
subtype 1..7 or 99 or blank, following approved clinical adjudication.
Run ONLY using authorized datasets. All outputs contain potentially sensitive data.
Software: PLINK 1.9, Rscript and data.table/ggplot2.
EOF
}
BFILE=''; CLINICAL=''; OUTDIR=''; REMOVE=''; MAF='0.05'; THREADS='4'
while (($#)); do
 case "$1" in
  --bfile) BFILE="$2";shift 2;;
  --clinical) CLINICAL="$2";shift 2;;
  --outdir) OUTDIR="$2";shift 2;;
  --remove) REMOVE="$2";shift 2;;
  --maf) MAF="$2";shift 2;;
  --threads) THREADS="$2";shift 2;;
  --help|-h) usage; exit 0;;
  *) echo "Unknown argument: $1" >&2;usage >&2;exit 2;;
 esac
done
[[ -n "$BFILE" && -n "$CLINICAL" && -n "$OUTDIR" ]] || { usage >&2; exit 2; }
for ext in bed bim fam; do [[ -r "${BFILE}.${ext}" ]] || { echo "Missing ${BFILE}.${ext}" >&2;exit 1; }; done
[[ -r "$CLINICAL" ]] || { echo "Missing protected clinical manifest" >&2;exit 1; }
[[ -z "$REMOVE" || -r "$REMOVE" ]] || { echo "Missing approved removal list" >&2;exit 1; }
command -v plink >/dev/null || { echo 'PLINK 1.9 required' >&2;exit 1; }
command -v Rscript >/dev/null || { echo 'Rscript required' >&2;exit 1; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
umask 077
mkdir -p "$OUTDIR"/{qc,pca,phenotypes,gwas,summary}
OUTDIR="$(cd "$OUTDIR" && pwd)"
# INPUT REQUIREMENT: --bfile has already undergone approved individual missingness,
# sex concordance review, relatedness/duplicate adjudication and phenotype reconciliation.
# This public pipeline must NOT guess which study participants to remove.
qc=(--bfile "$BFILE" --geno 0.05 --maf "$MAF" --hwe 0.0001
    --chr 1-22 --make-bed --threads "$THREADS" --out "$OUTDIR/qc/autosomal_qc")
if [[ -n "$REMOVE" ]]; then qc+=(--remove "$REMOVE"); fi
plink "${qc[@]}"
QC="$OUTDIR/qc/autosomal_qc"
# Historical LD-pruning parameters from the documentation: 200 50 0.2
plink --bfile "$QC" --indep-pairwise 200 50 0.2 \
  --threads "$THREADS" --out "$OUTDIR/pca/prune"
plink --bfile "$QC" --extract "$OUTDIR/pca/prune.prune.in" \
  --pca 5 --threads "$THREADS" --out "$OUTDIR/pca/pcs"
Rscript "$SCRIPT_DIR/gwas_analysis.R" prepare \
  --clinical "$CLINICAL" --fam "${QC}.fam" \
  --pca "$OUTDIR/pca/pcs.eigenvec" --outdir "$OUTDIR/phenotypes"
for n in 3 5; do
 covar="$OUTDIR/phenotypes/covar_PC${n}.tsv"
 covar_names="SEX";for ((j=1;j<=n;j++)); do covar_names+=",PC$j";done
 for subtype in ALL CL CLP CLPM CP; do
  plink --bfile "$QC" \
   --pheno "$OUTDIR/phenotypes/pheno_${subtype}.tsv" --pheno-name PHENO \
   --covar "$covar" --covar-name "$covar_names" \
   --logistic hide-covar --threads "$THREADS" \
   --out "$OUTDIR/gwas/${subtype}_PC${n}"
 done
done
Rscript "$SCRIPT_DIR/gwas_analysis.R" report \
 --gwas-dir "$OUTDIR/gwas" --outdir "$OUTDIR/summary"
echo "Private results: $OUTDIR"
