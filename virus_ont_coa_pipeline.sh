#!/usr/bin/env bash
# =============================================================================
# virus_ont_coa_pipeline.sh
# Wrapper pipeline: ONT virus minority variant calling + CoA visualisatie
# =============================================================================
#
# BESCHRIJVING
#   Combineert variant calling en visualisatie in één aanroep.
#   Stap 1: virus_ont_variants_coa.sh   — alignment, VarScan2, tabellen
#   Stap 2: virus_ont_variants_coa_viz.sh — coverage plots, SNP snippets
#
# GEBRUIK
#   bash virus_ont_coa_pipeline.sh -r <ref.fa> [opties]
#
# VERPLICHT
#   -r FILE   Referentie FASTA
#
# OPTIONEEL — VARIANT CALLING
#   -g FILE   GFF3 annotatie (voor gene/AA info en gen grenzen in plots)
#   -i DIR    Map met Nanopore FASTQ reads (default: 01_basecalled)
#   -o DIR    Output map (default: results)
#   -t INT    Aantal threads voor minimap2/samtools (default: 8)
#   -d INT    Minimale read depth voor variant calling (default: 5)
#   -f FLOAT  Minimale allele frequentie variant calling (default: 0.05)
#             ONT R10.4 SUP: 5% is aanbevolen ondergrens bij >500x coverage
#   -m FLOAT  Minimale AF voor minorities/CoA tabel (default: 0.05)
#   -q FLOAT  VarScan2 p-waarde drempel (default: 0.001)
#   -D INT    Max read depth per positie mpileup (default: 50000)
#
# OPTIONEEL — VISUALISATIE
#   -c INT    Context bp links/rechts van SNP in snippets (default: 10)
#   -V        Sla visualisatie over (alleen variant calling)
#
# OUTPUT
#   <outdir>/alignments/    Gesorteerde BAM bestanden
#   <outdir>/variants/      Gefilterde VCF bestanden (raw + CSQ geannoteerd)
#   <outdir>/tables/        TSV tabellen per sample:
#     *.variants.tsv          Alle varianten >= -f
#     *.variants.extended.tsv Uitgebreide CSQ annotatie tabel
#     *.minorities.tsv        Minority variants >= -m
#     *.coa.tsv               CoA rapport tabel
#     *.varscan_qc.tsv        AF distributie voor assay QC
#   <outdir>/viz/           Visualisaties per sample:
#     *.html                  Interactief Plotly rapport
#     *.pdf                   PDF voor CoA rapportage
#   pipeline_run_*.log       Pipeline run log met tool versies en timings
#
# REPO STRUCTUUR
#   virus_ont_coa/
#   ├── run_pipeline.sh               ← dit script
#   ├── scripts/
#   │   ├── virus_ont_variants_coa.sh
#   │   └── virus_ont_variants_coa_viz.sh
#   ├── references/                   ← niet in git, of git-lfs
#   ├── envs/
#   │   ├── coa_env.yml
#   │   └── coa_viz_env.yml
#   └── README.md
#
# INSTALLATIE
#   git clone <repo> virus_ont_coa && cd virus_ont_coa
#   mamba env create -f envs/coa_env.yml
#   mamba env create -f envs/coa_viz_env.yml
#
# VEREISTEN
#   Host    : lelycompute-02.wur.nl
#   Envs    : conda activate minimap2 (variant calling)
#             conda activate coa_viz  (visualisatie)
#   Tools   : minimap2, samtools, bcftools, varscan, tabix, bgzip, python3
#             pysam, plotly, pandas, biopython, kaleido
#
# VOORBEELDEN
#   # Standaard run met GFF3 annotatie:
#   bash virus_ont_coa_pipeline.sh \
#       -r references/KC776174.fa \
#       -g references/KC776174.gff3 \
#       -i 01_basecalled -o results -t 16
#
#   # Alleen variant calling, geen visualisatie:
#   bash virus_ont_coa_pipeline.sh \
#       -r references/KC776174.fa -g references/KC776174.gff3 -V
#
#   # Lagere AF drempel met meer context in snippets:
#   bash virus_ont_coa_pipeline.sh \
#       -r references/KC776174.fa -g references/KC776174.gff3 \
#       -f 0.03 -m 0.03 -c 15
#
# AUTEUR  : WBVR Bioinformatics (harde004)
# VERSIE  : 1.0
# DATUM   : 2026-06-10

set -Eeuo pipefail

############################################
# Scriptlocaties
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VARIANTS_SCRIPT="$SCRIPT_DIR/scripts/virus_ont_variants_coa.sh"
VIZ_SCRIPT="$SCRIPT_DIR/scripts/virus_ont_variants_coa_viz.sh"

############################################
# Usage
############################################
usage(){
    awk 'NR==1{next} /^#!/{next} /^[^#]/{exit} {sub(/^# ?/,""); print}' "$0"
    exit "${1:-0}"
}

############################################
# Defaults
############################################
REF=""
GFF_IN=""
READS_DIR="01_basecalled"
OUT="results"
THREADS=8
MINREADS=5
MINAF=0.05
MINORITY_AF=0.05
VARSCAN_PVAL=0.001
MAXDEPTH=50000
CONTEXT=10
SKIP_VIZ=0

while getopts ":r:g:i:o:t:d:f:m:q:D:c:Vh" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        g) GFF_IN="$OPTARG" ;;
        i) READS_DIR="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) MINREADS="$OPTARG" ;;
        f) MINAF="$OPTARG" ;;
        m) MINORITY_AF="$OPTARG" ;;
        q) VARSCAN_PVAL="$OPTARG" ;;
        D) MAXDEPTH="$OPTARG" ;;
        c) CONTEXT="$OPTARG" ;;
        V) SKIP_VIZ=1 ;;
        h) usage 0 ;;
        :) echo "ERROR: Optie -$OPTARG vereist een argument." >&2; usage 1 ;;
        \?) echo "ERROR: Onbekende optie -$OPTARG." >&2; usage 1 ;;
    esac
done

############################################
# Validatie
############################################
errors=0
[[ -z "$REF" ]]                        && { echo "ERROR: Referentie FASTA verplicht (-r)." >&2; errors=1; }
[[ -n "$REF" && ! -f "$REF" ]]         && { echo "ERROR: Referentie niet gevonden: $REF" >&2; errors=1; }
[[ -n "$GFF_IN" && ! -f "$GFF_IN" ]]   && { echo "ERROR: GFF3 niet gevonden: $GFF_IN" >&2; errors=1; }
[[ ! -d "$READS_DIR" ]]                && { echo "ERROR: Reads map niet gevonden: $READS_DIR" >&2; errors=1; }
[[ ! -f "$VARIANTS_SCRIPT" ]]          && { echo "ERROR: Script niet gevonden: $VARIANTS_SCRIPT" >&2; echo "       Verwacht in: $SCRIPT_DIR/scripts/" >&2; errors=1; }
[[ $SKIP_VIZ -eq 0 && ! -f "$VIZ_SCRIPT" ]] && { echo "ERROR: Script niet gevonden: $VIZ_SCRIPT" >&2; echo "       Verwacht in: $SCRIPT_DIR/scripts/" >&2; errors=1; }
[[ $errors -gt 0 ]] && usage 1

############################################
# Helpers
############################################
ts(){ date +"%F %T"; }
say(){ echo "[$(ts)] $*" >&2; }
elapsed(){ echo $(( $(date +%s) - $1 )); }

PIPELINE_START=$(date +%s)

say "================================================================"
say " virus_ont_coa_pipeline.sh v1.0 — WBVR Bioinformatics"
say "================================================================"
say "  Referentie : $REF"
say "  Reads map  : $READS_DIR"
say "  Output     : $OUT"
say "  Threads    : $THREADS"
say "  Min AF     : $MINAF  |  Minority AF: $MINORITY_AF"
[[ -n "$GFF_IN" ]] && say "  GFF3       : $GFF_IN" || say "  GFF3       : niet opgegeven"
[[ $SKIP_VIZ -eq 1 ]] && say "  Visualisatie: overgeslagen (-V)" || say "  Visualisatie: aan"

############################################
# STAP 1: Variant calling
############################################
say ""
say "════ STAP 1: Variant calling ════"

VARIANTS_CMD=(
    bash "$VARIANTS_SCRIPT"
    -r "$REF"
    -i "$READS_DIR"
    -o "$OUT"
    -t "$THREADS"
    -d "$MINREADS"
    -f "$MINAF"
    -m "$MINORITY_AF"
    -q "$VARSCAN_PVAL"
    -D "$MAXDEPTH"
)
[[ -n "$GFF_IN" ]] && VARIANTS_CMD+=(-g "$GFF_IN")

T1=$(date +%s)
"${VARIANTS_CMD[@]}"
say "════ Variant calling klaar in $(elapsed $T1)s ════"

############################################
# STAP 2: Visualisatie
############################################
if [[ $SKIP_VIZ -eq 0 ]]; then
    say ""
    say "════ STAP 2: Visualisatie ════"

    VIZ_CMD=(
        bash "$VIZ_SCRIPT"
        -o "$OUT"
        -r "$REF"
        -m "$MINORITY_AF"
        -c "$CONTEXT"
    )
    [[ -n "$GFF_IN" ]] && VIZ_CMD+=(-g "$GFF_IN")

    T2=$(date +%s)
    "${VIZ_CMD[@]}"
    say "════ Visualisatie klaar in $(elapsed $T2)s ════"
fi

############################################
# Samenvatting
############################################
TOTAL=$(elapsed $PIPELINE_START)
say ""
say "================================================================"
say " PIPELINE KLAAR in ${TOTAL}s"
say "================================================================"
say " Tabellen : $OUT/tables/"
[[ $SKIP_VIZ -eq 0 ]] && say " Rapporten: $OUT/viz/"
say " IGV      : $OUT/alignments/*.sorted.bam + $REF"
say "================================================================"
