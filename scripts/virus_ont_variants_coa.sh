#!/usr/bin/env bash
# =============================================================================
# virus_ont_variants_coa.sh
# Generiek virus Nanopore minority variant calling pipeline voor CoA rapportage
# =============================================================================
#
# BESCHRIJVING
#   Verwerkt Oxford Nanopore FASTQ reads naar variant calling tabellen geschikt
#   voor Certificate of Analysis (CoA) rapportage van virus batches.
#   Pipeline: minimap2 alignment → VarScan2 minority variant calling →
#             bcftools csq annotatie → TSV tabellen
#
# GEBRUIK
#   bash virus_ont_variants_coa.sh -r <referentie.fa> [opties]
#
# VERPLICHT
#   -r FILE   Referentie FASTA (Windows regeleindes verwijderen met: sed -i 's/\r//' ref.fa)
#
# OPTIONEEL
#   -g FILE   GFF3 annotatie voor gene/AA informatie in tabellen
#   -i DIR    Map met Nanopore FASTQ reads (default: 01_basecalled)
#   -o DIR    Output map (default: results)
#   -t INT    Aantal threads (default: 8)
#   -d INT    Minimale read depth (default: 5)
#   -f FLOAT  Minimale allele frequentie voor variant calling (default: 0.05)
#             ONT R10.4 SUP heeft ~0.5-1% systematische error rate bij hoge coverage.
#             Bij coverage >1000x is 5% de aanbevolen ondergrens. De QC tabel toont
#             altijd de volledige distributie vanaf 2% voor assay kwaliteitsbepaling.
#   -m FLOAT  Minimale allele frequentie voor minorities/CoA tabel (default: 0.05)
#   -q FLOAT  VarScan2 p-waarde drempel (default: 0.001)
#   -D INT    Max read depth per positie voor mpileup (default: 50000)
#             Let op: bij hoge ONT coverage (>500x) zijn zelfs 2% minorities
#             statistisch significant (p<E-20). P-waarde filtering heeft daardoor
#             weinig effect; gebruik -f voor ruisreductie.
#   -h        Toon deze helptext
#
# OUTPUT (per sample in <outdir>/tables/)
#   *.variants.tsv          Alle varianten >= -f, gesorteerd op positie
#   *.variants.extended.tsv Uitgebreide tabel met CSQ consequence annotatie
#   *.minorities.tsv        Minority variants >= -m met gen/AA info
#   *.coa.tsv               CoA tabel: leesbare gen namen, gestandaardiseerde AA notatie
#   *.varscan_qc.tsv        AF distributie 2-50% voor assay ruis karakterisatie
#
# VARIANT CALLING INSTELLINGEN (VarScan2 mpileup2snp)
#   -B              Geen BAQ correctie (aanbevolen voor ONT, voorkomt vals-negatieve calls
#                   bij posities nabij indels)
#   -Q 5            Lage base quality filter (ONT Q-scores zijn niet Illumina-equivalent)
#   --min-MQ 20     Mapping quality filter (verwijdert multi-mapping reads)
#   --strand-filter 0  Geen strand bias filter (niet van toepassing bij ONT)
#
# ANNOTATIE NOTATIE
#   AA changes worden gerapporteerd in standaard eenletter notatie: ref+pos+alt
#   Voorbeeld: G15D (Glycine positie 15 → Aspartaat)
#   Stop codons worden als X gerapporteerd (bijv. E336X, Q345X)
#   Synoniemen worden als "Silent mutation" gerapporteerd
#   Intergenic varianten als "Untranslated"
#   MNV (multi-nucleotide variant) annotaties worden automatisch geresolveerd
#   vanuit bcftools csq @pos pointers.
#
# VEREISTE TOOLS
#   minimap2, samtools, bcftools, varscan, tabix, bgzip, python3, awk, bc
#   Aanbevolen: conda environment 'minimap2' op lelycompute-02.wur.nl
#
# AUTEUR  : WBVR Bioinformatics (harde004)
# VERSIE  : 3.1
# DATUM   : 2026-06-09

set -Eeuo pipefail

############################################
# Hostname check
############################################
REQUIRED_HOST="lelycompute-02"
CURRENT_HOST=$(hostname -s)
if [[ "$CURRENT_HOST" != "$REQUIRED_HOST" ]]; then
    echo "ERROR: Dit script moet draaien op ${REQUIRED_HOST}.wur.nl" >&2
    echo "       Huidige host: ${CURRENT_HOST}" >&2
    echo "       Reden: conda environment 'minimap2' met alle vereiste tools" >&2
    echo "              is alleen beschikbaar op lelycompute-02." >&2
    exit 1
fi

############################################
# Conda environment activatie
############################################
CONDA_ENV="minimap2"
if [[ -z "${CONDA_DEFAULT_ENV:-}" || "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; then
    if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "${HOME}/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        source "/opt/conda/etc/profile.d/conda.sh"
    else
        echo "ERROR: conda niet gevonden. Activeer handmatig: conda activate ${CONDA_ENV}" >&2
        exit 1
    fi
    conda activate "$CONDA_ENV"
fi

############################################
# Usage / argument parsing
############################################
usage(){
    awk 'NR==1{next} /^#!/{next} /^[^#]/{exit} {sub(/^# ?/,""); print}' "$0"
    exit "${1:-0}"
}

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

while getopts ":r:g:i:o:t:d:f:m:q:D:h" opt; do
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
        h) usage 0 ;;
        :) echo "ERROR: Optie -$OPTARG vereist een argument." >&2; usage 1 ;;
        \?) echo "ERROR: Onbekende optie -$OPTARG." >&2; usage 1 ;;
    esac
done

# Validatie
errors=0
[[ -z "$REF" ]]              && { echo "ERROR: Referentie FASTA is verplicht (-r)." >&2; errors=1; }
[[ -n "$REF" && ! -f "$REF" ]] && { echo "ERROR: Referentie FASTA niet gevonden: $REF" >&2; errors=1; }
[[ -n "$GFF_IN" && ! -f "$GFF_IN" ]] && { echo "ERROR: GFF3 niet gevonden: $GFF_IN" >&2; errors=1; }
[[ ! -d "$READS_DIR" ]]      && { echo "ERROR: Reads map niet gevonden: $READS_DIR" >&2; errors=1; }
[[ $errors -gt 0 ]] && usage 1

# CSQ veld indexen (bcftools csq output: consequence|gene|transcript|biotype|strand|aa|dna)
CSQ_FIELD_GENE=2
CSQ_FIELD_AA=6

############################################
# Helpers
############################################
ts(){ date +"%F %T"; }
say(){ echo "[$(ts)] $*" >&2; }
die(){ echo "[$(ts)] ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "$1 niet gevonden in PATH."; }
hdr_id(){ head -1 "$1" | cut -d ' ' -f1 | sed 's/^>//'; }
elapsed(){ echo $(( $(date +%s) - $1 )) ; }

sname(){
    bn=$(basename "$1")
    bn=${bn%%.fastq.gz}; bn=${bn%%.fq.gz}; bn=${bn%%.fastq}; bn=${bn%%.fq}
    echo "$bn"
}

############################################
# Checks & versie logging
############################################
say "================================================================"
say " virus_ont_variants_coa.sh v3.1 - WBVR Bioinformatics"
say "================================================================"
say "STEP 0: Initialisatie en tool checks"
for t in minimap2 samtools bcftools varscan tabix bgzip python3 awk sed grep sort bc; do need "$t"; done

# Log tool versies naar main log
MAIN_LOG="$OUT/../pipeline_run_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$MAIN_LOG")" 2>/dev/null || MAIN_LOG="/tmp/pipeline_run_$(date +%Y%m%d_%H%M%S).log"

{
    echo "================================================================"
    echo " Pipeline: virus_ont_variants_coa.sh v3.1"
    echo " Datum   : $(date)"
    echo " Host    : $(hostname)"
    echo " User    : $(whoami)"
    echo " Werkmap : $(pwd)"
    echo "================================================================"
    echo ""
    echo "--- Tool versies ---"
    echo "minimap2  : $(minimap2 --version 2>&1 | head -1)"
    echo "samtools  : $(samtools --version 2>&1 | head -1)"
    echo "bcftools  : $(bcftools --version 2>&1 | head -1)"
    echo "varscan   : $(varscan 2>&1 | head -1)"
    echo "python3   : $(python3 --version 2>&1)"
    echo "bgzip     : $(bgzip --version 2>&1 | head -1)"
    echo "tabix     : $(tabix --version 2>&1 | head -1)"
    echo ""
    echo "--- Parameters ---"
    echo "Referentie  : $REF"
    echo "Reads map   : $READS_DIR"
    echo "Output      : $OUT"
    echo "Threads     : $THREADS"
    echo "Min DP      : $MINREADS"
    echo "Min AF      : $MINAF  (variant calling filter)"
    echo "Minority AF : $MINORITY_AF  (CoA/minorities drempel)"
    echo "VarScan p   : $VARSCAN_PVAL"
    [[ -n "$GFF_IN" ]] && echo "GFF3        : $GFF_IN" || echo "GFF3        : niet opgegeven"
    echo ""
} > "$MAIN_LOG"

say "  Referentie  : $REF"
say "  Reads map   : $READS_DIR"
say "  Output      : $OUT"
say "  Threads     : $THREADS"
say "  Min AF      : $MINAF  |  Minority AF: $MINORITY_AF  |  VarScan p: $VARSCAN_PVAL"
[[ -n "$GFF_IN" ]] && say "  GFF3        : $GFF_IN" || say "  GFF3        : niet opgegeven (geen CSQ annotatie)"

############################################
# Init
############################################
PIPELINE_START=$(date +%s)
mkdir -p "$OUT"
find "$OUT" -mindepth 1 -exec rm -rf -- {} + 2>/dev/null || true
mkdir -p "$OUT"/{alignments,variants,tables,logs,temp}
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

############################################
# STAP 1: GFF verwerken
############################################
CSQ_GFF=""
CSQ_OK=0

if [[ -n "$GFF_IN" ]]; then
    say "STEP 1: GFF-verwerking: normalisatie en synthese indien nodig"
    REFID="$(hdr_id "$REF")"
    GFF_NORM="$OUT/temp/input.relabel.gff3"

    say " - Normaliseren van GFF seqid naar FASTA ID: $REFID"
    AWK_SCRIPT='BEGIN{OFS="\t"} /^#/ {print; next} {$1=id; print}'
    if [[ "$GFF_IN" =~ \.gz$ ]]; then
        gzip -cd "$GFF_IN" | awk -v id="$REFID" "$AWK_SCRIPT" > "$GFF_NORM"
    else
        awk -v id="$REFID" "$AWK_SCRIPT" "$GFF_IN" > "$GFF_NORM"
    fi

    if ! grep -Eq $'\t(mRNA|transcript)\t' "$GFF_NORM"; then
        say " - Geen mRNA features: synthetiseren van gene/mRNA/exon uit CDS"
        CSQ_GFF="$OUT/temp/input.csqready.gff3"
        awk -v OFS="\t" -v refid="$REFID" '
            BEGIN{ print "##gff-version 3" }
            function get_attr(h,key, i,n,a){
                n=split(h,a,";");
                for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]==key) return kv[2]}
                return ""
            }
            {
                c=$1; s=$2; f=$3; L=$4; R=$5; sc=$6; st=$7; ph=$8; at=$9;
                if (L ~ /^[0-9]+$/ && R ~ /^[0-9]+$/) {
                    if (f=="gene") {
                        gid=get_attr(at,"ID"); if(gid!=""){ g_strand[gid]=st; g_seen[gid]=1; g_name[gid]=gid }
                    }
                    else if (f=="CDS") {
                        g=get_attr(at,"gene"); if(g=="") g=get_attr(at,"Parent"); if(g=="") g=get_attr(at,"GeneID");
                        if(g==""){ g="gene" ++gcount }
                        cdsL[g][++idx[g]]=L; cdsR[g][idx[g]]=R; cdsSC[g][idx[g]]=(sc==""?"." : sc); cdsPH[g][idx[g]]=(ph==""?"." : ph);
                        cdsST[g]=st;
                        if(!(g in minL) || L<minL[g]) minL[g]=L;
                        if(!(g in maxR) || R>maxR[g]) maxR[g]=R;
                        g_name[g]=g; cds_seen[g]=1;
                    }
                    else { rest=rest $0 "\n" }
                }
            }
            END{
                for(g in cds_seen){
                    strand = (cdsST[g]!="") ? cdsST[g] : "+";
                    L=minL[g]; R=maxR[g];
                    gene_attrs = "ID=" g ";Name=" g ";gene=" g
                    print refid, "synthetic", "gene", L, R, ".", strand, ".", gene_attrs
                    mr = g ".t1"
                    mr_attrs = "ID=" mr ";Parent=" g ";Name=" mr ";gene=" g ";biotype=protein_coding"
                    print refid, "synthetic", "mRNA", L, R, ".", strand, ".", mr_attrs
                    nseg = asorti(cdsL[g], ord_idx)
                    for(i=1;i<=nseg;i++){
                        k = ord_idx[i];
                        eL = cdsL[g][k]; eR = cdsR[g][k];
                        exon_attrs = "ID=" mr ".exon" i ";Parent=" mr
                        print refid, "synthetic", "exon", eL, eR, ".", strand, ".", exon_attrs
                        cds_attrs = "Parent=" mr ";gene=" g
                        phase = cdsPH[g][k]; if(phase=="") phase="."
                        score = cdsSC[g][k]; if(score=="") score="."
                        print refid, "synthetic", "CDS", eL, eR, score, strand, phase, cds_attrs
                    }
                }
                if(rest!="") printf "%s", rest
            }
        ' "$GFF_NORM" > "$CSQ_GFF"

        if ! grep -Eq $'\tCDS\t' "$CSQ_GFF"; then
            say " - WAARSCHUWING: Synthese mislukt. Annotatie wordt overgeslagen."
            CSQ_GFF=""
        else
            CSQ_OK=1
            say " - Synthese succesvol: $CSQ_GFF"
        fi
    else
        CSQ_GFF="$GFF_NORM"
        CSQ_OK=1
        say " - mRNA/transcript features gevonden: $CSQ_GFF"
    fi
else
    say "STEP 1: Geen GFF meegegeven → tabellen zonder gene/AA annotatie"
fi

# Gene-ID → leesbare naam lookup (voor CoA tabel)
GENE_LOOKUP_FILE="$OUT/temp/gene_lookup.tsv"
if [[ -n "$GFF_IN" ]]; then
    grep -v '^#' "$GFF_IN" | awk -F'\t' '$3=="gene"' | cut -f9 \
        | awk -F';' '{
            id=""; name="";
            for(i=1;i<=NF;i++){
                if($i~/^ID=/)   id=substr($i,4)
                if($i~/^Name=/) name=substr($i,6)
            }
            if(name=="") { name=id; sub(/^gene-/,"",name) }
            print id "\t" name " protein"
        }' > "$GENE_LOOKUP_FILE"
else
    touch "$GENE_LOOKUP_FILE"
fi

############################################
# Samples (stap 2 t/m 5)
############################################
shopt -s nullglob
READ_FILES=( "$READS_DIR"/*.fastq.gz "$READS_DIR"/*.fq.gz "$READS_DIR"/*.fastq "$READS_DIR"/*.fq )
TMP=(); for f in "${READ_FILES[@]}"; do [[ -f "$f" ]] && TMP+=("$f"); done; READ_FILES=( "${TMP[@]}" )
[[ ${#READ_FILES[@]} -gt 0 ]] || die "Geen FASTQs gevonden in $READS_DIR"

say "Gevonden: ${#READ_FILES[@]} sample(s) in $READS_DIR"
echo "--- Samples ---" >> "$MAIN_LOG"
for f in "${READ_FILES[@]}"; do echo "  $(basename "$f")" >> "$MAIN_LOG"; done
echo "" >> "$MAIN_LOG"

SAMPLE_NUM=0
for READ_FILE in "${READ_FILES[@]}"; do
    SAMPLE_NUM=$(( SAMPLE_NUM + 1 ))
    SAMPLE="$(sname "$READ_FILE")"
    LOG="$OUT/logs/${SAMPLE}.log"
    SAMPLE_START=$(date +%s)

    say "================================================================"
    say "SAMPLE $SAMPLE_NUM/${#READ_FILES[@]}: $SAMPLE"
    say "================================================================"

    # Header in sample log
    {
        echo "================================================================"
        echo " Sample: $SAMPLE"
        echo " Input : $READ_FILE"
        echo " Start : $(date)"
        echo "================================================================"
    } > "$LOG"

    # --- STAP 2: Alignment ---
    say "STEP 2: [$SAMPLE] minimap2 alignment → samtools sort/index"
    T2=$(date +%s)
    minimap2 -ax map-ont -t "$THREADS" "$REF" "$READ_FILE" 2>> "$LOG" \
        | samtools sort -@ "$THREADS" -o "$OUT/alignments/${SAMPLE}.sorted.bam" - 2>> "$LOG"
    samtools index -@ "$THREADS" "$OUT/alignments/${SAMPLE}.sorted.bam" 2>> "$LOG"

    # Alignment statistieken
    N_READS=$(samtools view -c "$OUT/alignments/${SAMPLE}.sorted.bam" 2>/dev/null || echo "?")
    N_MAPPED=$(samtools view -c -F 4 "$OUT/alignments/${SAMPLE}.sorted.bam" 2>/dev/null || echo "?")
    say "  -> Alignment: $N_MAPPED / $N_READS reads gemapt ($(elapsed $T2)s)"
    echo "Alignment: $N_MAPPED/$N_READS reads gemapt" >> "$LOG"

    BAM_FILE="$OUT/alignments/${SAMPLE}.sorted.bam"
    RAW_VCF="$OUT/temp/${SAMPLE}.varscan.vcf"
    FINAL_VCF="$OUT/variants/${SAMPLE}.raw.vcf.gz"

    # --- STAP 3: VarScan2 variant calling ---
    say "STEP 3: [$SAMPLE] VarScan2 minority variant calling"
    T3=$(date +%s)
    samtools mpileup -Q 5 -B -d "$MAXDEPTH" --min-MQ 20 \
        -f "$REF" "$BAM_FILE" 2>> "$LOG" \
        | varscan mpileup2snp \
            --min-coverage "$MINREADS" \
            --min-var-freq "$MINAF" \
            --min-avg-qual 5 \
            --strand-filter 0 \
            --p-value "$VARSCAN_PVAL" \
            --output-vcf 1 \
        > "$RAW_VCF" 2>> "$LOG"

    bgzip -f "$RAW_VCF"
    tabix -f "${RAW_VCF}.gz"
    mv "${RAW_VCF}.gz"     "$FINAL_VCF"
    mv "${RAW_VCF}.gz.tbi" "$FINAL_VCF.tbi"

    N_VARS=$(bcftools view -H "$FINAL_VCF" | wc -l)
    say "  -> $N_VARS varianten gevonden >= ${MINAF} ($(elapsed $T3)s)"

    CURRENT_SRC="$FINAL_VCF"
    VCF_CSQ_OK=0

    # --- STAP 4: CSQ annotatie + MNV pointer resolutie ---
    if [[ "$CSQ_OK" -eq 1 ]]; then
        say "STEP 4: [$SAMPLE] bcftools csq annotatie"
        T4=$(date +%s)
        CSQ_VCF_RAW="$OUT/temp/${SAMPLE}.csq_raw.vcf.gz"
        CSQ_VCF="$OUT/variants/${SAMPLE}.csq.vcf.gz"

        if bcftools csq -f "$REF" -g "$CSQ_GFF" -p a -c CSQ -O z -o "$CSQ_VCF_RAW" "$CURRENT_SRC" >> "$LOG" 2>&1; then
            tabix -f "$CSQ_VCF_RAW"
            bcftools view "$CSQ_VCF_RAW" > "$OUT/temp/${SAMPLE}.csq_raw.vcf"

            say "  -> CSQ OK. Resolven van MNV @pos pointers..."
            # Python resolver: bcftools csq annoteert MNV codon changes alleen op de eerste
            # positie; latere posities in hetzelfde codon krijgen CSQ=@refpos als pointer.
            # Dit script kopieert de annotatie van de referentiepositie naar de pointers.
            python3 - "$OUT/temp/${SAMPLE}.csq_raw.vcf" "$CSQ_VCF" << 'PYEOF_INLINE'
import sys
in_file, out_file = sys.argv[1], sys.argv[2]
cache = {}
lines = open(in_file).readlines()
# Pass 1: build cache van pos → CSQ annotatie
for line in lines:
    if line.startswith("#"): continue
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8: continue
    if "CSQ=" in fields[7] and "CSQ=@" not in fields[7]:
        for part in fields[7].split(";"):
            if part.startswith("CSQ="): cache[fields[1]] = part[4:]
# Pass 2: resolve @pos pointers en schrijf naar bgzip output
import subprocess
proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=open(out_file,"wb"))
for line in lines:
    if not line.startswith("#"):
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 8 and "CSQ=@" in fields[7]:
            new_info = []
            for part in fields[7].split(";"):
                if part.startswith("CSQ=@"):
                    ref_pos = part[5:]
                    part = "CSQ=" + cache.get(ref_pos, "@" + ref_pos)
                new_info.append(part)
            fields[7] = ";".join(new_info)
            line = "\t".join(fields) + "\n"
    proc.stdin.write(line.encode())
proc.stdin.close(); proc.wait()
PYEOF_INLINE
            tabix -f "$CSQ_VCF"
            rm -f "$OUT/temp/${SAMPLE}.csq_raw.vcf" "$CSQ_VCF_RAW" "$CSQ_VCF_RAW.tbi"
            VCF_CSQ_OK=1
            CURRENT_SRC="$CSQ_VCF"
            N_CSQ=$(bcftools view -H -i 'INFO/CSQ!=""' "$CSQ_VCF" | wc -l)
            say "  -> CSQ annotatie: $N_CSQ / $N_VARS varianten geannoteerd ($(elapsed $T4)s)"
        else
            say "  -> CSQ FAILED. Annotatie overgeslagen (zie $LOG)."
        fi
    fi

    # --- STAP 5: Tabellen genereren ---
    say "STEP 5: [$SAMPLE] tabellen genereren"
    T5=$(date +%s)

    # Gemeenschappelijke awk functies voor VarScan output parsing
    VARSCAN_AWK_FIELDS='
        function get_freq(freq_str) {
            gsub(/%/, "", freq_str); return freq_str+0
        }
        function format_aa(aa_raw, consequence,    p, ref_aa, pos_aa, alt_aa, result) {
            if (aa_raw == "") return ""
            # Formaat "1S>1P" of "15G>15D" (pos+AA > pos+AA) → ref_aa + pos + alt_aa
            if (aa_raw ~ /^[0-9]+[A-Z*]>[0-9]+[A-Z*]$/) {
                split(aa_raw, p, ">")
                ref_aa = substr(p[1], length(p[1]), 1)
                pos_aa = substr(p[1], 1, length(p[1])-1)
                alt_aa = substr(p[2], length(p[2]), 1)
                result = ref_aa pos_aa alt_aa
            } else if (aa_raw ~ /^[A-Z*][0-9]+>[A-Z*][0-9]+$/) {
                # Formaat "L6373>Q6373" (AA+pos > AA+pos) → ref_aa + pos + alt_aa
                split(aa_raw, p, ">")
                ref_aa = substr(p[1], 1, 1)
                pos_aa = substr(p[1], 2)
                alt_aa = substr(p[2], 1, 1)
                result = ref_aa pos_aa alt_aa
            } else {
                result = aa_raw
            }
            # Stop codons als X (HGVS-standaard voor rapportage)
            gsub(/\*/, "X", result)
            return result
        }
    '

    # --- SIMPLE TABLE: alle varianten >= MINAF ---
    SIMPLE_TSV="$OUT/tables/${SAMPLE}.variants.tsv"
    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        FMT_CSQ='%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t%INFO/CSQ\t1\n'
        FMT_NO_CSQ='%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n'
        {
            bcftools query -i 'INFO/CSQ!=""' -f "$FMT_CSQ"    "$CURRENT_SRC"
            bcftools query -i 'INFO/CSQ=""'  -f "$FMT_NO_CSQ" "$CURRENT_SRC"
        } \
        | awk -v OFS="\t" -v G_IDX="$CSQ_FIELD_GENE" -v AA_IDX="$CSQ_FIELD_AA" \
              "$VARSCAN_AWK_FIELDS"'
            {
                pos=$1; ref=$2; alt=$3; dp=$4; freq_raw=$5; csq_str=$6; annotated=$7;
                freq = get_freq(freq_raw)
                gene=""; cDNA=""; aa_raw=""; consequence=""; aa_fmt="";
                if (annotated == "1") {
                    split(csq_str, csq_vals, "|");
                    consequence=csq_vals[1]; gene=csq_vals[G_IDX]; aa_raw=csq_vals[AA_IDX];
                    cDNA = ref pos alt;
                    aa_fmt = (consequence ~ /synonymous/) ? "Silent mutation" : format_aa(aa_raw, consequence)
                } else {
                    gene="Intergenic"; cDNA=ref pos alt;
                }
                if(length(ref)==1&&length(alt)==1){type="SNP";len=1}
                else if(length(ref)<length(alt)){type="INS";len=length(alt)-length(ref)}
                else if(length(ref)>length(alt)){type="DEL";len=length(ref)-length(alt)}
                else{type="VAR";len=length(alt)}
                print type, tolower(ref)pos tolower(alt), dp, len, sprintf("%.2f%%",freq), gene, cDNA, aa_fmt
            }' \
        | ( echo -e "Variant Type\tVariant Position\tCoverage\tLength\tFrequency\tGene\tcDNA change\tAA change"
            sort -t $'\t' -k2,2V ) > "$SIMPLE_TSV"
    else
        bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\n' "$CURRENT_SRC" \
        | awk -v OFS="\t" "$VARSCAN_AWK_FIELDS"'
            {
                pos=$1; ref=$2; alt=$3; dp=$4; freq_raw=$5;
                freq = get_freq(freq_raw)
                if(length(ref)==1&&length(alt)==1){type="SNP";len=1}
                else if(length(ref)<length(alt)){type="INS";len=length(alt)-length(ref)}
                else if(length(ref)>length(alt)){type="DEL";len=length(ref)-length(alt)}
                else{type="VAR";len=length(alt)}
                print type, tolower(ref)pos tolower(alt), dp, len, sprintf("%.2f%%",freq), "Intergenic", ref pos alt, ""
            }' \
        | ( echo -e "Variant Type\tVariant Position\tCoverage\tLength\tFrequency\tGene\tcDNA change\tAA change"
            sort -t $'\t' -k2,2V ) > "$SIMPLE_TSV"
    fi

    # --- EXTENDED TABLE: volledige CSQ annotatie ---
    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        EXTENDED_TSV="$OUT/tables/${SAMPLE}.variants.extended.tsv"
        {
            bcftools query -i 'INFO/CSQ!=""' \
                -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t%INFO/CSQ\t1\n' "$CURRENT_SRC"
            bcftools query -i 'INFO/CSQ=""' \
                -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
        } \
        | awk -v OFS="\t" "$VARSCAN_AWK_FIELDS"'
            {
                chrom=$1; pos=$2; ref=$3; alt=$4; dp=$5; freq_raw=$6; csq_str=$7; annotated=$8;
                freq = get_freq(freq_raw)
                if (annotated == "1") {
                    split(csq_str, f, "|");
                    consequence=f[1]; gene=f[2]; transcript=f[3]; biotype=f[4];
                    aa_raw=f[6]; split(f[7], dna_parts, ","); dna_change=dna_parts[1];
                    aa_fmt = (consequence ~ /synonymous/) ? "Silent mutation" : format_aa(aa_raw, consequence)
                    printf "%s\t%s\t%s\t%s\t%.2f%%\t%s\t%s\t%s\t%s\t%s\n", \
                        chrom,pos,ref,alt,freq,consequence,gene,transcript,dna_change,aa_fmt
                } else {
                    printf "%s\t%s\t%s\t%s\t%.2f%%\tIntergenic\tIntergenic\t\t%s\t\n", \
                        chrom,pos,ref,alt,freq,ref pos alt
                }
            }' \
        | ( echo -e "CHROM\tPOS\tREF\tALT\tFrequency\tConsequence\tGene\tTranscript\tcDNA_change\tAA_change"
            sort -t $'\t' -k2,2n ) > "$EXTENDED_TSV"
        say " - Extended table: $(tail -n +2 "$EXTENDED_TSV" | wc -l) varianten"
    fi

    # --- MINORITY VARIANTS TABLE: >= MINORITY_AF ---
    MINORITY_PCT=$(echo "$MINORITY_AF * 100" | bc)
    MINORITY_TSV="$OUT/tables/${SAMPLE}.minorities.tsv"
    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        {
            bcftools query -i 'INFO/CSQ!=""' \
                -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t%INFO/CSQ\t1\n' "$CURRENT_SRC"
            bcftools query -i 'INFO/CSQ=""' \
                -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
        }
    else
        bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
    fi \
    | awk -v OFS="\t" -v min_af="$MINORITY_AF" -v G_IDX="$CSQ_FIELD_GENE" -v AA_IDX="$CSQ_FIELD_AA" \
          "$VARSCAN_AWK_FIELDS"'
        {
            pos=$1; ref=$2; alt=$3; dp=$4; freq_raw=$5; csq_str=$6; annotated=$7;
            freq = get_freq(freq_raw)
            if (freq < min_af*100) next
            gene=""; aa_raw=""; consequence="";
            if (annotated == "1") {
                split(csq_str, csq_vals, "|");
                consequence=csq_vals[1]; gene=csq_vals[G_IDX]; aa_raw=csq_vals[AA_IDX];
            } else { gene = "Intergenic" }
            if(length(ref)==1&&length(alt)==1) type="SNP";
            else if(length(ref)<length(alt))   type="INS";
            else if(length(ref)>length(alt))   type="DEL";
            else                                type="VAR";
            aa_display = (consequence ~ /synonymous/) ? "Silent mutation" : format_aa(aa_raw, consequence)
            print pos, type, ref, alt, sprintf("%.2f%%",freq), dp, gene, aa_display
        }' \
    | ( echo -e "Position\tType\tRef\tAlt\tFrequency\tCoverage\tGene\tAA change"
        sort -t $'\t' -k1,1n ) > "$MINORITY_TSV"

    N_MIN=$(tail -n +2 "$MINORITY_TSV" | wc -l)
    say " - Minorities (>= ${MINORITY_PCT}%): $N_MIN varianten"

    # --- CoA TABEL ---
    COA_TSV="$OUT/tables/${SAMPLE}.coa.tsv"
    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        {
            bcftools query -i 'INFO/CSQ!=""' \
                -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t%INFO/CSQ\t1\n' "$CURRENT_SRC"
            bcftools query -i 'INFO/CSQ=""' \
                -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
        }
    else
        bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
    fi \
    | awk -v OFS="\t" -v min_af="$MINORITY_AF" \
          -v G_IDX="$CSQ_FIELD_GENE" -v AA_IDX="$CSQ_FIELD_AA" \
          -v lookup_file="$GENE_LOOKUP_FILE" \
          "$VARSCAN_AWK_FIELDS"'
        BEGIN {
            while ((getline line < lookup_file) > 0) {
                split(line, kv, "\t")
                if (kv[1] != "") gene_name[kv[1]] = kv[2]
            }
            close(lookup_file)
        }
        {
            pos=$1; ref=$2; alt=$3; dp=$4; freq_raw=$5; csq_str=$6; annotated=$7;
            freq = get_freq(freq_raw)
            if (freq < min_af*100) next
            gene_id=""; aa_raw=""; consequence="";
            if (annotated == "1") {
                split(csq_str, csq_vals, "|");
                gene_id=csq_vals[G_IDX]; aa_raw=csq_vals[AA_IDX]; consequence=csq_vals[1];
            }
            # Gene display naam uit lookup
            if (gene_id == "" || gene_id == "Intergenic") {
                gene_display = "Intergenic"
            } else if (gene_id in gene_name) {
                gene_display = gene_name[gene_id]
            } else {
                gene_display = gene_id; sub(/^gene-/, "", gene_display)
                gene_display = gene_display " protein"
            }
            # AA mutatie display
            if (gene_id == "" || gene_id == "Intergenic") {
                aa_display = "Untranslated"
            } else if (consequence ~ /synonymous/) {
                aa_display = "Silent mutation"
            } else {
                aa_display = format_aa(aa_raw, consequence)
            }
            if(length(ref)==1&&length(alt)==1) type="SNP";
            else if(length(ref)<length(alt))   type="INS";
            else if(length(ref)>length(alt))   type="DEL";
            else                                type="VAR";
            vlabel = toupper(ref) pos toupper(alt)
            len = (type=="SNP") ? 1 : (type=="INS") ? length(alt)-length(ref) : \
                  (type=="DEL") ? length(ref)-length(alt) : length(alt)
            print type, vlabel, dp, len, sprintf("%.2f%%",freq), gene_display, aa_display
        }' \
    | ( echo -e "Variant Type\tVariant Position and Identified Alternative Base\tCoverage\tLength of Variant\tFrequency of Variant\tGene (Region)\tAmino Acid Mutation"
        sort -t $'\t' -k2,2V ) > "$COA_TSV"

    N_COA=$(tail -n +2 "$COA_TSV" | wc -l)
    say " - CoA tabel: $N_COA varianten >= ${MINORITY_PCT}%"

    # --- QC DISTRIBUTIE TABEL: altijd op basis van 2% run ---
    QC_TSV="$OUT/tables/${SAMPLE}.varscan_qc.tsv"
    QC_RAW_VCF="$OUT/temp/${SAMPLE}.varscan_qc_raw.vcf"
    samtools mpileup -Q 5 -B -d "$MAXDEPTH" --min-MQ 20 \
        -f "$REF" "$BAM_FILE" 2>> "$LOG" \
        | varscan mpileup2snp \
            --min-coverage "$MINREADS" \
            --min-var-freq 0.02 \
            --min-avg-qual 5 \
            --strand-filter 0 \
            --p-value "$VARSCAN_PVAL" \
            --output-vcf 1 \
        > "$QC_RAW_VCF" 2>> "$LOG"

    bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\n' "$QC_RAW_VCF" \
    | awk "$VARSCAN_AWK_FIELDS"'
        {
            freq = get_freq($5)
            if      (freq <  3) bin="02-03%"
            else if (freq <  5) bin="03-05%"
            else if (freq < 10) bin="05-10%"
            else if (freq < 20) bin="10-20%"
            else if (freq < 50) bin="20-50%"
            else                bin=">50%  "
            count[bin]++; total++
        }
        END {
            print "AF bin\tAantal varianten\tPercentage van totaal"
            bins[1]="02-03%"; bins[2]="03-05%"; bins[3]="05-10%"
            bins[4]="10-20%"; bins[5]="20-50%"; bins[6]=">50%  "
            for(i=1;i<=6;i++){
                b=bins[i]; n=(b in count)?count[b]:0
                printf "%s\t%d\t%.1f%%\n", b, n, (total>0?n/total*100:0)
            }
            printf "Totaal\t%d\t100.0%%\n", total
        }' > "$QC_TSV"
    rm -f "$QC_RAW_VCF"

    SAMPLE_ELAPSED=$(elapsed $SAMPLE_START)
    say " - QC tabel geschreven (totaal $(tail -n1 "$QC_TSV" | cut -f2) varianten bij 2% drempel)"
    say "DONE: [$SAMPLE] in ${SAMPLE_ELAPSED}s"

    # Log sample samenvatting
    {
        echo ""
        echo "--- Resultaten $SAMPLE ---"
        echo "Runtime    : ${SAMPLE_ELAPSED}s"
        echo "Gemapt     : $N_MAPPED / $N_READS reads"
        echo "Varianten  : $N_VARS (>= $MINAF)"
        [[ "$VCF_CSQ_OK" -eq 1 ]] && echo "CSQ annot  : $N_CSQ / $N_VARS"
        echo "CoA table  : $N_COA varianten (>= $MINORITY_AF)"
        echo "Minorities : $N_MIN"
    } >> "$MAIN_LOG"

done

TOTAL_ELAPSED=$(elapsed $PIPELINE_START)
say "================================================================"
say "ALL DONE in ${TOTAL_ELAPSED}s"
say "Resultaten : $OUT/"
say "Pipeline log: $MAIN_LOG"
say "IGV        : $OUT/alignments/*.sorted.bam + $REF"
say "================================================================"

echo "" >> "$MAIN_LOG"
echo "Pipeline voltooid in ${TOTAL_ELAPSED}s op $(date)" >> "$MAIN_LOG"
