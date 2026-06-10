#!/usr/bin/env bash
# =============================================================================
# virus_ont_variants_coa_viz.sh
# Visualisatie script voor virus ONT CoA variant calling resultaten
# Genereert interactief HTML rapport + PDF per sample
# =============================================================================
#
# BESCHRIJVING
#   Leest de output van virus_ont_variants_coa.sh en genereert per sample:
#   - Interactief HTML rapport (Plotly)
#   - PDF rapport voor CoA rapportage
#   Inhoud:
#   - Genome-brede coverage plot met SNP posities en gen namen
#   - Per-gen coverage plots (ingezoomd)
#   - IGV-stijl SNP snippets (10bp context, base frequenties)
#
# GEBRUIK
#   bash virus_ont_variants_coa_viz.sh -o <results_map> -r <referentie.fa> [opties]
#
# VERPLICHT
#   -o DIR    Output map van virus_ont_variants_coa.sh (bevat alignments/ en tables/)
#   -r FILE   Referentie FASTA (zelfde als gebruikt voor variant calling)
#
# OPTIONEEL
#   -g FILE   GFF3 annotatie (voor gen grenzen in coverage plot)
#   -m FLOAT  Minimale AF voor SNP snippets (default: 0.05)
#   -c INT    Context bp links/rechts van SNP in snippets (default: 10)
#   -h        Toon deze helptext
#
# OUTPUT (per sample in <results_map>/viz/)
#   <sample>.html   Interactief Plotly rapport
#   <sample>.pdf    PDF voor CoA rapportage
#
# VEREISTE ENV
#   conda activate coa_viz
#   (python=3.11, pysam, plotly, pandas, biopython, kaleido)
#
# AUTEUR  : WBVR Bioinformatics (harde004)
# VERSIE  : 1.0
# DATUM   : 2026-06-10

set -Eeuo pipefail

############################################
# Hostname check
############################################
REQUIRED_HOST="lelycompute-02"
CURRENT_HOST=$(hostname -s)
if [[ "$CURRENT_HOST" != "$REQUIRED_HOST" ]]; then
    echo "ERROR: Dit script moet draaien op ${REQUIRED_HOST}.wur.nl" >&2
    echo "       Huidige host: ${CURRENT_HOST}" >&2
    exit 1
fi

############################################
# Conda environment activatie
############################################
CONDA_ENV="coa_viz"
if [[ -z "${CONDA_DEFAULT_ENV:-}" || "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; then
    if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "${HOME}/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        source "/opt/conda/etc/profile.d/conda.sh"
    else
        echo "ERROR: conda niet gevonden." >&2; exit 1
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

OUT=""
REF=""
GFF_IN=""
MINORITY_AF=0.05
CONTEXT=10

while getopts ":o:r:g:m:c:h" opt; do
    case $opt in
        o) OUT="$OPTARG" ;;
        r) REF="$OPTARG" ;;
        g) GFF_IN="$OPTARG" ;;
        m) MINORITY_AF="$OPTARG" ;;
        c) CONTEXT="$OPTARG" ;;
        h) usage 0 ;;
        :) echo "ERROR: Optie -$OPTARG vereist een argument." >&2; usage 1 ;;
        \?) echo "ERROR: Onbekende optie -$OPTARG." >&2; usage 1 ;;
    esac
done

errors=0
[[ -z "$OUT" ]]               && { echo "ERROR: Results map verplicht (-o)." >&2; errors=1; }
[[ -n "$OUT" && ! -d "$OUT" ]] && { echo "ERROR: Results map niet gevonden: $OUT" >&2; errors=1; }
[[ -z "$REF" ]]               && { echo "ERROR: Referentie FASTA verplicht (-r)." >&2; errors=1; }
[[ -n "$REF" && ! -f "$REF" ]] && { echo "ERROR: Referentie FASTA niet gevonden: $REF" >&2; errors=1; }
[[ -n "$GFF_IN" && ! -f "$GFF_IN" ]] && { echo "ERROR: GFF3 niet gevonden: $GFF_IN" >&2; errors=1; }
[[ $errors -gt 0 ]] && usage 1

ts(){ date +"%F %T"; }
say(){ echo "[$(ts)] $*" >&2; }
die(){ echo "[$(ts)] ERROR: $*" >&2; exit 1; }

mkdir -p "$OUT/viz"

############################################
# Zoek samples
############################################
shopt -s nullglob
BAM_FILES=( "$OUT/alignments/"*.sorted.bam )
[[ ${#BAM_FILES[@]} -gt 0 ]] || die "Geen BAM bestanden gevonden in $OUT/alignments/"

say "================================================================"
say " virus_ont_variants_coa_viz.sh v1.0 - WBVR Bioinformatics"
say "================================================================"
say "Gevonden: ${#BAM_FILES[@]} sample(s)"

############################################
# Python visualisatie per sample
############################################
for BAM_FILE in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "$BAM_FILE" .sorted.bam)
    COA_TSV="$OUT/tables/${SAMPLE}.coa.tsv"
    OUT_HTML="$OUT/viz/${SAMPLE}.html"
    OUT_PDF="$OUT/viz/${SAMPLE}.pdf"

    [[ -f "$COA_TSV" ]] || { say "SKIP: $SAMPLE — geen CoA tabel gevonden"; continue; }

    say "Visualiseren: $SAMPLE"

    python3 << PYEOF
import sys
import pysam
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from Bio import SeqIO
import re

# ── Configuratie ──────────────────────────────────────────────
BAM_FILE    = "$BAM_FILE"
REF_FILE    = "$REF"
COA_TSV     = "$COA_TSV"
GFF_IN      = "$GFF_IN"
OUT_HTML    = "$OUT_HTML"
OUT_PDF     = "$OUT_PDF"
SAMPLE      = "$SAMPLE"
MINORITY_AF = float("$MINORITY_AF")
CONTEXT     = int("$CONTEXT")

# Kleurenschema IGV-stijl
BASE_COLORS = {"A": "#2ECC40", "T": "#CC2900", "C": "#0055CC", "G": "#CC7700", "N": "#AAAAAA"}

# ── Referentie inlezen ─────────────────────────────────────────
ref_seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(REF_FILE, "fasta")}
ref_id   = list(ref_seqs.keys())[0]
ref_seq  = ref_seqs[ref_id]
ref_len  = len(ref_seq)

# ── GFF3 gen grenzen inlezen ───────────────────────────────────
genes = []
if GFF_IN:
    with open(GFF_IN) as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            if parts[2] != "gene": continue
            start, end = int(parts[3]), int(parts[4])
            attrs = {k: v for k, v in (a.split("=") for a in parts[8].split(";") if "=" in a)}
            name = attrs.get("Name", attrs.get("ID", "unknown"))
            name = re.sub(r"^gene-", "", name)
            genes.append({"name": name, "start": start, "end": end})

# ── CoA tabel inlezen ─────────────────────────────────────────
coa = pd.read_csv(COA_TSV, sep="\t")
coa.columns = coa.columns.str.strip()
# Extraheer positie uit "A222G" kolom
coa["pos"] = coa["Variant Position and Identified Alternative Base"].str.extract(r"(\d+)").astype(int)
coa["freq_val"] = coa["Frequency of Variant"].str.replace("%","").astype(float)
snps = coa[coa["freq_val"] >= MINORITY_AF * 100].copy()

# ── Coverage berekenen ─────────────────────────────────────────
print(f"  Coverage berekenen voor {ref_id} ({ref_len} bp)...")
bam = pysam.AlignmentFile(BAM_FILE, "rb")

# Genome-brede coverage (gesamplede resolutie voor snelheid)
STEP = max(1, ref_len // 2000)
positions, coverages = [], []
for i in range(0, ref_len, STEP):
    cov = sum(1 for _ in bam.fetch(ref_id, i, i+1))
    positions.append(i + 1)
    coverages.append(cov)

# Exacte coverage op SNP posities
snp_covs = {}
for pos in snps["pos"].unique():
    c = sum(1 for _ in bam.fetch(ref_id, pos-1, pos))
    snp_covs[pos] = c

# Per-gen coverage (exacte resolutie)
gene_covs = {}
for g in genes:
    gpos, gcov = [], []
    step = max(1, (g["end"] - g["start"]) // 500)
    for p in range(g["start"], g["end"], step):
        c = sum(1 for _ in bam.fetch(ref_id, p-1, p))
        gpos.append(p)
        gcov.append(c)
    gene_covs[g["name"]] = (gpos, gcov)

# Pileup data voor SNP snippets
pileup_data = {}
for _, row in snps.iterrows():
    pos = int(row["pos"])
    start = max(0, pos - CONTEXT - 1)
    end   = min(ref_len, pos + CONTEXT)
    snippet = {}
    for col in bam.pileup(ref_id, start, end, truncate=True):
        p = col.reference_pos + 1
        bases = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}
        for pileup_read in col.pileups:
            if not pileup_read.is_del and not pileup_read.is_refskip:
                b = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                bases[b] = bases.get(b, 0) + 1
        snippet[p] = bases
    pileup_data[pos] = snippet

bam.close()

# ── Figuur opbouwen ───────────────────────────────────────────
n_genes    = len(genes)
n_snps     = len(snps)
# Rijen: 1 genome + n_genes gen plots + n_snps snippets
row_heights = [0.6] + [0.35] * n_genes + [0.6] * n_snps
n_rows      = 1 + n_genes + n_snps

subplot_titles = [f"{SAMPLE} — Genome coverage"] + \
                 [f"Coverage: {g['name']}" for g in genes] + \
                 [f"SNP {row['Variant Position and Identified Alternative Base']} | {row['Gene (Region)']} | {row['Amino Acid Mutation']} | {row['Frequency of Variant']}"
                  for _, row in snps.iterrows()]

# Spacing dynamisch berekenen op basis van aantal rijen
v_spacing = min(0.015, 0.6 / max(n_rows - 1, 1))

fig = make_subplots(
    rows=n_rows, cols=1,
    subplot_titles=subplot_titles,
    row_heights=row_heights,
    vertical_spacing=v_spacing
)

# ── Plot 1: Genome-brede coverage ─────────────────────────────
fig.add_trace(go.Scatter(
    x=positions, y=coverages,
    fill="tozeroy", fillcolor="rgba(70,130,180,0.3)",
    line=dict(color="steelblue", width=1),
    name="Coverage", showlegend=False,
    hovertemplate="Pos: %{x}<br>Coverage: %{y}<extra></extra>"
), row=1, col=1)

# Gen annotaties als achtergrond kleurvlakken
colors_genes = ["rgba(255,200,100,0.15)", "rgba(100,200,150,0.15)",
                "rgba(200,100,200,0.15)", "rgba(100,150,255,0.15)",
                "rgba(255,150,100,0.15)", "rgba(150,255,150,0.15)"]
for gi, g in enumerate(genes):
    fig.add_vrect(
        x0=g["start"], x1=g["end"],
        fillcolor=colors_genes[gi % len(colors_genes)],
        layer="below", line_width=0,
        annotation_text=g["name"],
        annotation_position="top left",
        annotation_font_size=9,
        row=1, col=1
    )

# SNP posities als verticale lijnen op genome plot
for _, row in snps.iterrows():
    pos = int(row["pos"])
    freq = row["freq_val"]
    gene = row["Gene (Region)"]
    aa   = row["Amino Acid Mutation"]
    fig.add_vline(
        x=pos, line_color="red", line_width=1, line_dash="dot",
        annotation_text=f"{pos} ({freq:.1f}%)",
        annotation_font_size=8, annotation_textangle=-90,
        row=1, col=1
    )

# ── Per-gen coverage plots ─────────────────────────────────────
for gi, g in enumerate(genes):
    row_idx = 2 + gi
    gpos, gcov = gene_covs.get(g["name"], ([], []))
    if not gpos:
        continue
    fig.add_trace(go.Scatter(
        x=gpos, y=gcov,
        fill="tozeroy", fillcolor=colors_genes[gi % len(colors_genes)].replace("0.15","0.4"),
        line=dict(color="steelblue", width=1),
        name=g["name"], showlegend=False,
        hovertemplate="Pos: %{x}<br>Coverage: %{y}<extra></extra>"
    ), row=row_idx, col=1)
    # SNPs in dit gen
    gene_snps = snps[snps["Gene (Region)"].str.contains(g["name"].replace(" protein",""), na=False)]
    for _, srow in gene_snps.iterrows():
        spos = int(srow["pos"])
        if g["start"] <= spos <= g["end"]:
            fig.add_vline(
                x=spos, line_color="red", line_width=1.5, line_dash="dot",
                annotation_text=srow["Amino Acid Mutation"],
                annotation_font_size=8, annotation_textangle=-90,
                row=row_idx, col=1
            )
    fig.update_xaxes(range=[g["start"], g["end"]], row=row_idx, col=1)

# ── SNP snippets: IGV-stijl pileup ────────────────────────────
for si, (_, srow) in enumerate(snps.iterrows()):
    row_idx = 2 + n_genes + si
    pos     = int(srow["pos"])
    snippet = pileup_data.get(pos, {})
    if not snippet: continue

    positions_snip = sorted(snippet.keys())
    ref_bases = [ref_seq[p-1] if 0 <= p-1 < ref_len else "N" for p in positions_snip]

    # Referentie sequentie als tekst annotaties (bovenste rij)
    fig.add_trace(go.Bar(
        x=positions_snip,
        y=[0] * len(positions_snip),
        text=ref_bases,
        textposition="outside",
        marker_color="rgba(0,0,0,0)",
        showlegend=False,
        hoverinfo="skip",
        name="Ref"
    ), row=row_idx, col=1)

    # Gestapelde balken per base
    all_bases = ["A", "T", "C", "G", "N"]
    for base in all_bases:
        counts = [snippet[p].get(base, 0) for p in positions_snip]
        total  = [sum(snippet[p].values()) for p in positions_snip]
        pcts   = [c/t*100 if t > 0 else 0 for c, t in zip(counts, total)]
        ref_b  = [ref_seq[p-1].upper() if 0 <= p-1 < ref_len else "N" for p in positions_snip]
        hover  = [f"Pos: {p}<br>Ref: {r}<br>{base}: {c} ({pct:.1f}%)"
                  for p, r, c, pct in zip(positions_snip, ref_b, counts, pcts)]
        fig.add_trace(go.Bar(
            x=positions_snip,
            y=counts,
            name=base,
            marker_color=BASE_COLORS[base],
            text=[base if c > 0 else "" for c in counts],
            textposition="inside",
            textfont=dict(size=9, color="white"),
            hovertext=hover,
            hoverinfo="text",
            showlegend=(si == 0),
            legendgroup=base
        ), row=row_idx, col=1)

    # Markeer de SNP positie
    fig.add_vline(
        x=pos, line_color="red", line_width=2,
        annotation_text=f"SNP {pos}",
        annotation_font_size=9,
        row=row_idx, col=1
    )
    # Referentie base annotatie
    fig.add_annotation(
        x=pos, y=0,
        text=f"REF: {ref_seq[pos-1]}",
        showarrow=False,
        font=dict(size=10, color="black"),
        bgcolor="white",
        bordercolor="red",
        borderwidth=1,
        row=row_idx, col=1
    )
    fig.update_xaxes(
        tickmode="array",
        tickvals=positions_snip,
        ticktext=[f"{p}<br>{ref_seq[p-1] if 0<=p-1<ref_len else 'N'}" for p in positions_snip],
        tickfont=dict(size=9, family="monospace"),
        row=row_idx, col=1
    )
    fig.update_yaxes(title_text="Reads", row=row_idx, col=1)

# ── Layout & export ───────────────────────────────────────────
total_height = 500 + n_genes * 350 + n_snps * 500
fig.update_layout(
    title=dict(
        text=f"<b>CoA Variant Rapport — {SAMPLE}</b>",
        font=dict(size=18)
    ),
    height=total_height,
    barmode="stack",
    template="plotly_white",
    legend=dict(
        title="Base",
        orientation="h",
        yanchor="bottom", y=1.01,
        xanchor="right", x=1
    ),
    font=dict(family="Arial", size=11),
    margin=dict(l=60, r=40, t=80, b=40)
)
fig.update_xaxes(title_text="Genomische positie (bp)", row=1, col=1)
fig.update_yaxes(title_text="Coverage (reads)", row=1, col=1)

# HTML opslaan
fig.write_html(OUT_HTML, include_plotlyjs="cdn")
print(f"  HTML: {OUT_HTML}")

# PDF opslaan
try:
    fig.write_image(OUT_PDF, format="pdf", width=1400, height=total_height)
    print(f"  PDF:  {OUT_PDF}")
except Exception as e:
    print(f"  PDF MISLUKT: {e}", file=sys.stderr)

print(f"  Klaar: {n_snps} SNP snippets, {n_genes} gen plots")
PYEOF

    say "  -> $SAMPLE klaar: $OUT_HTML"
done

say "================================================================"
say "ALL DONE. Rapporten in $OUT/viz/"
say "================================================================"