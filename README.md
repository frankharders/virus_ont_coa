# virus_ont_coa

**Nanopore ONT virus minority variant calling pipeline voor Certificate of Analysis (CoA) rapportage**

> *Your science, our sequencing, endless possibilities*

---

> ⚠️ **BELANGRIJK: Single-segment virussen only**
> 
> Deze pipeline is uitsluitend geschikt voor virussen met **één genomisch segment**.
> Voorbeelden: SARS-CoV-2, MERS-CoV, overige coronavirussen.
> 
> Voor **multi-segment virussen** (BTV, EHDV, NDV, SFTS, Avian Influenza)
> gebruik de aparte pipeline: [`virus_ont_coa_multiseg`](https://github.com/frankharders/virus_ont_coa_multiseg)

---

## Overzicht

`virus_ont_coa` is een complete bioinformatica pipeline voor het detecteren en rapporteren van minority variants in virus isolaten gesequenced met Oxford Nanopore Technology (ONT). De pipeline is specifiek ontwikkeld voor de productie van Certificate of Analysis (CoA) documenten bij WBVR en genereert zowel gestandaardiseerde tabellen als interactieve visualisaties.

### Wat doet de pipeline?

1. **Alignment** — Nanopore FASTQ reads worden gealigneerd tegen een virusreferentie met `minimap2`
2. **Variant calling** — Minority variants worden gedetecteerd met `VarScan2` (`mpileup2snp`), geoptimaliseerd voor ONT error profielen
3. **Annotatie** — Varianten worden geannoteerd met gene- en aminozuurinformatie via `bcftools csq`
4. **Tabellen** — Meerdere TSV tabellen per sample, waaronder een CoA-klare tabel
5. **Visualisatie** — Interactief HTML rapport met genome-brede coverage plot, per-gen coverage en IGV-stijl SNP snippets

### Waarom VarScan2 voor ONT?

ONT R10.4.1 SUP basecalling heeft een systematische error rate van ~0.5–1% per base. Bij hoge coverage (>500x, typisch voor virus isolaten) zijn zelfs 2% minorities statistisch significant (p<E-20), waardoor p-waarde filtering weinig effect heeft. VarScan2 `mpileup2snp` werkt zonder genotype aannames (geen ploidy bias) en rapporteert de werkelijke allele frequentie — essentieel voor quasi-species virus populaties. De aanbevolen afkapwaarde voor ONT is **5% AF**, onderbouwd door de meegeleverde QC distributie tabel per sample.

---

## Repo structuur

```
virus_ont_coa/
├── run_pipeline.sh                     ← Wrapper: roept beide scripts aan
├── scripts/
│   ├── virus_ont_variants_coa.sh       ← Variant calling pipeline
│   └── virus_ont_variants_coa_viz.sh   ← Visualisatie
├── envs/
│   ├── coa_env.yml                     ← Conda environment variant calling
│   └── coa_viz_env.yml                 ← Conda environment visualisatie
└── README.md
```

> **Let op:** De `references/` map (FASTA + GFF3) en `results/` map worden **niet** meegeleverd in de repo. Voeg ze toe aan `.gitignore` of gebruik `git-lfs` voor grote referentiebestanden.

---

## Vereisten

### Hardware
- **Host**: `lelycompute-02.wur.nl` (HPC node, 754 GB RAM)  
  De scripts bevatten een hostname check en weigeren te starten op andere machines.  
  Reden: de conda environments zijn alleen beschikbaar op lelycompute-02.

### Software
Alle tools zijn beschikbaar via de meegeleverde conda environments.

| Tool | Versie | Gebruik |
|------|--------|---------|
| minimap2 | ≥2.24 | ONT alignment (`-ax map-ont`) |
| samtools | ≥1.17 | BAM verwerking |
| bcftools | ≥1.17 | VCF filtering, CSQ annotatie |
| VarScan2 | ≥2.4.6 | Minority variant calling |
| tabix/bgzip | ≥1.17 | VCF indexering |
| python3 | ≥3.11 | MNV pointer resolutie, visualisatie |
| pysam | ≥0.21 | BAM coverage extractie |
| plotly | ≥5.x | Interactieve HTML visualisaties |
| pandas | ≥2.x | Data verwerking |
| biopython | ≥1.81 | FASTA sequence handling |

---

## Installatie

### 1. Repository clonen

```bash
git clone https://github.com/frankharders/virus_ont_coa.git
cd virus_ont_coa
```

### 2. Conda environments aanmaken

```bash
# Variant calling environment
mamba env create -f envs/coa_env.yml

# Visualisatie environment
mamba env create -f envs/coa_viz_env.yml
```

> `mamba` wordt aanbevolen boven `conda` voor snellere dependency resolutie. Installeer via: `conda install -n base -c conda-forge mamba`

### 3. Referentiebestanden voorbereiden

```bash
mkdir references/
# Kopieer of download je virus FASTA en GFF3:
# references/<virus>.fa
# references/<virus>.gff3
```

**Belangrijk:** Windows regeleindes veroorzaken problemen in Clair3/bcftools. Verwijder ze preventief:
```bash
sed -i 's/\r//' references/<virus>.fa
sed -i 's/\r//' references/<virus>.gff3
samtools faidx references/<virus>.fa
```

---

## Gebruik

### Snelstart — complete pipeline

```bash
conda activate coa
bash run_pipeline.sh \
    -r references/KC776174.fa \
    -g references/KC776174.gff3 \
    -i 01_basecalled \
    -o results \
    -t 16
```

### Alleen variant calling (geen visualisatie)

```bash
conda activate coa
bash run_pipeline.sh \
    -r references/KC776174.fa \
    -g references/KC776174.gff3 \
    -i 01_basecalled \
    -o results \
    -t 16 \
    -V
```

### Alleen visualisatie (op bestaande resultaten)

```bash
conda activate coa_viz
bash scripts/virus_ont_variants_coa_viz.sh \
    -o results \
    -r references/KC776174.fa \
    -g references/KC776174.gff3
```

### Alle CLI opties

```
run_pipeline.sh opties:

  VERPLICHT
    -r FILE   Referentie FASTA

  VARIANT CALLING
    -g FILE   GFF3 annotatie (voor gene/AA info en gen grenzen in plots)
    -i DIR    Map met Nanopore FASTQ reads          (default: 01_basecalled)
    -o DIR    Output map                            (default: results)
    -t INT    Threads voor minimap2/samtools         (default: 8)
    -d INT    Minimale read depth                   (default: 5)
    -f FLOAT  Minimale allele frequentie            (default: 0.05)
    -m FLOAT  Minimale AF voor CoA/minorities tabel (default: 0.05)
    -q FLOAT  VarScan2 p-waarde drempel             (default: 0.001)
    -D INT    Max read depth per positie mpileup    (default: 50000)

  VISUALISATIE
    -c INT    Context bp links/rechts van SNP       (default: 10)
    -V        Sla visualisatie over
```

---

## Output

Alle output wordt weggeschreven naar de opgegeven output map (`-o`).

### Tabellen (`results/tables/`)

Per sample worden de volgende TSV bestanden gegenereerd:

| Bestand | Inhoud |
|---------|--------|
| `*.variants.tsv` | Alle varianten ≥ MINAF, gesorteerd op positie |
| `*.variants.extended.tsv` | Uitgebreide tabel met consequence, gen, transcript, cDNA en AA change |
| `*.minorities.tsv` | Minority variants ≥ MINORITY_AF met gen/AA annotatie |
| `*.coa.tsv` | **CoA tabel**: leesbare gen namen, gestandaardiseerde AA notatie |
| `*.varscan_qc.tsv` | AF distributie 2–50% voor assay ruis karakterisatie |

#### CoA tabel kolommen

| Kolom | Beschrijving |
|-------|-------------|
| Variant Type | SNP / INS / DEL / VAR |
| Variant Position and Identified Alternative Base | bijv. `C24076T` |
| Coverage | Aantal reads op de positie |
| Length of Variant | Lengte in bp |
| Frequency of Variant | Allele frequentie in % |
| Gene (Region) | Gen naam uit GFF3, bijv. `S protein` |
| Amino Acid Mutation | bijv. `G15D`, `Silent mutation`, `Untranslated`, `Q345X` |

#### AA notatie

- Missense: `G15D` (ref_AA + positie + alt_AA)
- Synoniem: `Silent mutation`
- Stop codon introductie: `Q345X` (`X` = stop, HGVS-conventie voor rapportage)
- Intergenic: `Untranslated`

#### QC distributie tabel

De `*.varscan_qc.tsv` tabel toont altijd de volledige AF distributie vanaf 2% (ongeacht de `-f` instelling), zodat de assay ruis karakterisatie zichtbaar blijft:

```
AF bin    Aantal varianten    Percentage van totaal
02-03%    203                 54.9%
03-05%    105                 28.4%
05-10%    55                  14.9%
10-20%    4                   1.1%
20-50%    3                   0.8%
>50%      0                   0.0%
Totaal    370                 100.0%
```

Bij ONT R10.4.1 SUP is ~83% van de 2–5% variants systematische sequencing ruis. Dit onderbouwt de keuze voor 5% als afkapwaarde.

### Alignments (`results/alignments/`)

- `*.sorted.bam` — gesorteerde en geïndexeerde BAM bestanden
- Geschikt voor directe inspectie in IGV

### VCF bestanden (`results/variants/`)

- `*.raw.vcf.gz` — gefilterde VarScan2 output (BGZF + tabix index)
- `*.csq.vcf.gz` — bcftools csq geannoteerde VCF (inclusief MNV pointer resolutie)

### Visualisaties (`results/viz/`)

Per sample één interactief HTML rapport met:

- **Genome-brede coverage plot** — volledig genoom met gen annotaties als gekleurde achtergrondvlakken en SNP posities als rode stippellijnen
- **Per-gen coverage plots** — ingezoomd op elk gen met SNP labels (AA change)
- **IGV-stijl SNP snippets** — voor elke variant ≥ MINORITY_AF:
  - ±10 bp context (instelbaar via `-c`)
  - Gestapelde base frequentie balken (A=groen, T=rood, C=blauw, G=oranje)
  - Referentie base op x-as
  - SNP positie rood gemarkeerd

Het HTML rapport is interactief (zoom, hover, pan) en kan via de browser worden afgedrukt als PDF (Chrome: Ctrl+P → Save as PDF).

---

## Technische details

### Variant calling instellingen

VarScan2 wordt aangeroepen via `samtools mpileup | varscan mpileup2snp` met de volgende ONT-specifieke instellingen:

| Parameter | Waarde | Reden |
|-----------|--------|-------|
| `-B` | aan | Geen BAQ correctie — voorkomt vals-negatieve calls nabij indels bij ONT reads |
| `-Q 5` | 5 | Lage base quality filter — ONT Q-scores zijn niet Illumina-equivalent |
| `--min-MQ 20` | 20 | Mapping quality filter — verwijdert multi-mapping reads |
| `--strand-filter 0` | 0 | Geen strand bias filter — niet van toepassing bij ONT (asymmetrische strand distributie) |
| `-d 50000` | 50000 | Max depth — voorkomt afkapping bij hoge coverage virus isolaten |

### MNV pointer resolutie

`bcftools csq` annoteert multi-nucleotide variants (MNV, bijv. codon changes over meerdere posities) alleen op de eerste positie. Latere posities krijgen een `CSQ=@<refpos>` pointer. Een ingebouwde Python resolver kopieert de annotatie van de referentiepositie naar alle pointer-posities, zodat elke variant volledig geannoteerd is in de output tabellen.

### GFF3 vereisten

De pipeline accepteert zowel standaard NCBI GFF3 als minimale CDS-only GFF3 bestanden. Als er geen `mRNA`/`transcript` features aanwezig zijn, worden deze automatisch gesynthetiseerd uit de CDS features. Seqid normalisatie (aanpassen aan FASTA header) gebeurt automatisch.

---

## Bekende beperkingen

- **Indels**: VarScan2 `mpileup2snp` roept alleen SNPs aan. Indels worden niet gerapporteerd (ONT indel calling is onbetrouwbaar zonder speciaal getrainde modellen).
- **Single-segment only**: De pipeline is uitsluitend ontworpen voor virussen met één genomisch segment. Multi-segment virussen (BTV, EHDV, NDV, SFTS, Avian Influenza) worden **niet ondersteund** en geven onjuiste resultaten. Gebruik hiervoor [`virus_ont_coa_multiseg`](https://github.com/frankharders/virus_ont_coa_multiseg).
- **Ribosomal slippage**: ORF1ab-achtige genen met ribosomal slippage genereren een bcftools csq waarschuwing. De annotatie is correct maar de overlappende CDS regio wordt tweemaal gerapporteerd.

---

## Citatie

Als je deze pipeline gebruikt in een publicatie of rapport, citeer dan:

```
Harders, F. (2026). virus_ont_coa: Nanopore ONT virus minority variant calling
pipeline voor CoA rapportage. Veterinary Integrative Data Analytics Department,
Wageningen Bioveterinary Research, Lelystad.
https://github.com/frankharders/virus_ont_coa
```

---

## Contact

**Frank Harders**  
Veterinary Integrative Data Analytics Department  
Wageningen Bioveterinary Research  
PO Box 65, 8200 AB Lelystad, Nederland  
Houtribweg 39, 8221 RA Lelystad (ASG, Room 215.125)  
Tel: +31-(0)320-238273  
GitHub: [@frankharders](https://github.com/frankharders)

---

*Wageningen Bioveterinary Research — Your science, our sequencing, endless possibilities*
