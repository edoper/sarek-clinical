# Sarek Clinical Pipeline — Plain-Language Guide

A pipeline that takes raw DNA sequencing files (**FASTQ**) from a whole genome and
finds the genetic variants in them, using **four independent variant callers** and
keeping a variant when **at least two of them agree** (with DeepVariant given priority).
The heavy computing runs on **Google Cloud**, not on your laptop.

> This guide assumes almost no cloud experience. Read it top to bottom once.

---

## 1. The mental model (read this first)

Think of it like ordering food delivery:

| Real life | This pipeline |
|---|---|
| You, on your phone | **Your laptop** — only gives orders, does no cooking |
| The restaurant kitchen | **Google Cloud** — rents powerful computers that do the actual work |
| The fridge where ingredients live | **The bucket** (`gs://intergenica-sarek-clinical`) — cloud storage for your files |
| The recipe | **Sarek** — the standard, published analysis workflow |
| The waiter taking your order to the kitchen | **Nextflow** — the program on your laptop that sends jobs to the cloud |

Key idea: **your laptop stays light.** You type one command; Google spins up big
computers, runs the analysis, writes the results into the bucket, and shuts the
computers down automatically. You only pay while they run.

**Three names you'll see a lot:**
- **Bucket** = a cloud folder. Paths start with `gs://`. Ours is `gs://intergenica-sarek-clinical`.
- **Google Batch** = the service that rents computers on demand and turns them off when done.
- **Spot** = cheap "leftover" computers (~70–90% off). They can occasionally be taken back
  mid-job; the pipeline just automatically retries, so you save money safely.

---

## 2. What's already set up (you don't need to redo this)

- **Google project:** `intergenica`  **Billing:** your "Computacion-nube" account
- **Bucket:** `gs://intergenica-sarek-clinical` with three folders:
  - `fastq/` → where you put input files
  - `work/` → scratch space the pipeline uses while running (delete it afterwards)
  - `results/` → your output (variants + quality reports)
- **Your laptop (WSL):** Java + Nextflow installed; settings file `gcb.config` ready.

---

## 3. How to run it, step by step

### Step 0 — open the toolbox (every time you start a terminal)
```bash
source ~/sarek-clinical/env.sh
```
*Why:* loads Java + Nextflow so the commands below work.

### Step 1 — log in to Google (only the very first time)
```bash
gcloud auth application-default login
gcloud auth application-default set-quota-project intergenica
```
*Why:* gives Nextflow permission to rent cloud computers on your behalf. A browser opens; pick your Google account.

### Step 2 — upload your FASTQ files to the bucket
```bash
gcloud storage cp *_R1.fastq.gz *_R2.fastq.gz gs://intergenica-sarek-clinical/fastq/
```
*Why:* the cloud computers read from the bucket, not from your laptop.

### Step 3 — make a samplesheet (a small table listing your files)
Copy `samplesheet.example.csv`, rename it `samplesheet.csv`, and edit the file paths.
Columns: `patient,sample,lane,fastq_1,fastq_2`. One row per pair of FASTQ files.

### Step 4 — run the pipeline
```bash
cd ~/sarek-clinical
nextflow run nf-core/sarek -r 3.8.1 -profile docker -c gcb.config \
  --input  samplesheet.csv \
  --outdir gs://intergenica-sarek-clinical/results/run01 \
  --genome GATK.GRCh38 \
  --tools  deepvariant,strelka,freebayes,haplotypecaller
```
*What happens:* Nextflow prints a live list of steps. Each step runs on its own cloud
computer. A whole genome takes several hours. **You can close the laptop lid? No —**
keep the terminal open (or use `-bg` to run in the background). If it stops, just add
`-resume` to the same command and it continues where it left off.

### Step 5 — get your results
Results land in `gs://intergenica-sarek-clinical/results/run01/`. Download the variant
files (VCF) and the quality report (MultiQC) when ready:
```bash
gcloud storage cp -r gs://intergenica-sarek-clinical/results/run01 ./run01-results
```

### Step 6 — clean up to stop paying for storage
```bash
gcloud storage rm -r gs://intergenica-sarek-clinical/work
```
*Why:* `work/` is large scratch space. Deleting it after you have results saves storage cost.
Keep `results/`.

---

## 4. The variant-calling logic (why four callers)

No single variant caller is perfect. We run four and combine them:

- **DeepVariant** — Google's AI caller, very accurate → **given priority**
- **Strelka2** — fast, accurate
- **HaplotypeCaller** — the long-standing standard (GATK)
- **FreeBayes** — a different method, for breadth

**Rule:** keep a variant if **≥2 callers** find it, **OR** if **DeepVariant** finds it
(so we don't throw away DeepVariant's high-quality calls). This balances accuracy
(fewer false positives) with not missing real variants. *(The script that applies this
is `consensus.sh` — see Section 5.)*

---

## 5. Combining the four callers into one list — `consensus.sh`

Sarek gives you **one VCF per caller** (a DeepVariant file, a Strelka file, etc.). For a
clinical report you want **one** list of variants, where each variant is tagged with *how
much the callers agreed*. That is what `consensus.sh` does.

### What it does, in plain words
It builds a **union** list with two parts:

1. **Backbone — every DeepVariant call**, kept with DeepVariant's genotype details:
   zygosity (`GT`), quality (`GQ`), depth (`DP`), allelic depths (`AD`), allele fraction
   (`VAF`). *(Sarek's own built-in consensus throws these details away — that's why we use
   our own script.)*
2. **Rescue — variants DeepVariant missed but ≥2 of the other callers agreed on.** These are
   added back so we don't lose real variants. Since DeepVariant has no genotype for them, the
   genotype is **borrowed from Strelka2**, or from **HaplotypeCaller** if Strelka2 didn't call
   it (FreeBayes counts toward agreement but is never used as the genotype source). Rescued
   variants carry `GT/GQ/DP/AD` (no `VAF` — those callers don't report it; `AD` lets you
   compute allele fraction).

Every variant gets four tags:
- `NCALLERS` — how many of the four callers found it (1–4)
- `CALLERS`  — their names, e.g. `deepvariant,strelka`
- `CONF`     — confidence from agreement: **HIGH** (≥3 callers), **MEDIUM** (2), **LOW** (1)
- `GT_SOURCE` — which caller the genotype came from (`deepvariant`, `strelka`, or `haplotypecaller`)

Nothing concordant is silently lost, and **nothing is tiered away here** — you decide how
strict to be *afterwards* (see below). This keeps the decision visible and auditable.

> Before comparing, the script lines up the callers fairly: keeps only `PASS` variants,
> splits "two-variants-in-one-line" records apart, and shifts insertions/deletions to a
> standard position (so the same indel written differently by two callers still matches).

### What you need
- `bcftools`, `bgzip`, `tabix`, `samtools` available (e.g. `conda install -c bioconda bcftools htslib samtools`).
- The **reference genome FASTA** your data was aligned to, with its `.fai` index next to it
  (GATK.GRCh38 for this project). If the `.fai` is missing: `samtools faidx your_reference.fasta`.
- The four per-caller VCFs (`.vcf.gz`) for **one sample**.

### How to run it
```bash
source ~/sarek-clinical/env.sh        # puts tools on PATH (if you installed them there)

~/sarek-clinical/consensus.sh \
  -r /path/to/GATK.GRCh38.fasta \
  -d  results/variant_calling/deepvariant/SAMPLE/SAMPLE.deepvariant.vcf.gz \
  -o  results/consensus/SAMPLE \
  -c strelka=results/variant_calling/strelka/SAMPLE/SAMPLE.strelka.variants.vcf.gz \
  -c freebayes=results/variant_calling/freebayes/SAMPLE/SAMPLE.freebayes.vcf.gz \
  -c haplotypecaller=results/variant_calling/haplotypecaller/SAMPLE/SAMPLE.haplotypecaller.vcf.gz
```
- `-r` reference FASTA (with `.fai`)   · `-d` the DeepVariant VCF (the backbone)
- `-o` output prefix   · `-c name=path` for each of the **other** callers (repeat as needed)
- Run `~/sarek-clinical/consensus.sh -h` for all options (e.g. `-f` to change which FILTER
  values are kept; default `PASS,.`).

You can run this **on your laptop** — it only reads the small variant files, not the raw
sequencing data, so it does not need the cloud.

### What you get
- `SAMPLE.consensus.vcf.gz` (+ `.tbi`) — all DeepVariant calls **plus** the rescued ≥2-caller
  variants, genotypes intact, each tagged with `CALLERS` / `NCALLERS` / `CONF` / `GT_SOURCE`.
- `SAMPLE.consensus.log` — a record of exactly what was run (for clinical provenance).

### Choosing how strict to be (the actual filter)
The file is exactly the **"≥2 callers OR DeepVariant"** set from Section 4 — most sensitive,
nothing real thrown away. To tighten it later, filter on the tags:
```bash
# Sensitive (default): use the file as-is.

# Stricter — keep only variants 2+ callers agreed on (drops DeepVariant-only calls):
bcftools view -i 'NCALLERS>=2' SAMPLE.consensus.vcf.gz -Oz -o SAMPLE.concordant.vcf.gz

# Highest-confidence tier only:
bcftools view -i 'CONF="HIGH"' SAMPLE.consensus.vcf.gz

# See where each genotype came from:
bcftools query -f '%CHROM\t%POS\t%INFO/GT_SOURCE\t%INFO/CALLERS\n' SAMPLE.consensus.vcf.gz
```

> **Mixed genotype sources:** rescued variants carry a Strelka2/HaplotypeCaller genotype (and
> no `VAF`). Use `GT_SOURCE` if your downstream filter treats DeepVariant genotypes differently.
This is where the project's downstream **candidate-filtering** step takes over.

> **Clinical note:** `CONF` measures *agreement between callers*, not absolute correctness.
> A `LOW` (DeepVariant-only) variant is "unconfirmed by the others," **not** "wrong" —
> DeepVariant alone is still high quality.

---

## 6. Money & safety

- **Cost:** roughly **$50–100 for 4 genomes** (estimate), thanks to Spot computers.
- **You only pay while jobs run** + a small amount for files sitting in the bucket.
- The biggest waste is forgetting to delete the `work/` folder — do Step 6.
- **Clinical note:** before using results for patient care, the pipeline must be
  formally validated (run a known reference genome, "GIAB", and check the accuracy numbers).

---

## 7. Mini-glossary

| Term | Plain meaning |
|---|---|
| FASTQ | raw sequencing reads (the input) |
| VCF | a list of genetic variants (the output) |
| Bucket / `gs://` | a folder in Google Cloud Storage |
| Nextflow | the program that sends jobs to the cloud |
| Sarek | the published analysis recipe we run |
| Google Batch | the service that rents/returns cloud computers |
| Spot VM | a cheap, occasionally-interrupted cloud computer (auto-retried) |
| `-resume` | "continue where it stopped" — safe to reuse |
| GIAB | a reference genome with known answers, used to test accuracy |
| consensus | combining the four callers' results into one tagged list |
| `NCALLERS` / `CALLERS` | how many callers (and which) found a given variant |
| `CONF` | confidence from agreement: HIGH (≥3 callers) / MEDIUM (2) / LOW (1) |
| `GT_SOURCE` | which caller the genotype came from (`deepvariant` / `strelka` / `haplotypecaller`) |
| genotype (`GT`,`DP`,`AD`,`GQ`,`VAF`) | per-variant details: zygosity, depth, allele depths, quality, allele fraction |

---

## Files in this repo
- `README.md` — this guide
- `gcb.config` — the cloud settings for real runs (project, bucket, Spot, 4 callers)
- `gcb-smoke.config` — cloud settings for the tiny test run (no `params`, lets the `test` profile drive)
- `env.sh` — loads Java + Nextflow into your terminal (and sets `NXF_SYNTAX_PARSER=v1`, required for sarek 3.8.1 on Nextflow 26.x)
- `samplesheet.example.csv` — template for listing your input files
- `consensus.sh` — union consensus: all DeepVariant calls + variants ≥2 other callers agree on (genotype borrowed from Strelka2/HaplotypeCaller), tagged with `CALLERS`/`NCALLERS`/`CONF`/`GT_SOURCE` (Section 5)
