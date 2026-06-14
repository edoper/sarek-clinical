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
(fewer false positives) with not missing real variants. *(This consensus step is a
separate script — added next.)*

---

## 5. Money & safety

- **Cost:** roughly **$50–100 for 4 genomes** (estimate), thanks to Spot computers.
- **You only pay while jobs run** + a small amount for files sitting in the bucket.
- The biggest waste is forgetting to delete the `work/` folder — do Step 6.
- **Clinical note:** before using results for patient care, the pipeline must be
  formally validated (run a known reference genome, "GIAB", and check the accuracy numbers).

---

## 6. Mini-glossary

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

---

## Files in this repo
- `README.md` — this guide
- `gcb.config` — the cloud settings (project, bucket, Spot, 4 callers)
- `env.sh` — loads Java + Nextflow into your terminal
- `samplesheet.example.csv` — template for listing your input files
- `consensus.sh` — *(coming next)* builds the ≥2-caller / DeepVariant-priority variant list
