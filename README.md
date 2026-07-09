# Sarek Clinical Pipeline — Plain-Language Guide

Clinical germline variant calling on **Google Cloud** using **four independent variant callers**,
keeping a variant when **at least two agree** (DeepVariant given priority), then handing a single
consensus VCF to the `candidate-filtering` repo. The heavy computing runs in the cloud, not on your laptop.

> **Two entry points** (same callers, same `consensus.sh`, same Google Cloud setup):
> - **WGS from FASTQ** — this guide, below.
> - **BGE exome from CRAM** (Terra Blended Genome-Exome, exome arm) — see **[BGE.md](BGE.md)**.

This guide (the WGS path) assumes almost no cloud experience. Read it top to bottom once.

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

**Watch progress *and* live cost** (BGE cohort arm) without paying anything to look —
the monitors only *list* the bucket/Batch jobs (no compute, no egress):
```bash
watch -n 30 ~/sarek-clinical/bge_dashboard.sh          # progress bars + Spot cost/budget bar + projected total
BUDGET=30 watch -n 30 ~/sarek-clinical/bge_dashboard.sh # set your own budget ceiling for the bar
```
The cost bar reconstructs accrued Spot spend from the Batch job records and shows it
against a budget with a projected final cost — so a run can never quietly run past what
you expected. (`bge_cost.sh` is the cost bar alone; `bge_progress.sh` /
`bge_filter_progress.sh` are the calling and VEP/filter progress bars.)

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

## 6b. Running a big cohort reliably — checklist (hard-won)

Do these **before** launching a large parallel run (e.g. 20+ exomes) so it finishes in one go
without hand-holding:

1. **Raise the region CPU quota first.** Default `us-central1` `CPUS` and `N2D_CPUS` are **200** —
   far too low for a cohort. A run will still work (jobs queue/retry) but crawls and floods the log
   with `CODE_GCE_QUOTA_EXCEEDED`. Raise both to **≥1000** (free — a quota is a ceiling, not a charge):
   *Console → IAM & Admin → Quotas →* filter "CPUs" / "N2D CPUs", region us-central1 → Edit. CPU bumps
   usually auto-approve in minutes. (`gcloud` has no quota subcommand; use the Console or the Cloud
   Quotas REST API.)
2. **Spot vs on-demand.** Default is **Spot** (cheap, preemptible). Long **HaplotypeCaller** jobs can be
   reclaimed faster than they finish and, if a task is preempted more than `maxRetries` times, the whole
   run aborts. `gcb.config` now retries up to **5×** and caps concurrency (`queueSize=40`). If a run keeps
   getting preempted, switch to on-demand: **`SAREK_SPOT=false nextflow run … -c gcb.config …`** (≈3× the
   VM cost but zero preemption — worth it for the tail of a stuck run).
3. **Skip fragile QC on big runs.** `vcftools`/`multiqc` QC steps can themselves fail on Spot and abort
   an otherwise-complete run. Add `--skip_tools baserecalibrator,vcftools,multiqc` (BQSR is unnecessary
   with a DeepVariant consensus backbone; `bcftools stats` QC still runs). Drop `baserecalibrator` from
   that list if you want strict GATK-best-practice parity.
4. **Uploading FASTQ from Windows/WSL (`/mnt/c`):** disable parallel composite uploads or large files
   corrupt ("Temporary components were not uploaded correctly"):
   `export CLOUDSDK_STORAGE_PARALLEL_COMPOSITE_UPLOAD_ENABLED=False` before `gcloud storage rsync/cp`.
5. **If the driver dies, cancel orphan Batch jobs.** A killed Nextflow driver does **not** stop its
   Batch jobs — they keep running (and billing). Check `gcloud compute instances list --filter="name~^nf-"`
   and delete stragglers by exact job id before resuming, or they re-saturate the quota.
6. **Downstream naming for `candidate-filtering`.** Name samples `<FAMILY>-P/-M/-F` for trio/duo analysis.
   Plainly-named singleton cohorts (e.g. `EPIGEN01..20`) now auto-run as singletons in `filtering_r.pl`
   (it prints a `NOTE:`), so no `--proband` is needed — but they get `inheritance=NA`.

### Third entry point: **exome from FASTQ** (EPIGEN-style)
Besides WGS-from-FASTQ (this README) and BGE-exome-from-CRAM (`BGE.md`), the repo now supports **exome
from FASTQ** — full alignment + calling, exome-scoped. Launcher: `run_epigen_wes.sh`
(`--step mapping --wes --intervals <kit BED> --skip_tools baserecalibrator`, 4 callers). Provide the
capture-kit target BED (GRCh38, chr-prefixed; e.g. Agilent SureSelect V6 `S07604514` padded). Then
`consensus_from_results.sh` → `candidate-filtering`. Cost for 20 small exomes ≈ **$15–30** on Spot.

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

**Shared core**
- `consensus.sh` — union consensus: all DeepVariant calls + variants ≥2 other callers agree on (genotype borrowed from Strelka2/HaplotypeCaller), tagged with `CALLERS`/`NCALLERS`/`CONF`/`GT_SOURCE` (Section 5)
- `env.sh` — loads Java + Nextflow into your terminal (and sets `NXF_SYNTAX_PARSER=v1`, required for sarek 3.8.1 on Nextflow 26.x)

**WGS from FASTQ** (this guide)
- `README.md` — this guide
- `gcb.config` — cloud settings for real runs (project, bucket, Spot, 4 callers)
- `gcb-smoke.config` — cloud settings for the tiny test run (no `params`, lets the `test` profile drive)
- `samplesheet.example.csv` — template for listing your input FASTQs

**BGE exome from CRAM** (see [BGE.md](BGE.md))
- `gcb-bge-wes.config` — Google Batch profile: Sarek `--step variant_calling --wes` from CRAM
- `make_samplesheet.sh` — family table → Sarek samplesheet (`<family>-<role>` naming)
- `run_bge_wes.sh` — launch the BGE exome calling on Batch
- `consensus_from_results.sh` — pull per-caller VCFs and run `consensus.sh` per sample (sample column auto-detected by header)
- `run_bge_annotate_filter.sh` — VEP-annotate consensus VCFs + run candidate-filtering (`OUT_NAME` env sets the output folder)
- `families.example.tsv` — template family/CRAM table

**Exome from FASTQ** (align + call; Section 6b)
- `run_epigen_wes.sh` — launch exome calling from FASTQ (`--step mapping --wes`, kit BED, on-demand-capable via `SAREK_SPOT`)
- `upload_epigen_fastq.sh` — resumable FASTQ → bucket upload (parallel-composite disabled for WSL/`/mnt/c` reliability)
- Monitoring: `bge_cost.sh` (live Spot/on-demand cost vs budget), `bge_dashboard.sh`, `bge_progress.sh` / `bge_filter_progress.sh`

**Reliability defaults** (in `gcb.config`)
- `queueSize=40` concurrency cap · `maxRetries=5` (survives Spot preemption streaks) · `SAREK_SPOT=false` → on-demand VMs
