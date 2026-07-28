# Sarek Clinical Pipeline — Plain-Language Guide

Clinical germline variant calling on **Google Cloud** using **four independent variant callers**,
keeping a variant when **at least two agree** (DeepVariant given priority), then handing a single
consensus VCF to the `candidate-filtering` repo. The heavy computing runs in the cloud, not on your laptop.

> **Two entry points** (same callers, same `consensus.sh`, same Google Cloud setup):
> - **WGS from FASTQ** — this guide, below.
> - **BGE exome from CRAM** (Terra Blended Genome-Exome, exome arm) — see **[BGE.md](BGE.md)**.

This guide (the WGS path) assumes almost no cloud experience. Read it top to bottom once.

> **New here, setting this up from scratch?** Start at **[Section 0](#0-first-time-setup-starting-from-nothing)**
> — install, Google Cloud project, bucket, permissions, reference genome. Already using the
> provisioned `intergenica` setup? Skip straight to [Section 3](#3-how-to-run-it-step-by-step).

---

## 0. First-time setup (starting from nothing)

*Skip this entire section if you are using the already-provisioned setup in Section 2. It is here for
someone standing up their own copy — roughly 30 minutes plus a large download.*

**You will need:** a Google account with billing enabled, and ~10 GB of local disk for the reference
genome. Everything below is a one-time cost.

### 0.1 — Clone the repo (anywhere you like)

```bash
git clone https://github.com/edoper/sarek-clinical.git
cd sarek-clinical
```

Paths are derived from where the scripts actually are, so any location works — there is no required
directory.

### 0.2 — Install the tools

| Tool | Why | Install |
|---|---|---|
| **Java 17+** | Nextflow runs on the JVM | `sudo apt install openjdk-21-jdk` (or [SDKMAN](https://sdkman.io)) |
| **Nextflow** | sends jobs to the cloud | `curl -s https://get.nextflow.io \| bash && sudo mv nextflow /usr/local/bin/` |
| **gcloud CLI** | talks to Google Cloud | [cloud.google.com/sdk/docs/install](https://cloud.google.com/sdk/docs/install) |
| **bcftools, htslib, samtools** | the local consensus step | `conda install -c bioconda bcftools htslib samtools` |
| **python3** | `build_cohort.py` | usually already present |

Docker is **not** needed locally — containers run on the cloud VMs, not your laptop.

If your Java lives somewhere unusual, set `JAVA_HOME` in `site.env` (step 0.4); `env.sh` defaults it
to `$HOME/jdk21`.

### 0.3 — Create the Google Cloud project, bucket and permissions

Replace `MY-PROJECT` and `MY-BUCKET` throughout. These are the exact roles and APIs Google Batch
requires — no more, so you are not granting anything broad.

```bash
# a) project + APIs
gcloud projects create MY-PROJECT                  # or reuse an existing project
gcloud config set project MY-PROJECT
gcloud services enable batch.googleapis.com compute.googleapis.com \
                       logging.googleapis.com storage.googleapis.com
#    ...then link a billing account (Console → Billing) — nothing runs without it.

# b) the bucket, in the SAME region you will compute in (avoids egress charges)
gcloud storage buckets create gs://MY-BUCKET --location=US-CENTRAL1 --uniform-bucket-level-access

# c) let the Batch VMs report status, write logs, and read/write the bucket
PROJNUM=$(gcloud projects describe MY-PROJECT --format='value(projectNumber)')
SA="${PROJNUM}-compute@developer.gserviceaccount.com"
gcloud projects add-iam-policy-binding MY-PROJECT --member="serviceAccount:$SA" --role=roles/batch.agentReporter
gcloud projects add-iam-policy-binding MY-PROJECT --member="serviceAccount:$SA" --role=roles/logging.logWriter
gcloud storage buckets add-iam-policy-binding gs://MY-BUCKET --member="serviceAccount:$SA" --role=roles/storage.objectAdmin

# d) let YOUR account submit jobs on that service account's behalf
#    (skip if you are the project Owner — Owner already covers these)
gcloud projects add-iam-policy-binding MY-PROJECT --member="user:you@example.com" --role=roles/batch.jobsEditor
gcloud projects add-iam-policy-binding MY-PROJECT --member="user:you@example.com" --role=roles/iam.serviceAccountUser
gcloud projects add-iam-policy-binding MY-PROJECT --member="user:you@example.com" --role=roles/logging.viewer

# e) give Nextflow your credentials
gcloud auth application-default login
gcloud auth application-default set-quota-project MY-PROJECT
```

**Raise the CPU quota before your first cohort.** Default `us-central1` `CPUS` and `N2D_CPUS` are 200,
far too low; raise both to ≥1000 (free — a quota is a ceiling, not a charge). See
[Section 6b](#6b-running-a-big-cohort-reliably--checklist-hard-won), which also explains the
`IN_USE_ADDRESSES` limit that caps concurrency before CPUs do.

### 0.4 — Point the repo at your project

Create **`site.env`** next to `site.sh`. It is never committed — this is where your own settings live:

```bash
# site.env
SAREK_PROJECT=MY-PROJECT
SAREK_BUCKET=gs://MY-BUCKET
SAREK_REGION=us-central1
CF=$HOME/code/candidate-filtering       # downstream repo (default: beside this one)
WIN=                                    # WSL only: a Windows folder to copy results to; leave empty otherwise
```

Then:

```bash
source ./env.sh
```

It prints which project and bucket you are pointed at — check that line before your first run. Both
the shell scripts and the Nextflow configs read these values, so this one file retargets everything.
(Exporting the variables directly works too, and takes precedence over `site.env`.)

### 0.5 — Get the reference genome

Two copies are needed: one **in your bucket** (the cloud jobs read it) and one **on your laptop**
(the local consensus step reads it). Both come from Broad's public bucket — free, and the same files
this pipeline was developed against.

```bash
BROAD=gs://gcp-public-data--broad-references/hg38/v0

# a) into your bucket (~4.8 GB, server-side copy, takes a few minutes)
gcloud storage cp \
  $BROAD/Homo_sapiens_assembly38.fasta \
  $BROAD/Homo_sapiens_assembly38.fasta.fai \
  $BROAD/Homo_sapiens_assembly38.dict \
  $BROAD/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  $BROAD/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi \
  gs://MY-BUCKET/refs/GATK.GRCh38/

# b) onto your laptop for consensus.sh (~3.3 GB + index)
mkdir -p refs
gcloud storage cp $BROAD/Homo_sapiens_assembly38.fasta     refs/
gcloud storage cp $BROAD/Homo_sapiens_assembly38.fasta.fai refs/
```

> **It must be the full GATK GRCh38, with ALT and decoy contigs — not a GENCODE "primary assembly".**
> `bcftools norm` fails on any call landing on an ALT contig, and exome targets do include them. This
> is the single most common way to break the consensus step.

### 0.6 — Get `candidate-filtering` (the downstream half)

This repo stops at a consensus VCF. Turning that into a ranked candidate list is a separate repo:

```bash
git clone https://github.com/edoper/candidate-filtering.git
```

Clone it **beside** `sarek-clinical` and it is found automatically; otherwise set `CF` in `site.env`.
It has its own prerequisites (Ensembl VEP plus its cache and plugins) — see its README. You can run
everything in this guide without it, and only need it for Section 5's hand-off.

### 0.7 — Check it works

```bash
./test/test_consensus.sh     # ~4s, synthetic data, no cloud, no cost
```

`ALL TESTS PASSED` means the local half (bcftools, htslib, the reference logic) is sound. For the
cloud half, the cheapest real check is the smoke profile, which runs Sarek's own tiny built-in test
dataset on Batch for a few cents:

```bash
nextflow run nf-core/sarek -r 3.8.1 -profile test,docker -c gcb-smoke.config
```

If that finishes, your project, bucket, permissions and quota are all correct — go to Section 3.

---

## 1. The mental model (read this first)

Think of it like ordering food delivery:

| Real life | This pipeline |
|---|---|
| You, on your phone | **Your laptop** — only gives orders, does no cooking |
| The restaurant kitchen | **Google Cloud** — rents powerful computers that do the actual work |
| The fridge where ingredients live | **The bucket** (`$SAREK_BUCKET`) — cloud storage for your files |
| The recipe | **Sarek** — the standard, published analysis workflow |
| The waiter taking your order to the kitchen | **Nextflow** — the program on your laptop that sends jobs to the cloud |

Key idea: **your laptop stays light.** You type one command; Google spins up big
computers, runs the analysis, writes the results into the bucket, and shuts the
computers down automatically. You only pay while they run.

**Three names you'll see a lot:**
- **Bucket** = a cloud folder. Paths start with `gs://`. This deployment's is `gs://intergenica-sarek-clinical`; the commands below use `$SAREK_BUCKET` so they work whichever one is yours.
- **Google Batch** = the service that rents computers on demand and turns them off when done.
- **Spot** = cheap "leftover" computers (~70–90% off). They can occasionally be taken back
  mid-job; the pipeline just automatically retries, so you save money safely.

---

## 2. What's already set up

If you are working on **this** deployment, all of the below already exists and you can skip to
Section 3:

- **Google project:** `intergenica`  **Billing:** the "Computacion-nube" account
- **Bucket:** `gs://intergenica-sarek-clinical` with three folders:
  - `fastq/` → where you put input files
  - `work/` → scratch space the pipeline uses while running (delete it afterwards)
  - `results/` → your output (variants + quality reports)
- **The laptop (WSL):** Java + Nextflow installed; settings file `gcb.config` ready.

**Setting this up somewhere else for the first time? → [Section 0](#0-first-time-setup-starting-from-nothing).**
Those names above are only this deployment's defaults; nothing in the repo is tied to one machine,
one person, or one Google project.

---

## 3. How to run it, step by step

### Step 0 — open the toolbox (every time you start a terminal)
```bash
source ./env.sh          # from the repo directory (or use its full path)
```
*Why:* loads Java + Nextflow **and** your project/bucket settings, so the commands below work.

### Step 1 — log in to Google (only the very first time)
```bash
gcloud auth application-default login
gcloud auth application-default set-quota-project "$SAREK_PROJECT"
```
*Why:* gives Nextflow permission to rent cloud computers on your behalf. A browser opens; pick your Google account.

### Step 2 — upload your FASTQ files to the bucket
```bash
gcloud storage cp *_R1.fastq.gz *_R2.fastq.gz "$SAREK_BUCKET/fastq/"
```
*Why:* the cloud computers read from the bucket, not from your laptop.

### Step 3 — make a samplesheet (a small table listing your files)
Copy `samplesheet.example.csv`, rename it `samplesheet.csv`, and edit the file paths.
Columns: `patient,sample,lane,fastq_1,fastq_2`. One row per pair of FASTQ files.

### Step 4 — run the pipeline
```bash
cd /path/to/sarek-clinical
nextflow run nf-core/sarek -r 3.8.1 -profile docker -c gcb.config \
  --input  samplesheet.csv \
  --outdir "$SAREK_BUCKET/results/run01" \
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
watch -n 30 ./bge_dashboard.sh          # progress bars + Spot cost/budget bar + projected total
BUDGET=30 watch -n 30 ./bge_dashboard.sh # set your own budget ceiling for the bar
```
The cost bar reconstructs accrued Spot spend from the Batch job records and shows it
against a budget with a projected final cost — so a run can never quietly run past what
you expected. (`bge_cost.sh` is the cost bar alone; `bge_progress.sh` /
`bge_filter_progress.sh` are the calling and VEP/filter progress bars.)

### Step 5 — get your results
Results land in `$SAREK_BUCKET/results/run01/`. Download the variant
files (VCF) and the quality report (MultiQC) when ready:
```bash
gcloud storage cp -r "$SAREK_BUCKET/results/run01" ./run01-results
```

### Step 6 — clean up to stop paying for storage
```bash
gcloud storage rm -r "$SAREK_BUCKET/work"
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
source ./env.sh                       # puts tools on PATH (if you installed them there)

./consensus.sh \
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

> **If it is interrupted, it leaves nothing behind.** The VCF is built in a `.partial` file and only
> renamed into place once complete and indexed, so a Ctrl-C or shutdown leaves *no* output rather than
> a truncated one — and re-running redoes that sample instead of trusting a half-written file. See
> Section 8.

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

- **Cost: budget ~$5–6 of Spot compute per WGS genome** (align + 4 callers) — so ~**$25** for
  4 genomes, ~**$45–50** for 8. Both figures are **measured**, not estimated, reconstructed from
  Google Batch job records with `bge_cost.sh` (which prices any run from the control plane at no
  charge):

  | run type | samples | vCPU-hours | total | **per sample** |
  |---|---|---|---|---|
  | WGS from FASTQ, single | 1 | 518 | $5.67 | **$5.67** |
  | WGS from FASTQ, cohort | 8 | 4,194 | ~$45 | **~$5.46** |
  | Exome from FASTQ | 20 | 834 | $8.14 | **$0.41** |
  | Exome from CRAM (BGE, calling only) | 63 | 1,009 | $13.68 | **$0.22** |
  | Mixed exome/BGE cohort | 56 | 577 | $6.90 | **$0.12** |

  Read that last column: **an exome costs ~13x less than a genome, and starting from an aligned
  CRAM instead of FASTQ roughly halves what remains.** Choose the entry point before optimising
  anything else.

- **Budget per GENOME, not per gigabase — this is the easy mistake.** The 8-genome cohort used
  **8.03× the vCPU-hours of the single genome**, matching the *sample count* (8×) almost exactly,
  **not** the FASTQ volume (5.8×). Reason: the four callers each traverse the whole 3.1 Gb genome
  once per sample no matter how deep the reads are, so **calling is a fixed cost per genome** and
  only alignment scales with data volume. Estimating a cohort by scaling total FASTQ size
  under-budgets it by ~40%.
- **Spot preemption rework costs ~10–15%.** In the 8-genome cohort, 176 tasks were reclaimed mid-run
  and redone: 495 of 4,160 vCPU-hours, **$5.42 (13%)**, thrown away. That is the price of Spot being
  ~5× cheaper than on-demand, and it is a good trade — but include it when budgeting, and expect the
  share to rise with cohort size, since more concurrent VMs means more exposure to reclaims.
- *(Historical note: revisions of this file before 2026-07 estimated **$50–100 per 4 genomes**, then
  briefly claimed cost scales with FASTQ volume. Both were wrong — the first by ~4× high, the second
  in its scaling law. Budget from the measured per-genome table above.)*
- **Where the money actually goes (measured, one WGS genome, from the Nextflow execution trace):**

  | step | CPU-hours | share |
  |---|---|---|
  | DeepVariant | 124.3 | **52%** |
  | BWA-MEM alignment | 69.3 | **29%** |
  | HaplotypeCaller | 25.1 | 10% |
  | Strelka2 | 7.9 | 3% |
  | FreeBayes | 2.3 | **1%** |
  | everything else (QC, markdup, merges) | 10.0 | 4% |

  Two counter-intuitive consequences. **FreeBayes is essentially free** — dropping it to save money
  would buy ~1% while losing a caller that contributes to the >=2-caller rescue, so don't. And
  **DeepVariant alone is over half the bill**, but it is the consensus backbone, so it is the one
  caller you cannot drop. The genuinely large, zero-compromise lever is **alignment (29%)**: it
  disappears entirely if you can start from an aligned GRCh38 CRAM instead of FASTQ.
- **You only pay while jobs run** + a small amount for files sitting in the bucket.
- The biggest waste is forgetting to delete the `work/` folder — do Step 6. For WGS this
  dominates: `work/` runs ~200 GB **per genome**, so an 8-genome cohort left behind costs
  ~$30–60/month in storage — several times the compute that produced it.
- **Validation status — DONE (WGS 2026-07-27, EXOME 2026-07-28):** the pipeline has been run against the GIAB **HG002**
  reference genome and scored with RTG `vcfeval` against the **GIAB v4.2.1** truth set.
  **Genome-wide F1 = 0.9948** (SNV 0.9949, INDEL 0.9941) inside GIAB high-confidence regions for the
  default union consensus — and **0.9959** (SNV 0.9961, INDEL 0.9946) for the DeepVariant backbone
  alone. The ≥2-caller rescue arm is **net-negative genome-wide** (8.7 false positives per true variant
  recovered) but **neutral inside the g4e panel** (identical true variants found, 3 extra false
  positives across 1,066 genes) — so the default is fine as shipped. See the results file. Full numbers, per-tier comparison
  and limitations: **[validation/RESULTS-HG002.md](validation/RESULTS-HG002.md)**; how to re-run:
  [validation/README.md](validation/README.md). The validation cost **$4.68**.
  **The exome arm was validated separately and performs worse — quote the right number.** On a real
  capture library (GIAB HG002 exome), in the same panel genes: SNV F1 **0.9955**, but INDEL F1
  **0.9198** (vs 0.9934 for WGS reads merely restricted to the panel). Indel sensitivity is the
  limiting factor for exome/BGE reporting and belongs in report limitations — see
  **[validation/RESULTS-HG002-EXOME.md](validation/RESULTS-HG002-EXOME.md)**.
  This measures the *variant-calling* pipeline only, inside high-confidence regions, for SNVs and
  small indels — it does not cover structural variants, and it does not validate the downstream
  `candidate-filtering` interpretation step. Local clinical governance still applies before use in
  patient care.

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
   Plainly-named singleton cohorts (e.g. `SAMPLE01..NN`) now auto-run as singletons in `filtering_r.pl`
   (it prints a `NOTE:`), so no `--proband` is needed — but they get `inheritance=NA`.

### Third entry point: **exome from FASTQ**
Besides WGS-from-FASTQ (this README) and BGE-exome-from-CRAM (`BGE.md`), the repo now supports **exome
from FASTQ** — full alignment + calling, exome-scoped. Launcher: `run_epigen_wes.sh`
(`--step mapping --wes --intervals <kit BED> --skip_tools baserecalibrator`, 4 callers). Provide the
capture-kit target BED (GRCh38, chr-prefixed; e.g. Agilent SureSelect V6 `S07604514` padded). Then
`consensus_from_results.sh` → `candidate-filtering`. Cost ≈ **$1–2 per small exome** on Spot.

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

## 8. Checking the consensus step still works

`consensus.sh` is the only analysis logic this repo owns (everything else is standard Sarek), and a
mistake in it produces a **plausible-looking VCF rather than an error** — so nothing downstream would
notice. There is a test that pins its behaviour:

```bash
source ~/sarek-clinical/env.sh
~/sarek-clinical/test/test_consensus.sh        # ~4 seconds
```

It makes up its own tiny genome and four fake caller files — **no cloud, no patient data, no network** —
and checks two things:

1. **The rule is still the rule.** A variant all four callers found comes out `HIGH` with DeepVariant's
   genotype; a DeepVariant-only variant comes out `LOW`; a variant DeepVariant missed but three others
   found is rescued with Strelka2's genotype; one only HaplotypeCaller and FreeBayes found is rescued
   with HaplotypeCaller's; and variants only *one* non-DeepVariant caller found are correctly dropped.
2. **A crash can't produce a half-finished result.** The consensus VCF is written to a temporary
   `.partial` file and only renamed into place once it is complete and indexed. So if the run is
   interrupted (Ctrl-C, laptop shutdown), there is simply **no output file** — rather than a truncated
   one that a later re-run would mistake for finished work and pass on to VEP and candidate-filtering.
   The test proves this by killing a large run mid-write and confirming the final file never appears.

Run it after any change to `consensus.sh` or `consensus_from_results.sh`. It prints `ALL TESTS PASSED`
or names exactly what broke.

---

## Files in this repo

**Shared core**
- `LICENSE` — MIT (with a note that clinical use requires validation first)
- `site.sh` — one place for project/bucket/local paths (`SAREK_*`, `CF`, `WIN`); override in an untracked `site.env` (Section 2)
- `consensus.sh` — union consensus: all DeepVariant calls + variants ≥2 other callers agree on (genotype borrowed from Strelka2/HaplotypeCaller), tagged with `CALLERS`/`NCALLERS`/`CONF`/`GT_SOURCE` (Section 5)
- `env.sh` — sources `site.sh`, then loads Java + Nextflow into your terminal (and sets `NXF_SYNTAX_PARSER=v1`, required for sarek 3.8.1 on Nextflow 26.x)
- `test/test_consensus.sh` — regression test for `consensus.sh` (Section 8)

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

**Validation**
- `validation/` — GIAB HG002 accuracy validation (WGS **and** exome): `run_giab.sh` (end to end), `benchmark_giab.sh` (RTG vcfeval per confidence tier), `health_check.sh` (progress/cost + budget guard), `RESULTS-HG002.md` / `RESULTS-HG002-EXOME.md` (the measured numbers), `SOP.md`, `EFFICIENCY.md`, `../qc_gate.sh` (per-sample gate)

**Reliability defaults** (in `gcb.config`)
- `queueSize=40` concurrency cap · `maxRetries=5` (survives Spot preemption streaks) · `SAREK_SPOT=false` → on-demand VMs
