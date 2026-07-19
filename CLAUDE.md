# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Clinical **germline variant calling** on Google Cloud, feeding the separate `~/candidate-filtering`
repo. Heavy work (nf-core/sarek via Nextflow on **Google Batch**) runs in the cloud; only small VCFs
come local. There is no build/test suite — this is a collection of Bash/Python orchestration scripts
plus Nextflow config, driven from WSL. `README.md` (WGS) and `BGE.md` (BGE exome) are the
plain-language user guides; read them for the end-to-end story before changing behavior.

## The core idea: 4 callers → union consensus → filtering

Every arm runs the same four callers (**DeepVariant, Strelka2, FreeBayes, HaplotypeCaller**) and
combines them with **`consensus.sh`** — the single shared brain of the repo, reused (never duplicated)
by every arm. Its consensus rule is a **union, not a majority filter**:

- **Backbone** = *every* DeepVariant call, keeping DeepVariant's genotype fields (GT/GQ/DP/AD/VAF).
- **Rescue** = variants DeepVariant missed but **≥2 other callers** agreed on; genotype borrowed from
  Strelka2, else HaplotypeCaller (**FreeBayes counts toward concordance but is never a genotype donor**).
- Nothing is tiered away here. Each variant is tagged `CALLERS` / `NCALLERS` / `CONF`
  (HIGH≥3 / MEDIUM=2 / LOW=1) / `GT_SOURCE`. **Strictness is decided downstream** in candidate-filtering,
  by filtering on those tags. `CONF` measures cross-caller *agreement*, not correctness.

Pipeline inside `consensus.sh`: per-caller VCFs are normalized (keep `PASS`, `norm -m -both` split +
left-align against REF, `--rm-dup exact`), intersected with `bcftools isec -n +1`, then an awk pass
builds the presence/annotation VCF that drives backbone annotation + rescue. Two subtle invariants are
load-bearing and commented in-file: (1) `--rm-dup exact` **not** `all` (all collapses by position and
silently drops just-split alt alleles); (2) `harmonize()` computes an explicit drop-list from the
header instead of `-x ^keep,...` because bcftools 1.13's caret form errors when nothing needs removing.

## Three entry points (same callers, same consensus.sh, same cloud setup)

| Arm | Input | Sarek step | Config | Launcher | Guide |
|---|---|---|---|---|---|
| WGS from FASTQ | FASTQ | full (`mapping`) | `gcb.config` | manual `nextflow run` (README §4) | `README.md` |
| BGE exome from CRAM | Terra CRAM | `variant_calling --wes` (no re-align/BQSR) | `gcb-bge-wes.config` | `run_bge_wes.sh` | `BGE.md` |
| Exome from FASTQ | FASTQ | `mapping --wes` | `gcb.config` | `run_epigen_wes.sh` | README §6b |

After calling, all arms converge: **`consensus_from_results.sh`** (pull per-caller VCFs from bucket →
run `consensus.sh` per sample, resumable, sample column auto-detected) → **`run_bge_annotate_filter.sh`**
(VEP + candidate-filtering in an isolated workdir → `.candidatos` → copy to `$WIN`).

### Per-order scratch dirs (gitignored, but load-bearing)

A single order/cohort's working files — its samplesheet, launcher, unattended orchestrator, resume
scripts, and local consensus output — live in a **per-order scratch directory** that `.gitignore`
excludes (it holds real sample/family IDs and patient VCFs). Such a dir is **not a fourth arm**: it
reuses `gcb.config`/`gcb-bge-wes.config` + `consensus.sh` like everything else. The scripts below are a
**reusable template** — the same shapes recur every order and encode lessons worth preserving even
though the instances themselves stay private. Keep the private run-specific narrative (which cohort,
which IDs, what was found) in a per-run note under `$WIN`, never in tracked docs.

- **`run_<order>.sh {crams|fastq}`** — WGS variant: **no `--wes`, no `--intervals`**. CRAM-start and
  FASTQ-start must be **two separate `nextflow run`s** because sarek's `--step` is global (`crams` =
  `--step variant_calling`, `fastq` = `--step mapping`). `cd` to the repo root so the config resolves.
- **`orchestrate_<order>.sh`** — **unattended** end-to-end chain (stage/verify inputs → sarek → consensus
  → VEP + candidate-filtering → `$WIN/<order>-candidatos`). **Fail-closed**: any stage failure stops the
  chain. Every stage is **idempotent** (staging skips on exact byte-size match; sarek uses `-resume`),
  which is what makes re-running it the resume path. A `WAIT_SECS` env var can gate a pre-flight wait.
- **`resume_<order>.sh`** — full resume after a shutdown. Deletes orphan `nf-*` Batch jobs **first** (a
  dead driver leaves them billing), then re-execs the orchestrator with `WAIT_SECS=0`.
- **A local-only resume for the VEP/filter tail** — needed because of a real trap: a Stage-4 guard like
  `[[ ! -s VEPOUT ]]` treats a **partial** VEP VCF (killed mid-run) as "done" and would silently ship
  truncated candidatos. Such a script must **delete the partial VEP artifacts before re-annotating**.
  General rule: **any resume guard must invalidate partial outputs, not just test for existence.**
- **A `$WIN` log mirror** — `cp -f` (never append) the orchestrator log to `$WIN` on an interval, because
  drvfs append-caching hides live progress from the Windows side (see the WSL gotcha below).
- **Watch for placeholder samplesheets.** A committed `samplesheet-*.csv` example may hold `GS_CRAM_PATH_*`
  placeholders — never feed one to a real run.

## Fixed cloud environment (already provisioned — do not re-scaffold)

- GCP project `intergenica`, region **us-central1**, billing "Computacion-nube".
- Bucket `gs://intergenica-sarek-clinical` (`fastq/`, `work/`, `results/`; BGE under `bge-wes/`,
  EPIGEN under `epigen-wes/`). Reference pre-staged at `refs/GATK.GRCh38/` to skip the slow igenomes pull.
- Default execution is **Spot VMs** (~70–90% cheaper, preemptible, auto-retried). All work dirs live in
  the bucket, namespaced per arm so they never collide.

## Running things

**Always first:** `source ~/sarek-clinical/env.sh` (loads JDK21 + Nextflow and sets
`NXF_SYNTAX_PARSER=v1` — **required**; sarek 3.8.1 configs use legacy syntax the Nextflow 26.x parser
rejects). Pin **`-r 3.8.1`** and `-profile docker` on every `nextflow run`.

```bash
# BGE cohort, end to end (from BGE.md)
python3 build_cohort.py terra_export.tsv        # Terra table -> families.tsv + /tmp/cram_srcs.txt
gcloud storage cp -I gs://intergenica-sarek-clinical/bge-wes/crams-cohort/ < /tmp/cram_srcs.txt
./make_samplesheet.sh families.tsv > samplesheet-cohort.csv
nextflow run nf-core/sarek -r 3.8.1 -profile docker -c gcb-bge-wes.config \
  --step variant_calling --input samplesheet-cohort.csv --intervals <targets.bed> \
  --outdir gs://intergenica-sarek-clinical/bge-wes/results-cohort \
  --tools deepvariant,strelka,freebayes,haplotypecaller --genome GATK.GRCh38 \
  -work-dir gs://intergenica-sarek-clinical/bge-wes/work-cohort -ansi-log false -resume
SAMPLESHEET=samplesheet-cohort.csv OUTDIR=.../results-cohort \
  LOCAL_OUT=$HOME/sarek-clinical/consensus-cohort ./consensus_from_results.sh
./run_bge_annotate_filter.sh                     # -> $WIN/<OUT_NAME>/*.candidatos

# consensus for one sample, standalone
./consensus.sh -r refs/Homo_sapiens_assembly38.fasta -d DV.vcf.gz -o out/SAMPLE \
  -c strelka=STR.vcf.gz -c freebayes=FB.vcf.gz -c haplotypecaller=HC.vcf.gz

# zero-cost monitoring (control-plane listing only — no compute/egress)
watch -n 30 ~/sarek-clinical/bge_dashboard.sh    # progress + Spot cost/budget bar; BUDGET=30 to override
```

Scripts are env-var driven (override without editing): `SAMPLESHEET` `OUTDIR` `INTERVALS` `WORKDIR`
(launchers); `LOCAL_OUT` `REF` `CONSENSUS_SH` (consensus); `CONS_DIR` `CF` `WD` `WIN` `OUT_NAME`
(annotate/filter); `BUDGET` `SINCE` `REGION` `FRAC` (cost); `SAREK_SPOT` (see below).

## Non-obvious gotchas (would cost hours to rediscover)

- **Consensus REF must include ALT/decoy contigs.** Use `refs/Homo_sapiens_assembly38.fasta` (the full
  GATK.GRCh38), **not** a GENCODE primary-assembly FASTA — the latter fails `bcftools norm` on
  ALT-contig calls. This is why `consensus_from_results.sh` defaults `REF` to that file.
- **BGE: call coding-only (Twist ~35 Mb), never the broad 165 Mb BGE region.** BGE is deep (~110×) only
  on exome targets; the rest is low-pass (~3–13×) and direct-calling there yields noise that still
  passes the MANE/consequence gate as `lowDP`. The low-pass genome belongs to imputation (GLIMPSE, arm 2),
  which is **out of scope** here.
- **Spot preemption can abort a whole run — `errorStrategy` must cover the whole Batch `5000x` class.**
  Long jobs get reclaimed mid-run. All three configs retry the whole class `[8,10,14,50001..50007]` with
  `maxRetries` 6 (3 in smoke). **`50006` = "VM is recreated during task execution"** is the easy one to
  miss: with it absent, `errorStrategy` falls through to `finish` and a *single* reclaim aborts an entire
  cohort. Across a few hundred caller-runs a reclaim is near-certain, so a narrow list makes big Spot
  runs effectively impossible. Retrying `5000x` is safe: the task never executed (a failed workdir holds
  only `.command.sh`/`.command.run` — no `.command.err`, no output), so nothing partial can trip a retry. If a run still dies on the tail, rerun
  with **`SAREK_SPOT=false`** for on-demand VMs (~3× VM cost, zero preemption). Only `gcb.config` reads
  this env var; `gcb-bge-wes.config` hardcodes `spot=true`.
- **Two separate us-central1 quotas bite big cohorts, and `CODE_GCE_QUOTA_EXCEEDED` names both.**
  (1) CPU (`CPUS`, `N2D_CPUS`) — raise to ≥1000; default 200 starves parallelism.
  (2) **`IN_USE_ADDRESSES` (limit 69)** — *every* Batch VM takes an external IP, so concurrency is
  capped by IPs long before CPUs. This is why **every Batch config needs `executor { queueSize = 40 }`** —
  without it, one cohort can dispatch hundreds of jobs at once and saturate the IP quota. Quota is a free
  ceiling; there is no `gcloud` quota subcommand (use Console / Cloud Quotas REST API).
- **Big runs: `--skip_tools baserecalibrator,vcftools,multiqc`.** BQSR is unnecessary with a DeepVariant
  backbone; the fragile QC steps can themselves fail on Spot and abort an otherwise-complete run.
- **A killed Nextflow driver does NOT stop its Batch jobs** — they keep running and billing. Check
  `gcloud compute instances list --filter="name~^nf-"` and delete stragglers before `-resume`.
- **Uploading FASTQ from `/mnt/c` (WSL):** set `CLOUDSDK_STORAGE_PARALLEL_COMPOSITE_UPLOAD_ENABLED=False`
  or large files corrupt ("Temporary components were not uploaded correctly"). `upload_epigen_fastq.sh`
  already does this.
- **Sample naming is a downstream contract.** `<family>-P/-M/-F` drives candidate-filtering trio/duo
  auto-discovery. Plain singleton names (e.g. `EPIGEN01`) auto-run as singletons but get `inheritance=NA`.
  **Watch for Spanish role suffixes in source manifests**: a `…M`/`…P` pair can be *madre*/*padre* —
  a trio filed as three "probands". With no hyphen they run as unrelated singletons; rename to `-M`/`-F`
  (*padre*→`-F`, since `-P` is reserved for proband) and candidate-filtering computes real `inheritance`
  instead of `NA`. Always eyeball a manifest for families before launching.
- **Terra/TDR sources must be staged into our own bucket first.** TDR (`datarepo-*`) grants
  object-level read to the *user account* only — the Batch compute SA cannot read it, and even the user
  gets 403 on directory listing (exact-object `ls` works). So Batch can never read TDR directly; copy
  with user credentials into `gs://intergenica-sarek-clinical/<arm>/crams/`. Same-region, so ~$0 egress.
  **TDR exports often omit the `.crai`**: build them on Batch with a samtools
  container and a **gcsfuse volume** — no download, no egress. Mount under **`/mnt/disks/<name>`**;
  the Batch COS image has a read-only root, so any other `mountPath` fails with
  `Error while mounting gcsfuse: stat /mnt/bucket: no such file or directory`. `samtools quickcheck`
  in the same task doubles as an integrity gate on the staged copy.
- **`pkill -f <pattern>` self-kills** when the pattern appears in the very shell running it — killing
  `pkill -f nextflow` from a shell whose command line contains "nextflow" orphaned a live Nextflow JVM
  mid-run. Kill by PID (`ps -eo pid,cmd | grep …`). Same trap in `pgrep -f`-based liveness checks:
  they report a script as RUNNING when only the checking shell matches.
- **Nextflow's per-caller task counts are not file counts.** DeepVariant emits *both* `.vcf.gz` and
  `.g.vcf.gz` per sample, so counting `*.vcf.gz` double-counts it (a monitor read `76/56`). Count
  sample **directories** under `results/variant_calling/<caller>/` instead — exact and caller-agnostic.
- **`--skip_tools vcftools` does not skip `BCFTOOLS_STATS`.** The `VCF_QC_BCFTOOLS_VCFTOOLS` subworkflow
  still runs its bcftools half, which can fail on Spot and (without the `5000x` retry above) abort a
  finished run over a stats file.
- **`$WIN`** = `/mnt/c/Users/epere/Documents` — the Windows-side deliverable folder outputs are copied to.
