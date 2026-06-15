# BGE Exome Arm — Variant Calling from CRAM

Entry point of this repo for **Blended Genome-Exome (BGE)** CRAMs, focused on the **exome arm**
(the deep ~60× portion). The low-pass genome arm (imputation/GLIMPSE) is **out of scope**. The other
entry point — WGS from FASTQ — is the main [README](README.md); both share `consensus.sh`, `env.sh`,
the bucket, and the Google Cloud setup.

```
 Terra BGE CRAMs ─► Sarek WES variant-calling (Google Batch, from CRAM, 4 callers)
                       DeepVariant · Strelka2 · FreeBayes · HaplotypeCaller
                 ─► per-caller VCFs (bucket)
                 ─► consensus.sh  (DeepVariant backbone + ≥2-caller rescue)        ← small, local
                 ─► <family>-<role>.consensus.vcf.gz
                 ─► vep_annotate.sh + filtering_r.pl  (candidate-filtering repo)   ← local
```

**Design (this arm):** everything heavy runs **in the cloud** next to the Terra data; only small VCFs
come local. **Start from CRAM** (`--step variant_calling`, no re-align/BQSR — the Terra CRAMs are already
processed). **Call broadly across the confident exome regions, filter later** (intervals are the exome
calling targets, not MANE-restricted, so you can re-explore without re-calling).

---

## Files (this arm)

| File | Purpose |
|---|---|
| `gcb-bge-wes.config` | Google Batch profile: Sarek `--step variant_calling --wes`, 4 callers, GATK.GRCh38, Spot, `maxRetries=4`, **pre-staged reference** (`gs://…/refs/GATK.GRCh38/`). |
| `build_cohort.py` | Terra sample-table export (`*.tsv`) → `families.tsv` + a CRAM staging list. Role map: proband=no suffix(`-P`), `…M`=madre(`-M`), `…P`=padre(`-F`). |
| `make_samplesheet.sh` | `families.tsv` → Sarek samplesheet, naming samples `<family>-<role>` for candidate-filtering. |
| `run_bge_wes.sh` | Launch the calling on Batch (single sample / small set). |
| `consensus_from_results.sh` | Pull the small per-caller VCFs, run `consensus.sh`, emit `<sample>.consensus.vcf.gz`. |
| `bge_progress.sh` | Zero-cost progress bar (reads the local log + lists bucket objects). |
| `families.example.tsv` | Template family/CRAM table (for the manual `make_samplesheet.sh` path). |

---

## Cohort workflow (validated — Terra sample-table export → N samples)

Recommended for a real batch. Export the Terra `sample` table to a `.tsv` (it carries
`collaborator_sample_id`, `genome_cram_path`, `genome_crai_path`, `predicted_sex`).

```bash
source ~/sarek-clinical/env.sh
# 1) Terra export -> families.tsv (staged paths) + /tmp/cram_srcs.txt (sources to copy)
python3 build_cohort.py /path/to/terra_samples_export.tsv
# 2) Stage CRAMs TDR->bucket (parallel, server-side) and build the samplesheet
gcloud storage cp -I gs://intergenica-sarek-clinical/bge-wes/crams-cohort/ < /tmp/cram_srcs.txt
./make_samplesheet.sh families.tsv > samplesheet-cohort.csv
# 3) Launch all samples in ONE run (coding-only Twist 35Mb intervals; ref pre-staged via config)
nextflow run nf-core/sarek -r 3.8.1 -profile docker -c gcb-bge-wes.config \
  --step variant_calling --input samplesheet-cohort.csv \
  --intervals gs://intergenica-sarek-clinical/bge-wes/targets/twist_coding_targets.bed \
  --outdir   gs://intergenica-sarek-clinical/bge-wes/results-cohort \
  --tools deepvariant,strelka,freebayes,haplotypecaller --genome GATK.GRCh38 \
  -work-dir  gs://intergenica-sarek-clinical/bge-wes/work-cohort -ansi-log false -resume
# 4) Monitor (zero cost):  watch -n 60 ~/sarek-clinical/bge_progress.sh
# 5) Consensus for all samples, then VEP + candidate-filtering
SAMPLESHEET=samplesheet-cohort.csv OUTDIR=gs://intergenica-sarek-clinical/bge-wes/results-cohort \
  ./consensus_from_results.sh
```

**Intervals — call coding-only (Twist 35 Mb), not the broad 165 Mb BGE region.** BGE is deep
(~110×) only on the exome targets; the rest of the 165 Mb is low-pass (~3–13×), so direct-calling
there yields low-confidence noise that even passes the MANE/consequence gate as `lowDP`. The low-pass
genome is properly used by **imputation (GLIMPSE, arm 2)** — not by widening the calling region.
**Cost ≈ $0.10/sample** (coding-only, Spot). Raise the `CPUs` quota (us-central1) for faster wall-clock.

---

## Inputs you must provide

1. **CRAM locations** — Terra BGE CRAMs as `gs://` URIs + family/role/sex → `families.tsv` (see `families.example.tsv`).
2. **Exome target intervals** — the BGE **exome calling regions** (BED/`.interval_list`; the Broad BGE
   workspace ships one). Upload to `gs://intergenica-sarek-clinical/bge-wes/targets/`.
3. **Cloud read access** to the CRAMs — see Step 1.

---

## Cloud runbook

> Prereqs already done from the WGS setup: `gcloud` authenticated, ADC + quota project set, bucket +
> Batch SA exist. `source ~/sarek-clinical/env.sh` each shell.

```bash
# Step 1 — make Terra CRAMs reachable from Batch (copy into the intergenica bucket; same region).
gcloud storage cp -u intergenica \
  gs://fc-YOUR-WORKSPACE-BUCKET/path/FAM01_proband.cram* \
  gs://intergenica-sarek-clinical/bge-wes/crams/
#   (or grant 906897002329-compute@developer.gserviceaccount.com objectViewer on the Terra bucket)

# Step 2 — upload exome intervals
gcloud storage cp bge_calling_regions.bed gs://intergenica-sarek-clinical/bge-wes/targets/

# Step 3 — samplesheet
cp families.example.tsv families.tsv          # edit with your real gs:// CRAM paths
./make_samplesheet.sh families.tsv > samplesheet.csv

# Step 4 — run Sarek WES on Google Batch
export INTERVALS=gs://intergenica-sarek-clinical/bge-wes/targets/bge_calling_regions.bed
./run_bge_wes.sh

# Step 5 — consensus (small, local)
./consensus_from_results.sh                   # -> ./consensus/<family>-<role>.consensus.vcf.gz

# Step 6 — annotate + filter (candidate-filtering repo)
cd ~/candidate-filtering
bash vep_annotate.sh ~/sarek-clinical/consensus/FAM01-P.consensus.vcf.gz FAM01-P.germline.vep.vcf.gz
bash run_filtering.sh                          # -> FAM01-P.g4e-2025.candidatos

# Step 7 — clean up bucket scratch
gcloud storage rm -r gs://intergenica-sarek-clinical/bge-wes/work
```

---

## Notes
- **Reference match:** BGE CRAMs are GRCh38 (`Homo_sapiens_assembly38`); Sarek decodes them with
  `--genome GATK.GRCh38` (igenomes, cloud-side). The local consensus step uses the chr-named GRCh38
  already on disk, correct for exome candidates on the main chromosomes.
- **Cost:** Spot + exome-only + alignment-free keeps per-sample cost low; scale ≈ samples × 4 callers.
- **Local-only alternative** (not recommended at scale) needs a container engine (Apptainer), the GATK
  bundle locally (~30–50 GB), and more WSL disk (Windows-side `wsl --manage <distro> --resize`, which
  cannot be done from inside WSL). Validate against GIAB before clinical use.
