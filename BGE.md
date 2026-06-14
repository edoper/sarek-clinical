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
| `gcb-bge-wes.config` | Google Batch profile: Sarek `--step variant_calling --wes`, 4 callers, GATK.GRCh38, Spot. |
| `make_samplesheet.sh` | `families.tsv` → Sarek samplesheet, naming samples `<family>-<role>` for candidate-filtering. |
| `run_bge_wes.sh` | Launch the calling on Batch. |
| `consensus_from_results.sh` | Pull the small per-caller VCFs, run `consensus.sh`, emit `<sample>.consensus.vcf.gz`. |
| `families.example.tsv` | Template family/CRAM table. |

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
