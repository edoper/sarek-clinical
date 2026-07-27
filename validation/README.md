# Pipeline validation against GIAB

Accuracy validation of this pipeline against a **Genome in a Bottle** reference sample, using
**HG002 / NA24385** (the Ashkenazim son) and the **GIAB v4.2.1** benchmark truth set on GRCh38.

This is the step README section 6 calls for before results are used for patient care: run a sample
whose true genotypes are known, and report how close the pipeline gets.

## Why this design

- **The pipeline is run exactly as it is run clinically** — same `gcb.config`, same four callers,
  same `--skip_tools baserecalibrator,vcftools,multiqc`. A validation of a specially-tuned variant
  of the pipeline would not describe the pipeline you actually use.
- **Input is public and free to read.** HG002 30× PCR-free NovaSeq FASTQ lives in Google's
  `gs://brain-genomics-public/` bucket, which is `allUsers`-readable, so Batch VMs read it directly
  from the samplesheet — no staging copy, no duplicate storage.
- **Several confidence tiers are scored, not just one.** `consensus.sh` deliberately defers
  strictness downstream; this measures what that choice costs. The union (clinical default),
  `NCALLERS>=2` and `CONF=HIGH` are each compared to truth, so the precision/recall trade is a
  measurement rather than an assumption.
- **An exome-restricted row is reported too.** The clinical arms are exome/BGE, so the same calls
  are re-scored inside the capture BED ∩ GIAB high-confidence regions — closer to production than a
  genome-wide number.
- **Comparison uses RTG `vcfeval`**, which is Java and needs no Docker, and does proper
  haplotype-aware matching (an indel written two ways still matches).

## Files

| File | Purpose |
|---|---|
| `run_giab.sh` | End-to-end: sarek WGS from FASTQ → `consensus.sh` → benchmark |
| `benchmark_giab.sh` | RTG `vcfeval` of each confidence tier vs GIAB, genome-wide and exome-restricted |
| `health_check.sh` | One status snapshot (progress, VMs, reconstructed spend) + **hard budget guard** |
| `RESULTS-HG002.md` | The measured numbers from the most recent run |

## Running it

```bash
source ../env.sh
cd validation

# 1. one-time: RTG Tools + the GIAB truth set + an SDF of the reference
curl -sL -o rtg.zip https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-linux-x64.zip
unzip -q rtg.zip
./rtg-tools-3.12.1/rtg format -o GRCh38.sdf ../refs/Homo_sapiens_assembly38.fasta

G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
mkdir -p truth && for f in HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz{,.tbi} \
                           HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed; do
  curl -sL -o truth/$f $G/$f
done

# 2. the run (several hours on Spot; watch STATUS.md or run health_check.sh)
./run_giab.sh
```

`run_giab.sh` is resumable — re-running it picks up where Nextflow left off.

## Cost and footprint

Measured for one 30× WGS genome on Spot in `us-central1` (see README section 6 for the
per-genome cost model this follows):

| | |
|---|---|
| Compute | ~$5–7 (single WGS genome, four callers, Spot) |
| Input storage | $0 — the public FASTQ is read in place, never copied |
| Scratch (`work/`) | ~200 GB while running; **delete it afterwards** |
| Local disk | ~5 GB (truth set, RTG SDF, consensus VCF) |
| Wall clock | several hours |

`health_check.sh` reconstructs accrued spend from the Batch job records (control-plane listing
only — no compute, no egress) and **cancels the run** if it crosses `KILL_AT`, so an unattended
validation cannot walk past its budget.

## Interpreting the result

- **Recall** = of the true variants in GIAB's high-confidence regions, what fraction did we call.
- **Precision** = of what we called there, what fraction is real.
- **F1** = their harmonic mean; the single number usually quoted.
- Metrics are reported separately for **SNVs** and **INDELs** — indel accuracy is always the lower
  of the two and is the more discriminating number.

Scores apply **only inside GIAB's high-confidence regions**, which deliberately exclude the hardest
parts of the genome (segmental duplications, many repeats). Real-world sensitivity across the whole
genome is lower than these numbers, and this validation says nothing about structural variants,
copy number, or repeat expansions.
