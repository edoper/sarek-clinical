# Validation and quality control

Everything needed to (a) measure this pipeline's accuracy against a reference sample with known
genotypes, and (b) gate individual clinical samples on quality before anyone reads their results.

**If you are a new user or an agent picking this up cold, read §0 then follow §2 or §3 verbatim.**
Every command is copy-pasteable and every input is public.

| File | What it is |
|---|---|
| `SOP.md` | Draft standard operating procedure — scope, acceptance criteria, change control |
| `run_giab.sh` | **WGS** validation, end to end |
| `run_giab_exome.sh` | **Exome** validation, end to end (matches the clinical exome arm) |
| `benchmark_giab.sh` | RTG `vcfeval` scoring, per confidence tier and per region |
| `split_metrics.sh` | Per-class (SNV/INDEL) precision, recall, F1 from vcfeval output |
| `health_check.sh` | Live progress + reconstructed spend, with a hard budget guard |
| `RESULTS-HG002.md` | Measured WGS results |
| `RESULTS-HG002-EXOME.md` | Measured exome results |
| `../qc_gate.sh` | **Per-sample QC gate** — run this on every clinical sample |

---

## 0. Concepts, in one page

**Why validate at all.** A variant caller that silently misses 5% of true variants produces a
perfectly plausible result file. Nothing downstream reveals the gap. The only way to know is to run
a sample whose true genotypes are already established and count the differences.

**What we compare against.** [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
publishes reference samples with consensus genotypes derived from many technologies. We use
**HG002/NA24385** and the **v4.2.1** truth set. GIAB also ships a BED of *high-confidence regions*;
scores are only meaningful inside it, because outside it the "truth" is itself uncertain.

**The three numbers.**

- **Recall** (sensitivity) = TP / (TP + FN) — of the variants that are really there, how many did we find. Misses are missed diagnoses.
- **Precision** = TP / (TP + FP) — of the variants we called, how many are real. False positives are wasted curation.
- **F1** = harmonic mean of the two. One number, but it hides which of the two is weak — always read the pair.

Report SNVs and indels **separately**. Indels are harder and always score lower; a combined number
is dominated by SNVs and flatters the pipeline.

**Why a haplotype-aware comparator.** The same indel can be written correctly in several ways
(`AT→A` at one position vs `TA→T` at the next). Naive position matching scores those as both a false
positive and a false negative. RTG `vcfeval` compares *haplotypes*, so equivalent representations
match. Never benchmark variants with `bedtools intersect`.

**Validation ≠ QC.** Validation is done once per pipeline version and asks *"is this method good
enough?"*. QC is done on every sample and asks *"did this particular sample work?"*. You need both:
a validated pipeline still produces garbage from a contaminated library.

---

## 1. One-time setup

Needs ~10 GB local disk plus whatever the reference genome already occupies.

```bash
source ../env.sh            # toolchain + site settings
cd validation

# RTG Tools — Java, no Docker needed
curl -sL -o rtg.zip https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-linux-x64.zip
unzip -q rtg.zip

# A sequence index of the reference (vcfeval needs this format, ~1.2 GB, a few minutes)
./rtg-tools-3.12.1/rtg format -o GRCh38.sdf ../refs/Homo_sapiens_assembly38.fasta

# The GIAB truth set (~160 MB)
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
mkdir -p truth
for f in HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
         HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
         HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed; do
  curl -sL -o "truth/$f" "$G/$f"
done
bcftools index -n truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz   # expect ~4,048,342
```

---

## 2. WGS validation

```bash
./run_giab.sh        # several hours, ~$5 of Spot compute
```

Input is HG002 30× PCR-free FASTQ from `gs://brain-genomics-public/`. That bucket is
`allUsers`-readable, so Batch VMs read it **in place** — no staging copy, no duplicate storage. The
samplesheet points straight at the public URIs.

Watch it without spending anything:

```bash
watch -n 60 ./health_check.sh    # progress, VMs, reconstructed spend
```

`health_check.sh` also enforces a **hard budget guard**: if reconstructed spend crosses `KILL_AT`
(default $12) it cancels the driver and every live Batch job. An unattended validation cannot walk
past its ceiling.

Results → `RESULTS-HG002.md`.

## 3. Exome validation — do this one too

**The WGS validation does not certify the exome arm.** `run_giab.sh` reports an exome-restricted
row, but that is WGS reads filtered to a capture BED: it has none of the capture bias, duplicate
structure or uneven target coverage of a real hybrid-capture library, and those are precisely what
make exome calling harder. If your clinical work is exome or BGE, this is the validation that
matters.

```bash
./run_giab_exome.sh   # ~1-2 h, well under $1 of Spot compute
```

It streams GIAB's HG002 exome BAM (Oslo University Hospital, Agilent SureSelect V5, GRCh37),
converts it back to paired FASTQ, and lets the pipeline re-align to GRCh38 itself via
`--step mapping --wes`. Nothing from the original GRCh37 alignment is carried forward.

> **Kit caveat, read before interpreting.** The reference exome is captured with SureSelect **V5**;
> this deployment's clinical BED is **V6**. Regions in V6 but not V5 have no reads, and would score
> as false negatives caused by the *kit mismatch*, not by the pipeline. Judge the pipeline on the
> gene-panel numbers and treat the whole-capture numbers as a lower bound. If you capture with a
> different kit, substitute a reference exome captured with yours where one exists.

Results → `RESULTS-HG002-EXOME.md`.

## 4. Per-sample QC gate — run on EVERY clinical sample

Validation says the *method* works. The gate says *this sample* worked.

```bash
../qc_gate.sh <sample>.consensus.vcf.gz --sex M --targets panel.bed --json <sample>.qc.json
```

Exit codes: **0 = PASS** (proceed) · **2 = WARN** (a human reviews and records a decision) ·
**1 = FAIL** (**do not report** — investigate contamination, sample swap, or capture failure).

What it checks and why:

| Check | Catches |
|---|---|
| Variant count within bounds | Failed capture / collapsed library; artefact-laden call set |
| Ti/Tv ratio | A call set dominated by noise (random calls give ~0.5; real ~2.0 WGS, ~3.0 exome) |
| het/hom ratio | Contamination inflates it; wrong reference or consanguinity deflates it |
| Skewed allele-balance het fraction | **Contamination** — foreign alleles appear as low-AB hets |
| Mean depth at called sites | Under-sequenced sample |
| chrX het rate + chrY calls vs expected sex | **Sample swap** |

Thresholds are conservative defaults, overridable by environment variable. **A clinical laboratory
must set its own from its own data and record them in `SOP.md` §5.2** — that is a clinical decision,
not a technical one.

Verify the gate has teeth before trusting it (a gate that never fails is worse than none):

```bash
../qc_gate.sh good.vcf.gz --sex M      # expect PASS
../qc_gate.sh good.vcf.gz --sex F      # expect FAIL: SEX MISMATCH
```

---

## 5. Cost, time and footprint

| | WGS validation | Exome validation |
|---|---|---|
| Compute (Spot, us-central1) | ~$5 | < $1 |
| Wall clock | ~7 h calling + a QC tail | 1–2 h |
| Bucket scratch | ~135 GB — **delete `work/` afterwards** | ~15 GB |
| Local disk | ~5 GB | ~25 GB (FASTQ) |
| Input storage | $0 (public data read in place / streamed) | $0 |

## 6. Interpreting and reporting honestly

- Scores apply **only inside GIAB high-confidence regions**. Real-world sensitivity across the whole
  genome is lower — those regions deliberately exclude segmental duplications and many repeats.
- **One sample, one run is a point estimate**, not a validation dossier. Accreditation normally
  expects repeatability, reproducibility across days/operators, and more than one reference sample.
- **SNVs and small indels only.** Nothing here covers structural variants, CNVs, repeat expansions
  or mosaicism — and their absence from a report is not evidence of their absence in the patient.
- This validates **variant calling**. It says nothing about the downstream `candidate-filtering`
  interpretation step, which needs its own known-positive benchmark.
- Accuracy varies more by **region** than by filtering tier — see `RESULTS-HG002.md`. Quote the
  number for the region you actually report from, not the genome-wide one.
