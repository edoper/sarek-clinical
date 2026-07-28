# SOP-001 — Analytical validation and per-run quality control of the sarek-clinical germline pipeline

| | |
|---|---|
| **Document** | SOP-001 |
| **Version** | 1.0 (DRAFT — unapproved) |
| **Effective date** | *(to be set on approval)* |
| **Author** | *(to be completed)* |
| **Reviewed by** | *(to be completed)* |
| **Approved by** | *(to be completed — must be the laboratory director or delegate)* |
| **Review interval** | 12 months, or on any change to §4 |

> **STATUS: DRAFT TEMPLATE, NOT AN APPROVED PROCEDURE.** The acceptance criteria in §5 are
> *placeholders* carrying provisional values. They are a clinical decision for the laboratory,
> not a technical one, and must be set and signed off before this document governs anything.
> Until then this describes what was done; it does not certify it.

---

## 1. Purpose

To define how the analytical performance of the germline variant-calling pipeline is measured, what
performance is required, and what per-sample quality control every clinical run must pass.

## 2. Scope

**In scope:** single-nucleotide variants (SNVs) and small insertions/deletions (< 50 bp) called from
Illumina short-read WGS or exome data aligned to GRCh38, within regions where the reference truth
set is confident.

**Explicitly OUT of scope — this pipeline must not be relied on for these:**

- structural variants, copy-number variants, repeat expansions
- mosaicism and low-allele-fraction somatic variation
- mitochondrial and pharmacogenomic star-allele calling
- any region outside the reference truth set's high-confidence intervals (notably segmental
  duplications and many repeat classes)
- clinical interpretation — the downstream `candidate-filtering` step is a **triage aid**; every
  candidate requires review and sign-out by a qualified professional

## 3. Definitions

| Term | Meaning |
|---|---|
| Truth set | An independently established set of genotypes for a reference sample (here: GIAB) |
| High-confidence regions | The BED accompanying the truth set, where its genotypes are reliable |
| Recall (sensitivity) | TP / (TP + FN) — of true variants, the fraction found |
| Precision | TP / (TP + FP) — of variants called, the fraction real |
| F1 | Harmonic mean of precision and recall |
| Ti/Tv | Transition/transversion ratio; a call-set purity indicator |
| QC gate | A per-sample pass/fail decision made before any candidate list is read |

## 4. Materials

| Item | Value | Where recorded |
|---|---|---|
| Pipeline | `sarek-clinical`, pinned commit | `RUN-PROVENANCE.md` per run |
| Workflow | nf-core/sarek, pinned `-r` version | `gcb.config`, run log |
| Callers | DeepVariant, Strelka2, FreeBayes, GATK HaplotypeCaller | `gcb.config` |
| Consensus | `consensus.sh` (DeepVariant backbone + ≥2-caller rescue) | repo |
| Reference | GRCh38 (`Homo_sapiens_assembly38`, ALT/decoy included) | `gcb.config` |
| Reference sample | GIAB HG002 / NA24385 | `validation/README.md` |
| Truth set | GIAB v4.2.1, GRCh38 | `validation/README.md` |
| Comparator | RTG Tools `vcfeval`, haplotype-aware | `validation/benchmark_giab.sh` |

**Any change to a row in this table invalidates the current validation and requires §6 to be
repeated.** This is the change-control rule.

## 5. Acceptance criteria

### 5.0 Define the reportable range FIRST

A sensitivity threshold is meaningless without the region it applies to. This pipeline's measured
accuracy varies **more by region than by any filtering choice**: g4e panel F1 0.9977, genome-wide
0.9948, padded exome capture 0.9895. State the reportable range, then set thresholds for it.

**Reportable range for this laboratory:** ☐ ______________________________
*(recommended: the g4e panel — MANE exons ±20 bp of the panel genes — since that is what is
reported. Anything outside it is explicitly not covered by these criteria.)*

### 5.1 What the statistics allow you to claim

Thresholds cannot be tighter than the confidence interval of the measurement that verifies them.
95% Wilson intervals on recall from the 2026-07-27 HG002 run:

| region | class | true variants | recall | 95% CI | CI width |
|---|---|---|---|---|---|
| g4e panel | SNV | 3,521 | 0.9983 | [0.9963, 0.9992] | 0.29 pp |
| g4e panel | **INDEL** | **590** | 0.9932 | **[0.9827, 0.9974]** | **1.47 pp** |
| exome capture | INDEL | 10,666 | 0.9903 | [0.9883, 0.9920] | 0.37 pp |
| genome-wide | INDEL | 525,397 | 0.9930 | [0.9928, 0.9933] | 0.04 pp |

**The panel contains only 590 true indels.** One panel run therefore cannot distinguish 99.0% from
99.6% indel sensitivity — the interval spans both. **Setting a 99% indel threshold verified only by
a panel run would be aspirational, not evidenced.** Two honest ways forward:

1. **Verify on the wide region, confirm on the narrow one.** Establish the numeric threshold using
   the genome-wide or exome-wide measurement (hundreds of thousands of indels → CI ±0.04 pp), then
   use the panel run to confirm no region-specific degradation. This is the recommended route.
2. **Add samples.** GIAB publishes HG001–HG007 with the same truth-set machinery. Pooling indels
   across samples narrows the panel interval and simultaneously provides the reproducibility
   evidence §6 requires.

### 5.2 Analytical validation thresholds — two-tier

Two limits, deliberately. A single threshold you fail occasionally creates pressure to rationalise
failures; a floor plus a target keeps that honest.

**Thresholds must be set PER ASSAY.** The exome validation (2026-07-28) showed the same pipeline
achieves very different indel performance on WGS versus a real capture library. A single set of
numbers covering both would be wrong for one of them.

**WGS arm** (measured on HG002 30× PCR-free, panel region):

| Metric | **Floor** | **Target** | Measured |
|---|---|---|---|
| SNV sensitivity | ☐ *(suggest 99.0%)* | ☐ *(suggest 99.5%)* | 99.83% |
| SNV precision | ☐ *(suggest 99.0%)* | ☐ *(suggest 99.5%)* | 99.86% |
| INDEL sensitivity | ☐ *(suggest 98.0%)* | ☐ *(suggest 99.0%)* | 99.32% |
| INDEL precision | ☐ *(suggest 98.0%)* | ☐ *(suggest 99.0%)* | 99.35% |

**EXOME / BGE arm** (measured on a real HG002 capture library, panel ∩ calling intervals):

| Metric | **Floor** | **Target** | Measured | note |
|---|---|---|---|---|
| SNV sensitivity | ☐ *(suggest 98.5%)* | ☐ *(suggest 99.0%)* | 99.39% | CI ±0.4 pp |
| SNV precision | ☐ *(suggest 98.5%)* | ☐ *(suggest 99.0%)* | 99.72% | |
| INDEL sensitivity | ☐ *(suggest 88.0%)* | ☐ *(suggest 93.0%)* | 92.40% | **CI [87.4, 95.5] — 8 pp wide** |
| INDEL precision | ☐ *(suggest 88.0%)* | ☐ *(suggest 93.0%)* | 91.57% | |

The exome indel floor is suggested near the **lower confidence bound**, not the point estimate,
because 171 indels cannot support anything tighter, and because the reference exome used a
different capture kit (V5 vs V6) which depresses recall for reasons unrelated to the pipeline.
**Re-measure with a reference exome captured on your own kit before tightening these.**

Notes for whoever signs this off:

- **The suggested floors are what the current evidence supports**, not the best imaginable numbers.
  The panel indel CI lower bound is 98.27%, so a 98% floor is verified; a 99% floor is not — yet.
- **SNV thresholds must be region-aware.** Genome-wide SNV sensitivity is 99.37%, *below* a 99.5%
  target. If the reportable range is the panel, 99.5% is appropriate; if it is genome-wide, it is
  not. Do not copy one into the other.
- **Indels get lower thresholds than SNVs deliberately.** Every published benchmark shows the gap.
  For context, a 2025 software comparison put GATK-based tools at indel F1 0.89–0.93 and Illumina
  DRAGEN at 0.9699, against SNV F1 0.98–0.9969; this pipeline measured indel F1 0.9941 genome-wide.
  A single threshold for both classes is either unachievable for indels or vacuous for SNVs.
- Published clinical WES/WGS validations commonly claim **>99% sensitivity and specificity** for
  SNVs and small indels, with the strongest reporting >99.7%/99.8% in high-mapping-quality regions.
  The suggested targets sit inside that band.

### 5.3 Per-sample QC gate (§7)

| Metric | Suggested | Set by lab |
|---|---|---|
| Mean depth over the reportable range | ≥ 50× (exome/panel) | ☐ |
| Fraction of reportable range at ≥ 20× | ≥ 99% | ☐ |
| Variant count | within [20,000 – 8,000,000] | ☐ |
| Ti/Tv | ≥ 1.8 (WGS ~2.0–2.1; exome ~2.8–3.3) | ☐ |
| het/hom ratio | ≤ 3.0 | ☐ |
| Skewed-AB het fraction (contamination proxy) | ≤ 0.15 | ☐ |
| Sex concordance | must match the requisition | ☐ |

Coverage figures follow common clinical practice (mean >35× with 99.5% of the reportable range
above 20×); 50× is suggested here because capture is uneven and the margin costs little.

**A sample failing any criterion must not have its candidate list reported.**

## 6. Procedure — analytical validation

1. Confirm the environment: `source env.sh`, then `./test/test_consensus.sh` — must print
   `ALL TESTS PASSED`.
2. Obtain the reference sample data and truth set (`validation/README.md` §0.1–0.2).
3. Run `validation/run_giab.sh`. It executes the pipeline with the **same configuration used
   clinically** — no tuning specific to the reference sample is permitted, as that would invalidate
   the result.
4. On completion, record in `RESULTS-*.md`: pipeline commit, truth-set version, region definition,
   and precision/recall/F1 split by variant class.
5. Compare against §5.1. Record PASS or FAIL **per criterion**, not as an overall impression.
6. Any deviation from this procedure is recorded in §9 with its justification.

**Repeat this procedure when:** any row of §4 changes; at the defined review interval; or after any
change to `consensus.sh`, the caller set, or the reference.

**Frequency and replication.** A single run of a single reference sample is a point estimate and is
**not sufficient** for accreditation. A complete validation additionally requires repeatability
(the same sample processed repeatedly), reproducibility (different days/operators/batches), and
ideally a second reference sample. Record here what was actually done, and do not describe a
single-run validation as complete.

## 7. Procedure — per-run quality control

For **every** clinical sample, before any candidate list is reviewed:

```bash
./qc_gate.sh <sample>.consensus.vcf.gz --sex <M|F> --targets <panel.bed> --json <sample>.qc.json
```

- **exit 0 (PASS)** — proceed to interpretation.
- **exit 2 (WARN)** — a qualified person reviews the flagged metric and records a decision before
  proceeding.
- **exit 1 (FAIL)** — **stop.** Do not report. Investigate: contamination, sample swap, capture or
  library failure. Repeat from the wet lab if the cause is not resolved.

Retain the `.qc.json` with the run record. A sex mismatch or an elevated skewed-allele-balance
fraction is treated as a **suspected sample identity or contamination event** and escalated per
laboratory policy, not silently overridden.

## 8. Records

For each validation and each clinical run, retain: pipeline commit hash, `RUN-PROVENANCE.md`,
the QC gate JSON, the vcfeval summaries (validation only), and the reviewer's name and date.
Retention period per laboratory policy and applicable law.

## 9. Deviations

| Date | Deviation | Justification | Approved by |
|---|---|---|---|
| | | | |

## 10. Known limitations to carry into every report

- Performance is characterised only inside the truth set's high-confidence regions.
- Out-of-scope variant classes (§2) are not detected and their absence is not evidence of absence.
- The pipeline produces candidates; it does not make diagnoses.
- Cloud processing occurs outside the country of origin; confirm this is consistent with applicable
  data-protection law and institutional policy before processing patient data.

## 11. Revision history

| Version | Date | Change | Author |
|---|---|---|---|
| 1.0 | 2026-07-27 | Initial draft, unapproved | — |
