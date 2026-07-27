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

## 5. Acceptance criteria — SET THESE BEFORE RUNNING

### 5.1 Analytical validation (§6)

Measured on the reference sample, **within the reporting region** (the gene panel or capture target
actually used clinically) intersected with the truth set's high-confidence regions.

| Metric | Provisional threshold | Set by lab |
|---|---|---|
| SNV recall | ≥ 0.99 | ☐ |
| SNV precision | ≥ 0.99 | ☐ |
| INDEL recall | ≥ 0.95 | ☐ |
| INDEL precision | ≥ 0.95 | ☐ |

Rationale for the indel thresholds being lower: indel calling is intrinsically harder and every
published benchmark shows a gap to SNVs. Setting one threshold for both would either be
unachievable for indels or meaningless for SNVs.

### 5.2 Per-sample QC gate (§7)

| Metric | Provisional threshold | Set by lab |
|---|---|---|
| Variant count | within [20,000 – 8,000,000] | ☐ |
| Ti/Tv | ≥ 1.8 (WGS ~2.0–2.1; exome ~2.8–3.3) | ☐ |
| het/hom ratio | ≤ 3.0 | ☐ |
| Skewed-AB het fraction | ≤ 0.15 | ☐ |
| Mean depth at called sites | ≥ 20× | ☐ |
| Sex concordance | must match the requisition | ☐ |

**A sample failing any criterion in 5.2 must not have its candidate list reported.**

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
