# GIAB validation results — HG002 (NA24385)

**Run date:** 2026-07-27 · **Truth set:** GIAB v4.2.1, GRCh38 · **Input:** HG002 30× PCR-free
NovaSeq WGS (`gs://brain-genomics-public/`) · **Comparison:** RTG `vcfeval` 3.12.1, haplotype-aware,
restricted to GIAB's high-confidence regions.

The pipeline was run **exactly as it is run clinically** — `gcb.config`, all four callers,
`--skip_tools baserecalibrator,vcftools,multiqc`, Spot VMs — so these numbers describe the real
pipeline and not a tuned variant of it. Note this means CNNScoreVariants **did** run, and the
consensus consumed `HG002.haplotypecaller.filtered.vcf.gz`, exactly as production does.

## Headline

**Genome-wide F1 = 0.9948** (SNV 0.9949, INDEL 0.9941) for the clinical default (the full union
consensus), inside GIAB high-confidence regions.

## Consensus composition

| | variants |
|---|---|
| DeepVariant backbone | 4,927,041 |
| Rescued via Strelka2 | 46,641 |
| Rescued via HaplotypeCaller | 70,434 |
| **Union total (clinical default)** | **5,044,116** |
| `NCALLERS>=2` subset | 4,927,092 |
| `CONF="HIGH"` (≥3 callers) subset | 4,671,576 |

## Accuracy by confidence tier (genome-wide, GIAB high-confidence regions)

| tier | precision | recall | F1 | FP | FN |
|---|---|---|---|---|---|
| **union (clinical default)** | 0.9959 | **0.9936** | **0.9948** | 16,043 | 24,796 |
| `NCALLERS>=2` | 0.9961 | 0.9924 | 0.9942 | 15,260 | 29,609 |
| `CONF=HIGH` (≥3 callers) | **0.9982** | 0.9882 | 0.9932 | 7,022 | 45,920 |

### Split by variant class

| tier | class | TP | FP | FN | precision | recall | F1 |
|---|---|---|---|---|---|---|---|
| union | SNV | 3,345,197 | 13,278 | 21,177 | 0.9960 | 0.9937 | **0.9949** |
| union | INDEL | 521,739 | 2,723 | 3,658 | 0.9951 | 0.9930 | **0.9941** |
| `NCALLERS>=2` | SNV | 3,341,800 | 12,623 | 24,574 | 0.9962 | 0.9927 | 0.9945 |
| `NCALLERS>=2` | INDEL | 520,289 | 2,628 | 5,108 | 0.9953 | 0.9903 | 0.9928 |
| `CONF=HIGH` | SNV | 3,333,022 | 3,942 | 33,352 | 0.9988 | 0.9901 | 0.9944 |
| `CONF=HIGH` | INDEL | 512,589 | 3,080 | 12,808 | 0.9944 | 0.9756 | 0.9849 |

## Exome-restricted (Agilent SureSelect V6 capture BED ∩ GIAB high-confidence)

179,987 intersected regions.

| class | TP | FP | FN | precision | recall | F1 |
|---|---|---|---|---|---|---|
| SNV | 86,551 | 824 | 1,080 | 0.9906 | 0.9877 | 0.9891 |
| INDEL | 10,563 | 66 | 103 | 0.9941 | 0.9903 | 0.9922 |
| combined | 97,092 | 890 | 1,181 | 0.9910 | 0.9880 | 0.9895 |

## What this says about the consensus design

**The union — the pipeline's clinical default — has the best F1 of the three tiers, in both variant
classes.** That is the design decision validated: `consensus.sh` keeps every DeepVariant call plus
≥2-caller rescues and defers strictness downstream, and the measurement says that deferral is not
costing accuracy.

**Tightening buys precision at a disproportionate cost in recall.** Going from the union to
`CONF=HIGH` removes 9,021 false positives (16,043 → 7,022, −56%) but adds 21,124 false negatives
(24,796 → 45,920, +85%). For germline diagnosis a missed variant is worse than an extra one to
review, so the default is the right operating point — but the tags are there if a particular use
case wants the stricter tier.

**Indels are where the tiers diverge.** SNV F1 barely moves across tiers (0.9949 → 0.9944), but
indel F1 falls from 0.9941 to 0.9849, driven by recall collapsing from 0.9930 to 0.9756 — `CONF=HIGH`
misses 3.5× as many indels. Cross-caller agreement is weakest exactly where callers legitimately
disagree on indel representation, so requiring three callers discards real indels. **Do not use
`CONF=HIGH` for indel-sensitive questions.**

**The rescue arm earns its place.** 117,075 variants (2.3% of the consensus) came from the ≥2-caller
rescue rather than DeepVariant. The union's recall exceeds the `NCALLERS>=2` tier's (0.9936 vs
0.9924) while its precision is essentially identical, so those rescues are net-positive.

## Limitations — read before citing

- Scores apply **only inside GIAB's high-confidence regions**, which exclude segmental duplications
  and many repeats. Real-world genome-wide sensitivity is lower.
- **Single sample, single run.** One reference genome is a point estimate, not a distribution across
  operators, batches or runs.
- **SNVs and small indels only.** Nothing here validates structural variants, copy number, repeat
  expansions, or mosaicism.
- The input is **30× PCR-free WGS**. The exome row is the closest available proxy for the BGE/exome
  arms, but it is the same WGS reads restricted to a capture BED — not real capture data, with its
  own coverage profile and PCR duplicates.
- **Triage, not release.** This measures the variant-calling pipeline. It says nothing about the
  downstream `candidate-filtering` interpretation step.

## Reproducing

`validation/run_giab.sh` — end to end. Measured cost for this run: **$4.68** of Spot compute
(178 Batch jobs, 19 preempted and retried), 11h23m wall clock, ~$0 storage since the public FASTQ is
read in place. See `validation/README.md`.
