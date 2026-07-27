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
| **DeepVariant alone (backbone only)** | **0.9984** | 0.9933 | **0.9959** | **6,225** | 25,919 |
| union (clinical default) | 0.9959 | **0.9936** | 0.9948 | 16,043 | 24,796 |
| `NCALLERS>=2` | 0.9961 | 0.9924 | 0.9942 | 15,260 | 29,609 |
| `CONF=HIGH` (≥3 callers) | **0.9982** | 0.9882 | 0.9932 | 7,022 | 45,920 |

### Split by variant class

| tier | class | TP | FP | FN | precision | recall | F1 |
|---|---|---|---|---|---|---|---|
| **DeepVariant alone** | SNV | 3,344,578 | 4,586 | 21,796 | 0.9986 | 0.9935 | **0.9961** |
| **DeepVariant alone** | INDEL | 521,222 | 1,612 | 4,175 | 0.9971 | 0.9921 | **0.9946** |
| union | SNV | 3,345,197 | 13,278 | 21,177 | 0.9960 | 0.9937 | 0.9949 |
| union | INDEL | 521,739 | 2,723 | 3,658 | 0.9951 | 0.9930 | **0.9941** |
| `NCALLERS>=2` | SNV | 3,341,800 | 12,623 | 24,574 | 0.9962 | 0.9927 | 0.9945 |
| `NCALLERS>=2` | INDEL | 520,289 | 2,628 | 5,108 | 0.9953 | 0.9903 | 0.9928 |
| `CONF=HIGH` | SNV | 3,333,022 | 3,942 | 33,352 | 0.9988 | 0.9901 | 0.9944 |
| `CONF=HIGH` | INDEL | 512,589 | 3,080 | 12,808 | 0.9944 | 0.9756 | 0.9849 |

## Exome-restricted (Agilent SureSelect V6 padded capture BED ∩ GIAB high-confidence)

179,987 intersected regions, 97,631 true variants.

| tier | class | TP | FP | FN | precision | recall | **F1** |
|---|---|---|---|---|---|---|---|
| union | SNV | 86,551 | 824 | 1,080 | 0.9906 | 0.9877 | **0.9891** |
| union | INDEL | 10,563 | 66 | 103 | 0.9941 | 0.9903 | **0.9922** |
| union | combined | 97,092 | 890 | 1,181 | 0.9910 | 0.9880 | **0.9895** |
| DeepVariant alone | SNV | 86,511 | 220 | 1,120 | 0.9975 | 0.9872 | **0.9923** |
| DeepVariant alone | INDEL | 10,549 | 30 | 117 | 0.9973 | 0.9890 | **0.9931** |
| DeepVariant alone | combined | 97,060 | 250 | 1,237 | 0.9974 | 0.9874 | **0.9924** |

## Accuracy varies more by REGION than by tier

Same calls, three regions — this is the single most useful table here:

| region | size | true variants | union F1 | DeepVariant-alone F1 |
|---|---|---|---|---|
| **g4e panel** (exons ±20 bp, 1,066 genes) | 5.67 Mb | 4,111 | **0.9977** | **0.9981** |
| genome-wide (GIAB high-confidence) | ~2.5 Gb | 3,890,524 | 0.9948 | 0.9959 |
| **exome capture** (SureSelect V6 padded) | ~100 Mb | 97,631 | **0.9895** | **0.9924** |

**The padded exome capture region is the pipeline's *worst* territory — worse than the genome as a
whole** — while the curated gene panel is its best. A single genome-wide F1 hides a ~0.008 spread.
Plausible reason: a padded capture BED deliberately includes flanking, GC-rich promoter and first-exon
sequence, which is harder to call than average genome; the panel is 1,066 well-characterised,
mostly single-copy disease genes.

**The rescue arm's cost tracks region breadth**: negligible in the panel (0 TP gained, 3 FP added),
11.9 FP per TP in the padded exome (+54 TP, +640 FP), 8.7 FP per TP genome-wide. It is worth keeping
where you report from a tight panel, and worth reconsidering if you ever report across a broad
capture region.

## Restricted to the g4e epilepsy panel — the territory that actually matters

Panel BED built from GENCODE v38 exons (±20 bp for splice regions) of the 1,078 Genes4Epilepsy
symbols; 1,066 located (12 are newer HGNC symbols absent from v38, e.g. `GBA1`, `BMAL1`). 6.20 Mb,
intersected with GIAB high-confidence → **5.67 Mb, 16,671 intervals, 4,111 true variants**.

| tier | SNV TP / FP / FN | INDEL TP / FP / FN | total FP |
|---|---|---|---|
| **DeepVariant alone** | 3,515 / **4** / 6 | 586 / **2** / 4 | **6** |
| union (clinical default) | 3,515 / 5 / 6 | 586 / 4 / 4 | 9 |
| `NCALLERS>=2` | 3,515 / 5 / 6 | 586 / 3 / 4 | 8 |
| `CONF=HIGH` | 3,515 / 4 / 6 | 585 / 4 / 5 | 8 |

**Every tier recovers exactly the same 4,101 true variants.** The rescue arm recovered **zero**
additional true variants inside the panel, and added **3 false positives** (6 → 9). `CONF=HIGH` is
the only tier that loses anything: one indel.

**So inside the clinical territory the tier choice is practically irrelevant** — a 3-variant
difference across 1,066 genes, on calls that the downstream rare/damaging/MANE gates would very
likely drop anyway. The genome-wide 8.7:1 penalty is real but concentrated **outside** the regions
this pipeline reports on.

Two caveats. **Numbers are small** — 6 vs 9 false positives cannot be distinguished statistically
from one sample, and ~2–3 rescued true variants would have been the expected count in 5.67 Mb, so
observing zero is consistent with a small benefit as well as none. And **HG002 is a healthy
reference genome**: it says nothing about how the tiers behave on the pathogenic variants a patient
actually carries.

## What this says about the consensus design

> **CORRECTION (2026-07-27).** An earlier version of this file claimed the union tier had the best
> F1 and that "the rescue arm earns its place". That was wrong. It compared the union only against
> *stricter* tiers and never against the obvious baseline — DeepVariant on its own. Scoring that
> baseline reverses the conclusion, and the corrected finding is below.

**DeepVariant alone outperforms the union consensus**, in both variant classes:

| | SNV F1 | INDEL F1 | FP | FN |
|---|---|---|---|---|
| DeepVariant alone | **0.9961** | **0.9946** | **6,225** | 25,919 |
| Union (+ ≥2-caller rescue) | 0.9949 | 0.9941 | 16,043 | 24,796 |

**The ≥2-caller rescue arm is net-negative on this benchmark.** It recovered **1,123** true
variants that DeepVariant missed, and introduced **9,818** false positives doing so — a ratio of
**8.7 false positives for every true variant recovered**. Precision falls 0.9984 → 0.9959 to buy a
recall gain of 0.0003.

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

**What to do about the rescue arm — measured answer: leave it alone.** The panel-restricted
comparison above shows the genome-wide penalty does not reach the reported territory (0 true
variants gained, 3 false positives added, across 1,066 genes). Changing the default would be
optimising a number nobody reads. The options below remain available if the genome-wide behaviour
ever matters (e.g. a WGS arm reporting outside a panel):

1. **Report `GT_SOURCE != deepvariant` rows as a separate, lower-confidence tier** rather than mixing
   them into the primary list. Costs nothing, keeps the sensitivity, moves the false positives out of
   the main review burden.
2. **Tighten the rescue** — require the rescued variant to pass a quality floor (depth, GQ, allele
   balance) rather than caller agreement alone. The 8.7:1 ratio suggests most rescues are
   low-quality sites where two weaker callers agree *because* they share a failure mode.
3. **Drop it** and use the DeepVariant backbone alone, accepting 1,123 more missed variants
   genome-wide (~0.03% of true variants) for 9,818 fewer false ones.

Caveat before acting: this is measured **genome-wide on one sample**. In a panel/exome context the
absolute numbers are ~40× smaller, and the 1,123 recovered variants are not randomly distributed —
if any fall in clinically relevant genes the trade looks different. Re-measure inside your actual
panel before changing the default.

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
