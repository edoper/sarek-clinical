# GIAB validation results — HG002 EXOME (the assay actually delivered)

**Run date:** 2026-07-28 · **Truth set:** GIAB v4.2.1, GRCh38 · **Input:** GIAB HG002 exome,
Oslo University Hospital, **Agilent SureSelect V5**, Illumina, re-aligned to GRCh38 by the pipeline
from FASTQ (`--step mapping --wes`) · **Calling intervals:** SureSelect **V6** padded (this
deployment's clinical BED) · **Comparator:** RTG `vcfeval` · **Cost:** $0.78 · **QC gate:** PASS

## Why this run exists

The WGS validation reports an "exome-restricted" row — WGS reads filtered to a capture BED. That is
**not** an exome: it has none of the capture bias, duplicate structure or uneven target coverage of
a real hybrid-capture library. This run tests the real thing, and the difference turned out to be
large enough to matter clinically.

## Headline — the WGS proxy was optimistic, badly so for indels

Same panel genes, same truth set, same pipeline:

| | SNV F1 | **INDEL F1** |
|---|---|---|
| WGS reads restricted to the panel | 0.9984 | **0.9934** |
| **REAL exome, same panel region** | 0.9955 | **0.9198** |

**SNV performance holds up. Indel performance does not** — F1 falls from 0.9934 to 0.9198 in the
same genes. If you had relied on the WGS number to describe your exome and BGE arms, you would have
overstated indel sensitivity by a wide margin.

## Numbers

**Fair panel region** = g4e panel ∩ V6 calling intervals ∩ GIAB high-confidence = **3.50 Mb**,
1,792 true SNVs and 171 true indels.

| tier | class | TP | FP | FN | precision | recall | **F1** |
|---|---|---|---|---|---|---|---|
| union (default) | SNV | 1,781 | 5 | 11 | 0.9972 | 0.9939 | **0.9955** |
| union (default) | INDEL | 158 | 15 | 13 | 0.9157 | 0.9240 | **0.9198** |
| DeepVariant alone | SNV | 1,781 | 3 | 11 | 0.9983 | 0.9939 | **0.9961** |
| DeepVariant alone | INDEL | 155 | 7 | 16 | 0.9581 | 0.9064 | **0.9315** |

**Whole V6 capture region** ∩ high-confidence (97,631 true variants):

| tier | class | precision | recall | F1 |
|---|---|---|---|---|
| union | SNV | 0.9839 | 0.9461 | 0.9646 |
| union | INDEL | 0.9242 | 0.8584 | 0.8901 |

**DeepVariant alone again beats the union**, consistent with the WGS finding — and by more on indels
here (0.9315 vs 0.9198), because the rescue arm adds 8 indel false positives to gain 3 true ones.

## Read these caveats before using any number above

1. **Kit mismatch — every recall figure is a LOWER BOUND.** The reference exome was captured with
   **V5**; calling used the **V6** padded BED. Regions present in V6 but not V5 have no reads at
   all, so they generate guaranteed false negatives that are the kit's doing, not the pipeline's.
   The only way to remove this confound is a reference exome captured with your own kit.
2. **The indel estimate is imprecise.** It rests on 171 true indels: recall 0.9240, 95% CI
   **[0.8743, 0.9550]** — 8.1 pp wide. The *direction* (exome indels are much worse than WGS
   indels) is solid; the *magnitude* is not. SNV recall 0.9939 has a CI of [0.9890, 0.9966], 0.8 pp.
3. **The reference exome is from 2015.** Older chemistry and read length than anything you run now.
   Some of the gap is the data, not the pipeline.
4. **A naive panel evaluation is meaningless here and was discarded.** Only **61%** of the g4e panel
   BED lies inside the V6 calling intervals; the other 39% was never called, which produced an
   apparent panel recall of 0.4717. That number is a BED artefact, not performance. All panel
   figures above are restricted to the region where calling was actually attempted.

## What to do about it

- **Do not quote the WGS exome-restricted row as evidence for your exome or BGE arms.** Quote this
  file. Update any method text that cites a single genome-wide F1.
- **Treat indels as the limiting factor in exome reporting.** For an epilepsy panel where
  frameshifts matter, a plausible indel sensitivity in the low-to-mid 90s is a material clinical
  limitation and belongs in the report's limitations section.
- **Re-measure with your own kit.** A reference exome captured with SureSelect V6 (or the BGE
  protocol) would remove caveat 1 and is the single highest-value next validation.
- **Consider orthogonal confirmation for reported indels**, which is what ACMG guidance already
  advises for indels specifically until a laboratory has extensive validation showing it is safe to
  drop.

## Reproducing

`validation/run_giab_exome.sh` — streams the reference exome BAM from GIAB, converts to FASTQ,
runs the pipeline in `--wes` mode, applies the QC gate, and benchmarks. ~2 h, **$0.78**.
