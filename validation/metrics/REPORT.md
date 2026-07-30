# GIAB validation of the sarek-clinical pipeline — run report

**Date:** 2026-07-27 · **Sample:** HG002 / NA24385 (GIAB Ashkenazim son) · **Truth:** GIAB v4.2.1, GRCh38
**Authorised budget:** $15 · **Actual cost: $5.49 total** (WGS $4.71 + exome $0.78)

---

## 0. The finding that matters most

**The exome arm is materially worse than the WGS arm, and only for indels.** Running GIAB HG002
through the *real* exome path — not WGS reads restricted to a capture BED — gives, in the same
panel genes:

| | SNV F1 | **INDEL F1** |
|---|---|---|
| WGS reads restricted to the panel | 0.9984 | **0.9934** |
| **REAL exome, same panel region** | 0.9955 | **0.9198** |

Your clinical arms are exome and BGE. **Do not describe them with the WGS number** — it overstates
indel sensitivity by roughly 0.07 F1. For an epilepsy panel, where frameshifts matter, an indel
sensitivity plausibly in the low-to-mid 90s is a real clinical limitation that belongs in report
limitations.

Bounded honestly: the reference exome used SureSelect **V5** while calling used your **V6** BED, so
V6-only regions have no reads and every recall figure is a **lower bound**; and the estimate rests
on 171 indels (95% CI [0.8743, 0.9550], 8 pp wide). The direction is solid, the magnitude is not.
Full detail and caveats: `validation/RESULTS-HG002-EXOME.md`.

## 1. Result (WGS arm)

**The pipeline passes.** Genome-wide, inside GIAB's high-confidence regions, the default union
consensus achieves:

| | precision | recall | **F1** |
|---|---|---|---|
| **SNV** | 0.9960 | 0.9937 | **0.9949** |
| **INDEL** | 0.9951 | 0.9930 | **0.9941** |
| combined | 0.9959 | 0.9936 | **0.9948** |

Good, but **not best-in-class** — see §2, where DeepVariant on its own scores higher. The full per-tier table,
the exome-restricted numbers, and the limitations are in
`sarek-clinical/validation/RESULTS-HG002.md` (now on GitHub).

## 2. The finding that matters for how you run it

Four variants of the *same* call set were scored, because `consensus.sh` deliberately defers
strictness downstream. Among the **tiers of the union**, the default is the best operating point —
but the union itself loses to its own backbone (see the correction below):

| tier | precision | recall | F1 | FP | FN |
|---|---|---|---|---|---|
| **union (default)** | 0.9959 | **0.9936** | **0.9948** | 16,043 | 24,796 |
| `NCALLERS>=2` | 0.9961 | 0.9924 | 0.9942 | 15,260 | 29,609 |
| `CONF=HIGH` | **0.9982** | 0.9882 | 0.9932 | 7,022 | 45,920 |

Tightening to `CONF=HIGH` removes 9,021 false positives but adds **21,124 false negatives**. For
germline diagnosis that is a bad trade — a missed variant costs more than one extra to review.

**Indels are where this bites.** Indel recall falls 0.9930 → 0.9756 at `CONF=HIGH`: 3.5× more missed
indels. Requiring three callers to agree discards real indels, because callers legitimately disagree
on indel representation. **Do not use `CONF=HIGH` for indel-sensitive questions.**

**CORRECTION — the rescue arm is NOT validated.** My first pass compared the union only against
*stricter* tiers and concluded the rescue "earns its place". I then scored the obvious missing
baseline, DeepVariant on its own, and it reverses the conclusion:

| | SNV F1 | INDEL F1 | FP | FN |
|---|---|---|---|---|
| **DeepVariant alone** | **0.9961** | **0.9946** | **6,225** | 25,919 |
| Union (+ ≥2-caller rescue) | 0.9949 | 0.9941 | 16,043 | 24,796 |

The rescue arm recovered **1,123** true variants and introduced **9,818** false positives —
**8.7 FP per TP**. DeepVariant alone is the better caller on this benchmark. See
`RESULTS-HG002.md` for what to do about it (report rescues as a separate tier, tighten them with a
quality floor, or drop them) and why you should re-measure inside your panel before changing
anything.

## 3. What happened overnight — why there were no results when you woke up

The run **did not fail**. It completed successfully, but ~4 hours later than it should have.

- The analytical work — alignment plus all four callers — finished at **09:50 UTC (~7.4h)**.
- It then spent **3 h 55 m on a single `CNNSCOREVARIANTS` task** (2 cores), running after every
  other stage had finished. That is **34% of the whole run's wall clock for 2.7% of its CPU cost**.
- It runs because `haplotypecaller_filter` is absent from the usual `--skip_tools` list.

> **CORRECTION.** I first told you the stall was `BCFTOOLS_STATS`, inferred from the last line of
> the log. The execution trace disproves that: `BCFTOOLS_STATS` ran four times, 6.6 s / 6.9 s /
> 23.1 s / 13.6 s. `CNNSCOREVARIANTS` was the culprit. `CLAUDE.md` and `validation/EFFICIENCY.md`
> now carry the corrected version.

Cost impact was small — the tail was one small VM, ~$0.08 over four hours. The impact was **time**:
my 4–8h estimate was for the analytical work and did not account for the QC tail.

**Fix for next time:** treat the serial tail as expendable — the run is scientifically done once
`results/variant_calling/<caller>/<sample>/` holds the merged VCFs, so judge completion on that, not
on the driver exiting. Note `haplotypecaller_filter` (CNNScoreVariants) is **not** free to skip:
`consensus_from_results.sh` prefers the CNN-filtered HaplotypeCaller VCF, so skipping it changes what
feeds the consensus and would invalidate these numbers. Both points are recorded in `CLAUDE.md`.

There was **no "execution aborted"** in the log. The two `[ERROR] ... exit: 50001` entries you may
have seen were Spot preemptions of individual HaplotypeCaller shards — code 50001 is in the retry
class, they were retried, and they succeeded. 19 of 178 Batch jobs were preempted and recovered;
that is normal Spot behaviour and is why the run is cheap.

## 4. A bug this exercise caught in your pipeline

`consensus_from_results.sh` was **broken for relative invocation** — it did `cd "$(dirname "$0")"`
*before* sourcing `site.sh`, so calling it as `../consensus_from_results.sh` looked for `site.sh` in
the wrong directory and aborted. I introduced this yesterday when adding the portable-paths layer;
it is fixed and pushed, and I verified it now resolves correctly from absolute, relative and
symlinked invocations.

## 5. Cost and footprint (actual vs projected)

| | projected | actual |
|---|---|---|
| Compute | ~$6.40 | **$4.68** |
| Storage/transfer | ~$0.15 | ~$0.05 (public FASTQ read in place, never copied) |
| **Total** | ~$6.50 | **$4.73** — 32% of the $15 authorised |
| Wall clock | 4–8h | 11h 23m (QC tail, see §3) |
| Batch jobs | — | 178 (159 clean, 19 preempted+retried) |

The `work/` scratch directory has been deleted; `results/` (per-caller VCFs) is kept so the
consensus can be re-derived without re-calling.

## 6. What is now in the repo

`sarek-clinical/validation/` — tracked and pushed:

- `run_giab.sh` — end-to-end reproduction
- `benchmark_giab.sh` — RTG vcfeval per confidence tier, genome-wide and exome-restricted
- `health_check.sh` — progress/cost snapshot with a hard budget guard that cancels the run if spend
  crosses a ceiling
- `README.md` — how to obtain the truth set and re-run
- `RESULTS-HG002.md` — the measured numbers

README section 6 no longer says *"the pipeline must be formally validated"* — it now carries the
measured F1 and links to the results.

## 7. Honest limitations

- Valid **only inside GIAB high-confidence regions** — segmental duplications and many repeats are
  excluded, so real genome-wide sensitivity is lower.
- **One sample, one run** — a point estimate, not a reproducibility study across batches/operators.
- **SNVs and small indels only** — no structural variants, CNVs, repeat expansions or mosaicism.
- The exome row is WGS reads restricted to a capture BED, **not real capture data** — it does not
  reproduce the BGE coverage profile or duplicate structure.
- This validates **variant calling**, not the `candidate-filtering` interpretation step.
- Accreditation bodies normally expect replicates and a documented SOP; this is the technical
  evidence, not a completed clinical validation dossier.
