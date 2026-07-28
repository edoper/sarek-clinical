# Efficiency evaluation — measured from the GIAB WGS run

Source: the Nextflow execution trace of the 2026-07-27 HG002 validation
(`results/pipeline_info/execution_trace_*.txt`), 178 tasks, one 30× WGS genome.
Everything below is measured, not modelled.

## Headline

| | |
|---|---|
| Wall clock | **11.4 h** |
| Total task realtime | 51.7 h (summed across parallel tasks) |
| Total CPU-hours | **281** |
| Average parallelism achieved | **4.5×** (concurrency cap was 40) |
| Cost | **$4.68** Spot |
| Preemption waste | 3.93 h realtime, **7.6%** — 19 of 178 attempts, all recovered |

**Cost efficiency is excellent. Wall-clock efficiency is poor.** Those are different problems with
different fixes, and conflating them is the main trap here.

## Where the MONEY goes (CPU-hours)

| stage | CPU-h | share | cores used |
|---|---|---|---|
| BWA-MEM alignment | 109.5 | **39.0%** | 22.3 |
| DeepVariant | 94.9 | **33.8%** | 11.2 |
| HaplotypeCaller | 43.7 | 15.5% | 2.0 |
| Strelka2 | 11.9 | 4.2% | 5.4 |
| CNNScoreVariants | 7.5 | 2.7% | 1.9 |
| FreeBayes | 5.2 | 1.9% | 0.9 |
| MarkDuplicates | 3.1 | 1.1% | 1.3 |

Alignment plus DeepVariant is **73% of the bill**. The single biggest cost lever remains the one
README §6 already identifies: **start from an aligned CRAM and 39% of the cost disappears.**

## Where the TIME goes (task realtime) — a different ranking

| stage | realtime h | share of realtime | cores used |
|---|---|---|---|
| **HaplotypeCaller** | **21.53** | **41.6%** | **2.0** |
| DeepVariant | 8.46 | 16.3% | 11.2 |
| FreeBayes | 6.11 | 11.8% | 0.9 |
| BWA-MEM | 4.92 | 9.5% | 22.3 |
| **CNNScoreVariants** | **3.88** | **7.5%** | 1.9 |

**The two stages that dominate wall clock are the two that barely use their VMs.** HaplotypeCaller
occupies 41.6% of all task-time on 2 cores; FreeBayes runs effectively single-threaded (0.9 cores).
Compare BWA at 22.3 cores and DeepVariant at 11.2 — those finish fast because they actually
parallelise.

## The three concrete inefficiencies

**1. CNNScoreVariants is a 4-hour serial tail.** One task, 3 h 55 m, 2 cores — **34% of the entire
run's wall clock for 2.7% of its cost**, executing after every other stage had finished. This is
what makes a completed run look hung. It runs because `haplotypecaller_filter` is absent from the
usual `--skip_tools` list.

*It is not free to remove:* `consensus_from_results.sh` prefers `*.haplotypecaller.filtered.vcf.gz`,
so skipping it feeds raw HaplotypeCaller calls into the consensus. The validated F1 was measured
**with** it. Removing it is a legitimate ~4 h saving that requires re-validation.

**2. HaplotypeCaller is under-parallelised.** 23 interval shards at 2 cores each. DeepVariant
covers the same genome in 8.5 h of realtime using 11 cores per shard. More shards, or more threads
per shard, would cut the longest pole — and HaplotypeCaller's role in the consensus is only as a
rescue genotype donor, which the validation showed contributes little.

**3. Parallelism sat at 4.5× against a cap of 40.** One sample cannot fill the cluster: sarek shards
calling ~26 ways, so that is the ceiling for a single genome. **This pipeline is far more efficient
per sample on cohorts than on singletons** — the fixed serial tail is amortised across samples.
Validating or reprocessing one genome at a time is the expensive way to use it.

## Unclaimed win: the reference is pre-staged but only one config uses it

`gcb-bge-wes.config` points its `fasta` / `fasta_fai` / `dict` / `dbsnp` at
`$SAREK_BUCKET/refs/GATK.GRCh38/` — same-region, instant. **`gcb.config` does not**, so both
from-FASTQ arms (WGS and the EPIGEN exome) fall back to igenomes and pull the reference over the
network from `s3://ngi-igenomes` on every run.

Verified during the exome validation: the bucket **already contains everything needed**, including
the BWA index the from-FASTQ arms require —

```
refs/GATK.GRCh38/Homo_sapiens_assembly38.fasta          3.03 GB
refs/GATK.GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz 1.45 GB
refs/GATK.GRCh38/BWAIndex/Homo_sapiens_assembly38.fasta.64.{amb,ann,bwt,pac,sa}
```

…yet the exome run's log shows it staging `Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta`
and `Sequence/BWAIndex/` from S3 anyway. The cost is wall clock, not dollars — but it is paid again
on **every Spot preemption retry** of an alignment task, which is exactly when a run already looks
stalled.

**Suggested fix** (one config change, needs one verification first):

```groovy
// in gcb.config params { }
fasta     = "${gcpBucket}/refs/GATK.GRCh38/Homo_sapiens_assembly38.fasta"
fasta_fai = "${gcpBucket}/refs/GATK.GRCh38/Homo_sapiens_assembly38.fasta.fai"
dict      = "${gcpBucket}/refs/GATK.GRCh38/Homo_sapiens_assembly38.dict"
bwa       = "${gcpBucket}/refs/GATK.GRCh38/BWAIndex/"
```

**Verify before adopting:** confirm the bucket FASTA is byte-identical to the igenomes one
(`gcloud storage ls -l` both and compare size/checksum). If it differs, alignment differs, and the
GIAB numbers — which were measured with the igenomes reference — would no longer describe the
pipeline. This was deliberately NOT applied automatically for that reason.

## What NOT to change

- **FreeBayes** looks wasteful (11.8% of realtime, 0.9 cores) and is only 1.9% of cost. But it feeds
  the ≥2-caller concordance logic. Given the validation found the rescue arm neutral in-panel, its
  removal is defensible — but it is a science change, not an efficiency fix, and needs re-validation.
- **Spot.** 7.6% rework for roughly 5× cheaper VMs is a clearly good trade. All 19 preemptions
  recovered automatically.

## Efficiency of the validation exercise itself

| | projected | actual |
|---|---|---|
| Cost | ~$6.50 | **$4.68** (31% of the $15 authorised) |
| Wall clock | 4–8 h | 11.4 h |

Two decisions that paid off: reading the public HG002 FASTQ **in place** from `gs://brain-genomics-public/`
after verifying it is `allUsers`-readable (saved a 48.7 GB copy, its storage and its egress), and
extracting the DeepVariant-alone comparison from the `CALLERS` tag in the existing consensus VCF
rather than re-running anything (that comparison cost **$0** and overturned a published conclusion).

The estimate missed because it modelled the analytical work and ignored the serial tail — the
lesson recorded in `CLAUDE.md`: judge completion by whether
`results/variant_calling/<caller>/<sample>/` holds the merged VCFs, not by the driver exiting.
