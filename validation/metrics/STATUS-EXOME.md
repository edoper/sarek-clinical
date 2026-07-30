# GIAB EXOME validation — live status

_updated 2026-07-28 07:10:20Z (checks run hourly)_

| | |
|---|---|
| Nextflow driver | not running |
| Tasks submitted | 66 |
| Batch jobs | 33 succeeded / 0 running / 6 failed |
| VMs up now | 0 |
| Spend so far | $0.78 of $15 authorised |
| Driver errors | 0 |
| Pipeline finished | yes |

```
cost Spot [#...................] $0.78 / $15 budget (5%) | 39 jobs, 0 run, 6 fail
```

Recent pipeline lines:
```
[2026-07-28T07:09:15Z] callers (presence-string order): deepvariant strelka freebayes haplotypecaller
[2026-07-28T07:09:16Z] isec sites: 118973
[2026-07-28T07:09:17Z] backbone (DeepVariant): 105851 variants
[2026-07-28T07:09:17Z] rescue via strelka: 2626 variants
[2026-07-28T07:09:18Z] rescue via haplotypecaller: 591 variants
[2026-07-28T07:09:18Z] FINAL consensus: /home/edo/sarek-clinical/validation/exome/consensus/HG002EX.consensus.vcf.gz  [109068 variants total]
[2026-07-28T07:09:18Z] consensus.sh done
ERROR: missing ./rtg-tools-3.12.1/rtg
```
