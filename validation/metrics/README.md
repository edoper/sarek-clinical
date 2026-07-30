# Raw validation metrics — HG002 / GIAB v4.2.1, GRCh38

The vcfeval output that `../RESULTS-HG002.md` and `../RESULTS-HG002-EXOME.md` summarise.
Committed so every published accuracy figure has a checkable source. HG002 is a public
reference sample, so nothing here is patient data.

## Which file answers which question

| File | Run | Use it for |
|---|---|---|
| `metrics-genomewide.txt` | WGS | The WGS arm, all four consensus tiers |
| `metrics-panel.txt` | WGS, panel-restricted | Panel behaviour of the WGS arm |
| `metrics-exome.txt` | **Real exome** | **The exome and BGE clinical arms** |
| `metrics-all-tiers-WGS-PROXY.txt` | WGS proxy | Superseded — see below |
| `vcfeval-*.txt` | both | Per-tier vcfeval summaries, unaggregated |
| `HG002.consensus.log` | WGS | Caller inputs and rescue counts for the consensus |
| `REPORT.md` | both | Run report: cost, what was run, what it means |
| `STATUS-EXOME.md` | exome | Live-status snapshot from the end of the exome run |

## The one trap in here

`metrics-all-tiers-WGS-PROXY.txt` reports the exome as indel F1 **0.9922**. That is the WGS
run restricted to a capture BED, not the exome path, and it overstates indel performance by
roughly 0.07 F1. The real exome run gives **0.9198** in the same panel genes
(`metrics-exome.txt`, `fairpanel-union`). Commit `70ce62b` retracted the proxy method; the
file is kept for provenance, with a warning header, not for quoting.

The `union-exome` rows in `metrics-genomewide.txt` come from that same proxy — they are rows
of the genomewide table, so read them as WGS-restricted, never as exome performance.

For the clinical arms, quote `metrics-exome.txt`.
