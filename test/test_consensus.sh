#!/usr/bin/env bash
#
# test_consensus.sh — regression test for the two load-bearing behaviours of the
# consensus stage. Synthetic data only (no patient data, no cloud, no network);
# runs in ~30s on a laptop. Requires bcftools, bgzip, tabix, samtools, python3.
#
#   1. CONSENSUS RULE — DeepVariant backbone + rescue of variants >=2 OTHER callers
#      agreed on, genotype borrowed from Strelka2 else HaplotypeCaller, single-caller
#      non-DeepVariant sites dropped. Guards the invariants documented in consensus.sh
#      (notably `--rm-dup exact`, which a well-meaning switch to `all` would break by
#      collapsing just-split multiallelics).
#
#   2. CRASH SAFETY — consensus.sh must never leave a truncated VCF at its FINAL output
#      path, and the resume guard in consensus_from_results.sh must reject anything that
#      is not whole+indexed. A partial file at the final path is the dangerous failure:
#      a resume reads it as "done" and ships an incomplete variant list to VEP and
#      candidate-filtering, where nothing downstream would reveal the missing variants.
#
# Usage:  ./test/test_consensus.sh          (or: CONSENSUS_SH=/path/to/consensus.sh ...)
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONSENSUS_SH="${CONSENSUS_SH:-$REPO/consensus.sh}"
[[ -x "$CONSENSUS_SH" ]] || { echo "ERROR: consensus.sh not executable: $CONSENSUS_SH"; exit 1; }
for t in bcftools bgzip tabix samtools python3; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH (source env.sh / conda activate)"; exit 1; }
done

TD="$(mktemp -d)"; trap 'rm -rf "$TD"' EXIT
cd "$TD"
fails=0
ok()   { echo "  PASS  $*"; }
bad()  { echo "  FAIL  $*"; fails=$((fails+1)); }

# The completeness test must stay identical to complete_vcf() in consensus_from_results.sh.
BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
complete_vcf() {
    [[ -s "$1" && -s "$1.tbi" ]] || return 1
    [[ "$(tail -c 28 "$1" | od -An -tx1 | tr -d ' \n')" == "$BGZF_EOF" ]]
}

# ───────────────────────── test 1: consensus rule ─────────────────────────
echo "== test 1: consensus rule (backbone + rescue + drops) =="
{ echo ">chr1"
  python3 -c "
import random
random.seed(7)
s=''.join(random.choice('ACGT') for _ in range(1000))
print('\n'.join(s[i:i+60] for i in range(0,1000,60)))"; } > ref.fa
samtools faidx ref.fa
b() { samtools faidx ref.fa "chr1:$1-$1" | tail -1; }                              # REF base
a() { case "$(b "$1")" in A) echo T;; T) echo A;; C) echo G;; G) echo C;; esac; }   # ALT != REF
hdr() {
  cat <<EOF
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FILTER=<ID=PASS,Description="passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
$1#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
EOF
}
VAF='##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Alt allele fraction">
'
row() { printf 'chr1\t%s\t.\t%s\t%s\t50\tPASS\t.\t%s\t%s\n' "$1" "$(b "$1")" "$(a "$1")" "$2" "$3"; }

# 100 = all four | 200 = DeepVariant only | 300 = STR+HC+FB (rescue via strelka)
# 400 = HC+FB (rescue via haplotypecaller) | 500 = FB only, 600 = STR only  -> both DROPPED
{ hdr "$VAF"; row 100 "GT:GQ:DP:AD:VAF" "0/1:99:30:15,15:0.5"
              row 200 "GT:GQ:DP:AD:VAF" "1/1:80:25:0,25:1.0"; } | bgzip -c > dv.vcf.gz
{ hdr "";     row 100 "GT:GQ:DP:AD" "0/1:70:28:14,14"
              row 300 "GT:GQ:DP:AD" "0/1:60:20:10,10"
              row 600 "GT:GQ:DP:AD" "0/1:55:18:9,9"; }  | bgzip -c > str.vcf.gz
{ hdr "";     row 100 "GT:GQ:DP:AD" "0/1:65:27:13,14"
              row 300 "GT:GQ:DP:AD" "0/1:58:19:9,10"
              row 400 "GT:GQ:DP:AD" "0/1:52:16:8,8"; }  | bgzip -c > hc.vcf.gz
{ hdr "";     row 100 "GT:GQ:DP:AD" "0/1:40:26:13,13"
              row 300 "GT:GQ:DP:AD" "0/1:38:18:9,9"
              row 400 "GT:GQ:DP:AD" "0/1:35:15:7,8"
              row 500 "GT:GQ:DP:AD" "0/1:30:12:6,6"; }  | bgzip -c > fb.vcf.gz
for f in dv str hc fb; do tabix -f -p vcf "$f.vcf.gz"; done

mkdir -p out
"$CONSENSUS_SH" -r ref.fa -d dv.vcf.gz -o out/S1 \
  -c strelka=str.vcf.gz -c freebayes=fb.vcf.gz -c haplotypecaller=hc.vcf.gz >/dev/null 2>&1
bcftools query -f '%POS\t%INFO/NCALLERS\t%INFO/CONF\t%INFO/GT_SOURCE\t[%GT]\n' out/S1.consensus.vcf.gz > got.txt
cat > want.txt <<'EOF'
100	4	HIGH	deepvariant	0/1
200	1	LOW	deepvariant	1/1
300	3	HIGH	strelka	0/1
400	2	MEDIUM	haplotypecaller	0/1
EOF
if diff -u want.txt got.txt >/dev/null; then ok "backbone+rescue correct; single-caller sites 500/600 dropped"
else bad "consensus content mismatch:"; diff -u want.txt got.txt || true; fi
[[ -s out/S1.consensus.vcf.gz.tbi ]]      && ok "output indexed"          || bad "no .tbi written"
[[ ! -e out/S1.consensus.vcf.gz.partial ]] && ok "no .partial left behind" || bad ".partial leaked on success"

# ───────────────────── test 2: completeness guard ─────────────────────
echo "== test 2: resume completeness guard =="
complete_vcf out/S1.consensus.vcf.gz && ok "whole+indexed VCF accepted" || bad "complete VCF rejected"
cp out/S1.consensus.vcf.gz T.vcf.gz; cp out/S1.consensus.vcf.gz.tbi T.vcf.gz.tbi
truncate -s -40 T.vcf.gz
complete_vcf T.vcf.gz && bad "truncated VCF ACCEPTED (would ship partial data)" \
                      || ok "truncated VCF rejected (bare [[ -s ]] would have accepted it)"
cp out/S1.consensus.vcf.gz U.vcf.gz     # data present, never indexed
complete_vcf U.vcf.gz && bad "unindexed VCF accepted" || ok "unindexed VCF rejected"

# ─────────────────── test 3: no partial at the final path ───────────────────
# Needs enough variants that the final concat|sort write window is observable.
# Kill the whole process group the instant the first output file appears — that is
# mid-write by construction, so this is deterministic rather than a timing race.
echo "== test 3: crash mid-write leaves the final path clean =="
python3 - <<'PY'
import random
random.seed(11); N=600_000
s=''.join(random.choice('ACGT') for _ in range(N))
with open('bigref.fa','w') as f:
    f.write('>chr1\n')
    for i in range(0,N,60): f.write(s[i:i+60]+'\n')
flip={'A':'T','T':'A','C':'G','G':'C'}
hdr=('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=%d>\n##FILTER=<ID=PASS,Description="p">\n'
     '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
     '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
     '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
     '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')%N
pos=list(range(100,N-100,4))
for name,frac in (('bigdv',1.0),('bigstr',0.9),('bighc',0.85),('bigfb',0.8)):
    with open(name+'.vcf','w') as f:
        f.write(hdr)
        for p in pos:
            if frac>=1.0 or random.random()<frac:
                f.write(f"chr1\t{p}\t.\t{s[p-1]}\t{flip[s[p-1]]}\t50\tPASS\t.\tGT:GQ:DP:AD\t0/1:60:20:10,10\n")
PY
samtools faidx bigref.fa
for f in bigdv bigstr bighc bigfb; do bgzip -f "$f.vcf" && tabix -f -p vcf "$f.vcf.gz"; done

mkdir -p out3
setsid "$CONSENSUS_SH" -r bigref.fa -d bigdv.vcf.gz -o out3/S1 \
  -c strelka=bigstr.vcf.gz -c freebayes=bigfb.vcf.gz -c haplotypecaller=bighc.vcf.gz >/dev/null 2>&1 &
p=$!; i=0
while (( i < 4000 )); do
    compgen -G "out3/S1.consensus.vcf.gz*" >/dev/null && { kill -9 -- -"$p" 2>/dev/null || true; break; }
    sleep 0.005; i=$((i+1))
done
wait "$p" 2>/dev/null || true
sleep 0.3
if   compgen -G "out3/S1.consensus.vcf.gz*" >/dev/null && [[ ! -e out3/S1.consensus.vcf.gz ]]; then
     ok "killed mid-write: only .partial present, final path clean"
elif [[ -e out3/S1.consensus.vcf.gz ]] && complete_vcf out3/S1.consensus.vcf.gz; then
     echo "  SKIP  run completed before the kill landed — inconclusive (rerun on a slower box)"
elif [[ -e out3/S1.consensus.vcf.gz ]]; then
     bad "TRUNCATED VCF AT FINAL PATH ($(stat -c%s out3/S1.consensus.vcf.gz) bytes) — a resume would ship it"
else bad "no output produced at all — test did not exercise the write window"; fi

echo
if (( fails == 0 )); then echo "ALL TESTS PASSED"; else echo "$fails TEST(S) FAILED"; exit 1; fi
