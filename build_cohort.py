#!/usr/bin/env python3
# Build families.tsv + a CRAM staging list from a Terra sample-table export.
# Role: proband = no suffix; <id>M = madre(-M); <id>P = padre(-F).
import csv, re, sys
TSV = sys.argv[1]
DST = "gs://intergenica-sarek-clinical/bge-wes/crams-cohort"
def classify(cid):
    m = re.match(r"^(.+\d)M$", cid)
    if m: return m.group(1), "M"
    m = re.match(r"^(.+\d)P$", cid)
    if m: return m.group(1), "F"
    return cid, "P"
rows = list(csv.DictReader(open(TSV), delimiter="\t"))
fam_tsv = open("families.tsv", "w")
fam_tsv.write("# auto-built from %s — family\\trole\\tsex\\tcram_uri (staged)\n" % TSV.split("/")[-1])
copy_list = open("/tmp/cram_srcs.txt", "w")
n = 0
for r in rows:
    cid = (r.get("collaborator_sample_id") or "").strip()
    cram = (r.get("genome_cram_path") or "").strip()
    crai = (r.get("genome_crai_path") or "").strip()
    sex  = (r.get("predicted_sex") or "NA").strip()
    if not cid or not cram.startswith("gs://"):
        sys.stderr.write("SKIP (no cram): %s\n" % cid); continue
    fam, role = classify(cid)
    staged = "%s/%s.cram" % (DST, cid)        # basename preserved on copy
    fam_tsv.write("%s\t%s\t%s\t%s\n" % (fam, role, sex, staged))
    copy_list.write(cram + "\n"); copy_list.write(crai + "\n")
    n += 1
fam_tsv.close(); copy_list.close()
print("wrote families.tsv (%d samples) and /tmp/cram_srcs.txt (%d files to copy)" % (n, n*2))
