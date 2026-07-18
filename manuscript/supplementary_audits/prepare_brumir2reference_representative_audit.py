#!/usr/bin/env python3

Purpose
-------
Generate audit tables summarizing BrumiR2Reference annotation
status for every CD-HIT representative sequence.

This script was used during manuscript preparation to inspect
core clusters, representative annotations, passfilter status,
and DESeq2 statistics. It is not part of the nf-Sarcopipe
workflow.

from __future__ import annotations
import argparse, csv, math, re
from collections import defaultdict
from pathlib import Path

CLSTR_RE = re.compile(r">(.+?)\.\.\.")
FIELDS = ["chr","start","stop","MFE","B","E","H","I","K","M","S","X","SEGMENTS","Precursor_Seq"]

def norm(s):
    return re.sub(r"\s+", "", s).upper().replace("T", "U")

def fasta(path):
    out, name, buf = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None: out[name] = norm("".join(buf))
                name, buf = line[1:].split()[0], []
            else: buf.append(line)
        if name is not None: out[name] = norm("".join(buf))
    return out

def clstr(path):
    members, reps, current = defaultdict(list), {}, None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith(">Cluster"):
                current = f"cluster_{line.split()[-1]}"; continue
            if current is None: continue
            m = CLSTR_RE.search(line)
            if not m: continue
            member = m.group(1).strip(); members[current].append(member)
            if line.endswith("*"): reps[current] = member
    for cid, vals in members.items():
        reps.setdefault(cid, vals[0])
    return dict(members), reps

def core_sets(path):
    out = defaultdict(set)
    with open(path, encoding="utf-8", errors="replace") as fh:
        next(fh, None)
        for raw in fh:
            p = raw.rstrip("\n").split("\t")
            if len(p) >= 2: out[p[1].strip()].add(p[0].strip())
    return dict(out)

def group(labels):
    a, s = "athlete_core" in labels, "sedentary_core" in labels
    return "shared" if a and s else "athlete" if a else "sedentary" if s else "not_core"

def read_b2r(path, status):
    out = defaultdict(list)
    if not path or not Path(path).exists(): return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        idx = {c.strip(): i for i,c in enumerate(header)}
        needed = {"miRNA", *FIELDS}
        miss = needed - set(idx)
        if miss: raise SystemExit(f"Missing columns in {path}: {sorted(miss)}")
        for raw in fh:
            if not raw.strip() or raw.startswith("#"): continue
            p = raw.rstrip("\n").split("\t")
            if len(p) <= max(idx.values()): continue
            member = p[idx["miRNA"]].strip()
            row = {"status":status, "brumir2ref_miRNA":member}
            for f in FIELDS:
                value = p[idx[f]].strip()
                row["precursor_seq" if f=="Precursor_Seq" else f] = norm(value) if f=="Precursor_Seq" else value
            out[member].append(row)
    return out

def deseq(path):
    if not path: return {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return {r["feature"]:r for r in csv.DictReader(fh)}

def num(x):
    try:
        y=float(x); return y if math.isfinite(y) else None
    except: return None

def best(pass_hits, nonpass_hits):
    hits = pass_hits if pass_hits else nonpass_hits
    if not hits: return None
    return min(hits, key=lambda r:(num(r["MFE"]) is None, num(r["MFE"]) or 0.0))

def write(path, rows, cols):
    with open(path,"w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--clustered_fasta",required=True)
    ap.add_argument("--clstr",required=True)
    ap.add_argument("--core_sets_tsv",required=True)
    ap.add_argument("--passfilter_tsv",required=True)
    ap.add_argument("--nonpassfilter_tsv")
    ap.add_argument("--brumir_deseq2")
    ap.add_argument("--out_representatives",default="brumir2ref_representative_report.tsv")
    ap.add_argument("--out_all_hsps",default="brumir2ref_representative_all_hsps.tsv")
    ap.add_argument("--out_summary",default="brumir2ref_representative_summary.tsv")
    a=ap.parse_args()

    seqs=fasta(a.clustered_fasta); members,reps=clstr(a.clstr); cores=core_sets(a.core_sets_tsv)
    pf=read_b2r(a.passfilter_tsv,"passfilter"); npf=read_b2r(a.nonpassfilter_tsv,"nonpassfilter")
    de=deseq(a.brumir_deseq2)
    rows=[]; allrows=[]; summ=defaultdict(int); gs=defaultdict(lambda:defaultdict(int))
    cols=["cluster","group","core_sets","is_core","representative_member","representative_sequence","representative_length","cluster_size","status","n_passfilter_hsps","n_nonpassfilter_hsps","brumir2ref_miRNA","chr","start","stop","MFE","B","E","H","I","K","M","S","X","SEGMENTS","precursor_len","representative_in_precursor","mature_start","mature_end","precursor_seq","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","is_DE_padj_lt_0.05"]

    for cid in sorted(reps, key=lambda x:int(x.split("_")[1])):
        rep=reps[cid]; seq=seqs.get(rep,""); labels=cores.get(cid,set()); grp=group(labels); iscore=grp!="not_core"
        ph=list(pf.get(rep,[])); nh=list(npf.get(rep,[])); hit=best(ph,nh)
        status="passfilter" if ph else "nonpassfilter" if nh else "unmapped"
        d=de.get(cid,{}); padj=num(d.get("padj","")); isde=padj is not None and padj<0.05
        summ["clusters_total"]+=1; summ[f"status_{status}"]+=1
        if iscore:
            summ["core_clusters_total"]+=1; summ[f"core_status_{status}"]+=1
            gs[grp]["core_total"]+=1; gs[grp][status]+=1
        base={"cluster":cid,"group":grp,"core_sets":",".join(sorted(labels)) if labels else "NA","is_core":"1" if iscore else "0","representative_member":rep,"representative_sequence":seq or "NA","representative_length":str(len(seq)) if seq else "NA","cluster_size":str(len(members.get(cid,[]))),"status":status,"n_passfilter_hsps":str(len(ph)),"n_nonpassfilter_hsps":str(len(nh)),"baseMean":d.get("baseMean","NA"),"log2FoldChange":d.get("log2FoldChange","NA"),"lfcSE":d.get("lfcSE","NA"),"stat":d.get("stat","NA"),"pvalue":d.get("pvalue","NA"),"padj":d.get("padj","NA"),"is_DE_padj_lt_0.05":"1" if isde else "0"}
        def merged(h):
            r=dict(base)
            if h is None:
                r.update({k:"NA" for k in ["brumir2ref_miRNA","chr","start","stop","MFE","B","E","H","I","K","M","S","X","SEGMENTS","precursor_len","mature_start","mature_end","precursor_seq"]}); r["representative_in_precursor"]="0"
            else:
                pre=h["precursor_seq"]; pos=pre.find(seq) if seq else -1
                r.update(h); r["precursor_len"]=str(len(pre)); r["representative_in_precursor"]="1" if pos>=0 else "0"; r["mature_start"]=str(pos+1) if pos>=0 else "NA"; r["mature_end"]=str(pos+len(seq)) if pos>=0 else "NA"
            return r
        rows.append(merged(hit))
        for h in ph+nh: allrows.append(merged(h))

    write(a.out_representatives,rows,cols); write(a.out_all_hsps,allrows,cols)
    with open(a.out_summary,"w",encoding="utf-8") as fh:
        fh.write("scope\tgroup\tmetric\tcount\n")
        for k in sorted(summ): fh.write(f"global\tall\t{k}\t{summ[k]}\n")
        for grp in ("athlete","sedentary","shared"):
            for k in ("core_total","passfilter","nonpassfilter","unmapped"):
                fh.write(f"core_group\t{grp}\t{k}\t{gs[grp][k]}\n")
    print(f"[OK] {a.out_representatives}")
    print(f"[OK] {a.out_all_hsps}")
    print(f"[OK] {a.out_summary}")

if __name__=="__main__": main()
