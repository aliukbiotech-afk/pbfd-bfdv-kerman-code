from Bio import SeqIO, Align
kerman = P.read_fasta(DATA / "kerman_10_unaligned.fasta")
hosts = pd.read_csv(META / "host_mapping.csv", dtype=str, keep_default_na=False)
al = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1)
rows = []
for ab in sorted((DATA / "ab1").glob("*.ab1")):
    rec = SeqIO.read(ab, "abi"); q = np.array(rec.letter_annotations["phred_quality"])
    tr = SeqIO.read(ab, "abi-trim")
    stem = ab.stem; no = int(stem.split("-")[0])
    h = hosts[hosts.seq_no == str(no)].iloc[0]
    sub = next(v for k, v in kerman.items() if k.startswith(h.accession))
    best = max((al.align(str(x), sub)[0] for x in (tr.seq, tr.seq.reverse_complement())), key=lambda a: a.score)
    cnt = best.counts()
    cov = sum(b1 - b0 for b0, b1 in best.aligned[1]) / len(sub)
    rows.append(dict(seq_no=no, voucher=h.voucher, accession=h.accession, raw_length=len(rec), trimmed_length=len(tr),
                     mean_phred=q.mean(), q20_fraction=(q >= 20).mean(), q30_fraction=(q >= 30).mean(),
                     n_ambiguous_calls=sum(c not in "ACGT" for c in str(rec.seq)),
                     submitted_length=len(sub), submitted_covered_by_read=cov,
                     identity_read_vs_submitted=cnt.identities / max(1, cnt.identities + cnt.mismatches),
                     gaps_read_vs_submitted=cnt.gaps))
qc = pd.DataFrame(rows).sort_values("seq_no"); qc.to_csv(RES / "tables" / "sanger_qc.csv", index=False)
raw = P.read_fasta(DATA / "kerman_10_raw_basecalls.fasta")
qc["lowercase_edited_bases_in_raw"] = [sum(c.islower() for c in next((v for k, v in raw.items() if k.startswith(f"{n}-")), "")) for n in qc.seq_no]
qc.to_csv(RES / "tables" / "sanger_qc.csv", index=False); qc.round(3)