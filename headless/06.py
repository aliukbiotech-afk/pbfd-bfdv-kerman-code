thesis_aln = P.read_fasta(DATA / "thesis_alignment_598.fasta")
REFS = [k for k in thesis_aln if "_Kerman_" not in k]
print(len(REFS), "reference accessions used in the thesis tree")
gbmd = P.entrez_fetch(REFS, EMAIL, RES / "tables")
old = pd.read_csv(META / "reference_metadata.csv")
mdm = gbmd.merge(old[["accession", "isolation_source"]], on="accession", how="left")
shah = pd.read_csv(META / "shah2023_clades.csv")
mdm = mdm.merge(shah[["accession", "clade"]].rename(columns={"clade": "shah2023_clade"}), on="accession", how="left")
mdm.to_csv(RES / "tables" / "reference_metadata_merged.csv", index=False)
print("records:", len(mdm), "| complete genomes:", mdm.complete_genome.sum(), "| with collection year:", mdm.year.notna().sum(),
      "| with Shah-2023 clade:", mdm.shah2023_clade.notna().sum())
print(mdm.country_simple.value_counts().head(15))