"""Build Supplementary_Tables.xlsx (Tables S1-S7) from the saved results.

    python scripts/make_supplement.py [results_dir] [metadata_dir] [out.xlsx]
"""
import json, sys
from pathlib import Path
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parents[1]
RES = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results"
META = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "metadata"
OUT = Path(sys.argv[3]) if len(sys.argv) > 3 else ROOT / "Supplementary_Tables.xlsx"

groups = pd.read_csv(META / "species_groups.csv", dtype=str, keep_default_na=False)
hosts = pd.read_csv(META / "host_mapping.csv", dtype=str, keep_default_na=False)
gname = {r.group_id: (f"{r.FINAL_COMMON} ({r.FINAL_SPECIES})" if r.FINAL_COMMON else r.FINAL_SPECIES) for r in groups.itertuples()}
host_of = {r.voucher: r.FINAL_HOST for r in hosts.itertuples()}


def iso(tip):  # 'PX634653_Kerman_Pb24' -> 'Pb24'
    return tip.split("_Kerman_")[-1]


def with_host(df, col):
    df = df.copy(); v = df[col].map(iso)
    df.insert(df.columns.get_loc(col) + 1, "isolate_host", v.map(host_of)); return df


wb = Workbook(); wb.remove(wb.active)
HEAD = PatternFill("solid", fgColor="D9E2F3")


def sheet(name, title, note, df, round_=4):
    ws = wb.create_sheet(name)
    ws["A1"] = title; ws["A1"].font = Font(bold=True, size=12)
    ws["A2"] = note; ws["A2"].alignment = Alignment(wrap_text=False)
    r0 = 4
    for j, c in enumerate(df.columns, 1):
        x = ws.cell(r0, j, str(c)); x.font = Font(bold=True); x.fill = HEAD
        x.alignment = Alignment(wrap_text=True, vertical="top")
    for i, row in enumerate(df.itertuples(index=False), r0 + 1):
        for j, v in enumerate(row, 1):
            if isinstance(v, float):
                v = None if pd.isna(v) else round(v, round_)
            elif v is not None and not isinstance(v, (int, str, bool)):
                v = str(v)
            ws.cell(i, j, v)
    for j, c in enumerate(df.columns, 1):
        w = max([len(str(c))] + [len(str(x)) for x in df.iloc[:, j - 1].head(200)])
        ws.column_dimensions[get_column_letter(j)].width = min(max(10, w + 2), 60)
    ws.freeze_panes = ws.cell(r0 + 1, 1)
    return ws


# S1 Sanger QC
q = pd.read_csv(RES / "tables" / "sanger_qc.csv")
q.insert(2, "host", q.voucher.map(host_of))
q = q.rename(columns={"seq_no": "sequencing_order_no", "voucher": "isolate", "raw_length": "raw_read_length_nt",
                      "trimmed_length": "trimmed_read_length_nt", "mean_phred": "mean_Phred", "q20_fraction": "fraction_Q20",
                      "q30_fraction": "fraction_Q30", "n_ambiguous_calls": "ambiguous_base_calls",
                      "submitted_length": "submitted_length_nt", "submitted_covered_by_read": "fraction_of_submitted_sequence_covered_by_read",
                      "identity_read_vs_submitted": "identity_read_vs_submitted", "gaps_read_vs_submitted": "gaps_read_vs_submitted",
                      "lowercase_edited_bases_in_raw": "lowercase_bases_in_raw_FASTA"})
sheet("S1_Sanger_QC", "Table S1. Quality of the forward Sanger reads and agreement with the submitted sequences",
      "Forward reads only. Coverage <0.95 or identity <0.97 (S7, S27, Pb64, Pb30) is flagged in the manuscript.", q)

# S2 survey statistics
prev = pd.read_csv(RES / "stats" / "prevalence_by_species.csv")
prev["host_group"] = prev.group.map(gname); prev = prev[["group", "host_group", "n", "positive", "prevalence", "ci_low", "ci_high", "low", "moderate", "high"]]
sex = pd.read_csv(RES / "stats" / "prevalence_by_sex.csv")
fr = pd.read_csv(RES / "stats" / "firth_logistic_regression.csv")
fr["term"] = fr.term.map(lambda t: gname.get(t, t) if t.startswith("SP") else t)
at = json.load(open(RES / "stats" / "association_tests.json"))
rows = []
for k, v in at.items():
    if not isinstance(v, dict):
        continue
    rows.append(dict(analysis=k, test=v.get("test", ""), statistic=v.get("chi2", v.get("odds_ratio")),
                     p_value=v.get("p", v.get("p_mc")), p_asymptotic=v.get("p_asymptotic"), cramers_v=v.get("cramers_v"),
                     min_expected_count=v.get("min_expected"), pct_cells_expected_lt5=v.get("pct_cells_expected_lt5")))
tests = pd.DataFrame(rows)
ws = sheet("S2_Survey_statistics", "Table S2. Prevalence, association tests and Firth logistic regression (136 birds)",
           "Monte Carlo exact tests: 100,000 permutations, margins fixed, seed 20261103. Wilson 95% CIs. Reference level for regression: Green-cheeked conure and females.", prev)
r = ws.max_row + 3
for title, df in (("Prevalence by sex", sex), ("Association tests", tests), ("Firth penalised logistic regression (positive ~ host group + sex)", fr)):
    ws.cell(r, 1, title).font = Font(bold=True); r += 1
    for j, c in enumerate(df.columns, 1):
        x = ws.cell(r, j, str(c)); x.font = Font(bold=True); x.fill = HEAD
    for row in df.itertuples(index=False):
        r += 1
        for j, v in enumerate(row, 1):
            ws.cell(r, j, None if (isinstance(v, float) and pd.isna(v)) else (round(v, 4) if isinstance(v, float) else v))
    r += 3

# S3 placement, identities, nearest references
pl = pd.read_csv(RES / "tables" / "placement_best.csv"); pl = with_host(pl, "query")
lw = pd.read_csv(RES / "tables" / "placement_lwr_by_clade.csv")
lw = lw.pivot_table(index="query", columns="clade", values="like_weight_ratio", aggfunc="sum").reset_index()
pl = pl.merge(lw, on="query", how="left", suffixes=("", "_LWR"))
cmp_ = pd.read_csv(RES / "tables" / "kerman_clade_assignment_comparison.csv").rename(columns={"tip": "query"})
pl = pl.merge(cmp_[["query", "amplicon_tree", "constrained_tree"]], on="query", how="left")
idc = pd.read_csv(RES / "tables" / "kerman_mean_identity_by_clade.csv").rename(columns={"Unnamed: 0": "query"})
nnb = pd.read_csv(RES / "tables" / "kerman_nearest_neighbours.csv"); nnb = with_host(nnb, "isolate")
nnb = nnb.rename(columns={"neighbour": "nearest_reference", "host": "reference_host", "country": "reference_country", "year": "reference_year", "clade": "reference_clade_Shah2023"})
ws = sheet("S3_Placement_identity", "Table S3. Phylogenetic placement, identity to sub-clades and nearest references of the Kerman sequences",
           "LWR = likelihood weight ratio (EPA-ng, 429-genome backbone). 'unassigned' = no uniformly labelled clade (see Methods).", pl)
for title, df in (("Mean % identity of each Kerman sequence to each sub-clade", idc), ("Ten nearest references per Kerman sequence (uncorrected p-distance)", nnb)):
    r = ws.max_row + 3; ws.cell(r, 1, title).font = Font(bold=True); r += 1
    for j, c in enumerate(df.columns, 1):
        x = ws.cell(r, j, str(c)); x.font = Font(bold=True); x.fill = HEAD
    for row in df.itertuples(index=False):
        r += 1
        for j, v in enumerate(row, 1):
            ws.cell(r, j, None if (isinstance(v, float) and pd.isna(v)) else (round(v, 4) if isinstance(v, float) else v))

# S4 diversity and trait association
dv = pd.read_csv(RES / "tables" / "diversity_haplotypes.csv"); ta = pd.read_csv(RES / "tables" / "trait_association.csv")
ws = sheet("S4_Diversity_association", "Table S4. Genetic diversity and phylogeny-trait association", "1,000 random permutations of tip states; the smallest attainable p is 0.001.", dv)
r = ws.max_row + 3; ws.cell(r, 1, "Phylogeny-trait association (BaTS-like)").font = Font(bold=True); r += 1
for j, c in enumerate(ta.columns, 1):
    x = ws.cell(r, j, str(c)); x.font = Font(bold=True); x.fill = HEAD
for row in ta.itertuples(index=False):
    r += 1
    for j, v in enumerate(row, 1):
        ws.cell(r, j, round(v, 4) if isinstance(v, float) else v)

# S5 selection and motifs
sel = pd.read_csv(RES / "tables" / "selection_sites.csv")
fel = sel[sel.method == "FEL"][["codon", "alpha", "beta", "LRT", "p-value"]].rename(columns={"alpha": "FEL_alpha(dS)", "beta": "FEL_beta(dN)", "LRT": "FEL_LRT", "p-value": "FEL_p"})
mem = sel[sel.method == "MEME"][["codon", "&alpha;", "&beta;<sup>+</sup>", "p<sup>+</sup>", "LRT", "p-value", "# branches under selection"]]
mem.columns = ["codon", "MEME_alpha", "MEME_beta_plus", "MEME_p_plus", "MEME_LRT", "MEME_p", "MEME_branches_under_selection"]
site = fel.merge(mem, on="codon", how="left")
site["selection_FEL_p<=0.1"] = ((site.FEL_p <= 0.1) & (site["FEL_beta(dN)"] < site["FEL_alpha(dS)"])).map({True: "purifying", False: ""})
mot = pd.read_csv(RES / "tables" / "rep_motifs_kerman.csv"); mot.insert(1, "host", mot.sequence.map(iso).map(host_of))
ws = sheet("S5_Selection_motifs", "Table S5. Site-wise selection (HyPhy FEL and MEME) and Rep motifs", "217 codons of the Rep fragment. Significance threshold p <= 0.1. SLAC global dN/dS = 0.14. GARD: no breakpoint.", site)
r = ws.max_row + 3; ws.cell(r, 1, "Conserved Rep motifs (@ = position in the 217-aa fragment)").font = Font(bold=True); r += 1
for j, c in enumerate(mot.columns, 1):
    x = ws.cell(r, j, str(c)); x.font = Font(bold=True); x.fill = HEAD
for row in mot.itertuples(index=False):
    r += 1
    for j, v in enumerate(row, 1): ws.cell(r, j, v)

# S6 isolates and host names
h = hosts[["seq_no", "voucher", "accession", "FINAL_HOST", "host_genbank_sheet", "host_thesis", "s_code_match", "evidence"]].copy()
h.columns = ["sequencing_order_no", "isolate", "GenBank_accession", "host_final", "host_in_GenBank_sheet", "host_in_thesis", "match_to_survey_code", "evidence_from_thesis"]
sheet("S6_Isolates", "Table S6. Sequenced isolates and host species", "Final host names supplied by the authors. The Pb vouchers could not be matched individually to the 136 surveyed birds.", h)

# S7 references
ref = pd.read_csv(RES / "tables" / "reference_metadata_merged.csv")
tc = pd.read_csv(RES / "tables" / "tip_clades_backbone.csv").rename(columns={"tip": "accession", "clade": "backbone_clade"})[["accession", "backbone_clade"]]
ref = ref.merge(tc, on="accession", how="left")
cols = [c for c in ["accession", "length", "complete_genome", "host", "country_simple", "year", "collection_date", "isolate", "isolation_source", "shah2023_clade", "backbone_clade", "definition"] if c in ref.columns]
sheet("S7_References", "Table S7. The 589 reference sequences", "Metadata from GenBank (NCBI Entrez). Sub-clade labels transcribed from Shah et al. (2023) and propagated as described in Methods.", ref[cols])

# README
ws = wb.create_sheet("README", 0)
for i, t in enumerate(["Supplementary Tables for: Subclinical beak and feather disease virus infection in pet psittacines in Kerman, Iran",
                       "S1 Sanger read quality | S2 survey statistics | S3 placement and identity | S4 diversity and association | S5 selection and motifs | S6 isolates | S7 references",
                       "All tables are generated by scripts/make_supplement.py from the files in results/ (code repository included with the submission)."], 1):
    ws.cell(i, 1, t)
ws["A1"].font = Font(bold=True, size=12); ws.column_dimensions["A"].width = 140
wb.save(OUT); print("wrote", OUT)
