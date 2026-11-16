"""Phase 1: build clean, analysis-ready inputs from the original thesis files.

Originals under ../PBFD are only read, never modified.
Run from the project root:  python analysis/scripts/prepare_inputs.py
"""
import csv, re, json, collections
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "PBFD"
OUT = ROOT / "analysis"
DATA, META = OUT / "data", OUT / "metadata"
DATA.mkdir(parents=True, exist_ok=True); META.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- Kerman isolates
# Sequencing order number (from .ab1 file prefixes) -> voucher -> accession
ISOLATES = [
    # seq_no, voucher, accession, ab1 stem
    (19, "Pb24", "PX634653", "19-Pb24  PBFD F  F"),
    (20, "Pb25", "PX634654", "20-Pb25  PBFD F  F"),
    (21, "Pb30", "PX634655", "21-p30  PBFD F  F"),
    (22, "S7",   "PX634656", "22-s7  PBFD F  F"),
    (23, "S27",  "PX634657", "23-s27  PBFD F  F"),
    (24, "Pb34", "PX634658", "24-pb34  PBFD F  F"),
    (25, "Pb64", "PX634659", "25-p64  PBFD F  F"),
    (26, "Pb23", "PX634660", "26-pb23  PBFD F  F"),
    (27, "Pb28", "PX634661", "27-p28  PBFD F  F"),
    (28, "Pb22", "PX634662", "28-pb22  PBFD F  F"),
]
# tip labels used in the 598-sequence PHYLIP alignment (10-char truncated)
PHY_TIP = {"19Pb24": "Pb24", "20Pb25_PBF": "Pb25", "21p30_PBFD": "Pb30", "22s7_PBFD_": "S7",
           "23s27_PBFD": "S27", "24pb34_PBF": "Pb34", "25p64_PBFD": "Pb64", "26pb23_PBF": "Pb23",
           "27p28_PBFD": "Pb28", "28pb22_PBF": "Pb22"}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].strip(); seqs[name] = []
        else:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(seqs, path):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")


# 1. Kerman unaligned sequences (GenBank-submitted versions)
raw = read_fasta(SRC / "10 seq.fas")
kerman = {}
for hdr, s in raw.items():
    acc = hdr.split("-")[0]
    voucher = next(v for _, v, a, _ in ISOLATES if a == acc)
    kerman[f"{acc}_Kerman_{voucher}"] = s.upper().replace("-", "")
write_fasta(kerman, DATA / "kerman_10_unaligned.fasta")

# 2. Raw (pre-edit) sequences keep lowercase = manually edited/low-quality calls
raw_edit = {}
for p in sorted((SRC / "DVM Daini" / "Seq RAW edition").glob("*.fsta")):
    s = "".join(l.strip() for l in open(p, encoding="utf-8", errors="replace") if not l.startswith(">"))
    raw_edit[p.stem] = s
write_fasta(raw_edit, DATA / "kerman_10_raw_basecalls.fasta")

# 3. 598-taxon thesis alignment (PHYLIP sequential, 10-char names)
lines = open(SRC / "new final PBFD.phy", encoding="utf-8").read().splitlines()
ntax, nchar = map(int, lines[0].split())
aln = collections.OrderedDict()
for l in lines[1:]:
    if not l.strip():
        continue
    name, seq = l[:10].strip(), l[10:].strip().replace(" ", "")
    aln[name] = seq.upper()
assert len(aln) == ntax, (len(aln), ntax)
assert all(len(s) == nchar for s in aln.values())
renamed = collections.OrderedDict()
for k, s in aln.items():
    if k in PHY_TIP:
        v = PHY_TIP[k]; acc = next(a for _, vv, a, _ in ISOLATES if vv == v)
        renamed[f"{acc}_Kerman_{v}"] = s
    else:
        renamed[k] = s
write_fasta(renamed, DATA / "thesis_alignment_598.fasta")

# 4. Reference metadata (host / source / country)
ws = openpyxl.load_workbook(SRC / "DVM Daini" / "extracted_features_final Darini.xlsx", read_only=True).active
refmeta = [dict(zip(["accession", "host", "isolation_source", "country"], r)) for r in list(ws.iter_rows(values_only=True))[1:] if r[0]]
for r in refmeta:
    r["country_simple"] = (r["country"] or "").split(":")[0].strip()
with open(META / "reference_metadata.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(refmeta[0])); w.writeheader(); w.writerows(refmeta)
in_aln = set(renamed) - set(kerman)
missing_meta = sorted(in_aln - {r["accession"] for r in refmeta})

# ---------------------------------------------------------------- 136-bird dataset
rows = list(openpyxl.load_workbook(SRC / "dataset_figures (1).xlsx", data_only=True).active.iter_rows(values_only=True))[1:]
EXCEL_TO_SUGGESTED = {
    # excel label : (suggested scientific name, common name, evidence, confidence)
    "Pyrrhura molinae": ("Pyrrhura molinae", "Green-cheeked conure", "Excel=thesis Table 4-1 (n=73, 48/13/11/1)", "high"),
    "Green-cheeked parakeet(Pyrrhura molinae)": ("Myiopsitta monachus", "Monk parakeet",
        "Excel label duplicates P. molinae, but thesis Table 4-1 and manuscript list M. monachus n=21 with identical band counts 19/1/1/0", "high"),
    "Lovebird(Agapornis spp.)": ("Forpus spp. OR Agapornis spp. [UNRESOLVED]", "Parrotlet OR lovebird",
        "Thesis Table 4-1 = Forpus spp. (n=20, 14/4/2/0); thesis isolate text calls Pb34 'کوتوله' (parrotlet); manuscript/Excel = Agapornis", "medium (favours Forpus)"),
    "African grey parrot(Psittacus erithacus)": ("Psittacus erithacus", "African grey parrot", "all sources agree", "high"),
    "Cockatiel(Nymphicus hollandicus)": ("Nymphicus hollandicus", "Cockatiel", "all sources agree", "high"),
    "Rose-ringed parakeet(Psittacula krameri)": ("Psittacula krameri", "Rose-ringed parakeet", "all sources agree", "high"),
    "Other spp.": ("Other spp. (Pionites/Pionus?, P. eupatria, Cacatua, Amazona, A. solstitialis)", "Other",
        "thesis lists Pionites spp.; manuscript lists Pionus maximiliani; individual species per bird not recorded", "low"),
}
BAND = {0: "None", 1: "Low", 2: "Moderate", 3: "High"}
species_rows = []
for code, name, sex, band in rows:
    lab = (name or "").strip()
    sug, common, ev, conf = EXCEL_TO_SUGGESTED[lab]
    species_rows.append(dict(sample=code, host_excel=lab, suggested_species=sug, suggested_common=common,
                             evidence=ev, confidence=conf, FINAL_SPECIES="", sex=sex, pcr_band=band,
                             band_label=BAND[band], positive=int(band > 0)))
with open(META / "species_mapping.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(species_rows[0])); w.writeheader(); w.writerows(species_rows)

# group-level placeholder table (used for figure labels before final names)
SHORT = {"Pyrrhura molinae": "P. molinae", "Green-cheeked parakeet(Pyrrhura molinae)": "M. monachus",
         "Lovebird(Agapornis spp.)": "Forpus/Agapornis spp.", "African grey parrot(Psittacus erithacus)": "P. erithacus",
         "Cockatiel(Nymphicus hollandicus)": "N. hollandicus", "Rose-ringed parakeet(Psittacula krameri)": "P. krameri",
         "Other spp.": "Other spp."}
groups = collections.OrderedDict()
for i, lab in enumerate(EXCEL_TO_SUGGESTED, 1):
    sug, common, ev, conf = EXCEL_TO_SUGGESTED[lab]
    groups[lab] = dict(group_id=f"SP{i:02d}", host_excel=lab, suggested_species=sug, suggested_short=SHORT[lab],
                       suggested_common=common, confidence=conf, FINAL_SPECIES="", FINAL_COMMON="")
with open(META / "species_groups.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(next(iter(groups.values())))); w.writeheader(); w.writerows(groups.values())

# ---------------------------------------------------------------- isolate host mapping
ws = openpyxl.load_workbook(SRC / "for correction.xlsx", data_only=True).active
host_excel = {str(r[3]).strip().lower(): r[1] for r in list(ws.iter_rows(values_only=True))[1:] if r[3]}
THESIS_HOST = {  # thesis pp.22-23, by sequencing number
    19: ("", "not described in thesis text; header in '10 seq.fas' says 'Missing'"),
    20: ("Nymphicus hollandicus", "thesis Fig 4-9 text: isolates 20,22,23,27 = عروس هلندی (cockatiel)"),
    21: ("Psittacus erithacus", "thesis Fig 4-9 text: isolate 21 = African grey"),
    22: ("Nymphicus hollandicus", "thesis Fig 4-9 text"),
    23: ("Nymphicus hollandicus", "thesis Fig 4-9 text"),
    24: ("Forpus sp.", "thesis Fig 4-9 text: isolate 24 = کوتوله (parrotlet)"),
    25: ("Pyrrhura molinae", "thesis Fig 4-10 text: isolate 25 = green-cheek"),
    26: ("Psittacus erithacus", "thesis Fig 4-10 text: isolate 26 = African grey"),
    27: ("Nymphicus hollandicus", "thesis Fig 4-9 text"),
    28: ("Psittacula krameri", "thesis Fig 4-8 text: isolate 28 = ملنگو (rose-ringed parakeet)"),
}
bird = {r["sample"].upper(): r for r in species_rows}
iso_rows = []
for no, v, acc, stem in ISOLATES:
    hx = host_excel.get(v.lower(), "")
    ht, tev = THESIS_HOST[no]
    s_match = bird.get(v.upper().replace("PB", "S")) if v.upper().startswith("S") else None
    s_note = ""
    if s_match:
        s_note = f"{s_match['sample']}: {s_match['host_excel']}, band {s_match['pcr_band']}"
    else:
        alt = bird.get("S" + re.sub(r"\D", "", v))
        if alt:
            s_note = f"no direct match (numeric twin {alt['sample']} = {alt['host_excel']}, band {alt['pcr_band']}; Pb series appears to be a separate numbering)"
    if ht and hx and ht.split()[0] == hx.split()[0]:
        sug, conf = ht, "high" if s_match else "medium"
        short = ht.split()[0][0] + ". " + ht.split()[1]
    elif ht and hx:
        sug, conf = f"{ht} (thesis) vs {hx} (GenBank sheet)", "conflict"
        short = f"{ht.split()[0]}/{hx.split()[0]}"
    else:
        sug, conf = hx or ht, "low"
        short = (sug.split()[0][0] + ". " + sug.split()[1]) if sug else "unknown"
    iso_rows.append(dict(seq_no=no, voucher=v, accession=acc, placeholder=f"HOST_{v}",
                         host_genbank_sheet=hx, host_thesis=ht, s_code_match=s_note,
                         suggested_host=sug, suggested_short=short, evidence=tev, confidence=conf, FINAL_HOST=""))
with open(META / "host_mapping.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(iso_rows[0])); w.writeheader(); w.writerows(iso_rows)

summary = dict(n_birds=len(species_rows), n_positive=sum(r["positive"] for r in species_rows),
               alignment_taxa=ntax, alignment_sites=nchar, kerman_in_alignment=sum(1 for k in renamed if "_Kerman_" in k),
               ref_metadata_rows=len(refmeta), refs_without_metadata=missing_meta)
json.dump(summary, open(META / "prep_summary.json", "w"), indent=2)
print(json.dumps(summary, indent=2))
