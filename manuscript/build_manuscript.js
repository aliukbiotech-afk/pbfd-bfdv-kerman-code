// Builds manuscript_PBFD_v2.docx  (node build_manuscript.js)
// Markup inside strings:  *italic*   **bold**   [[PENDING: ...]] / [[HOST_...]] -> yellow highlight
const fs = require("fs");
const path = require("path");
const D = require(process.env.DOCX_PATH || "docx");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, AlignmentType, HeadingLevel,
        BorderStyle, WidthType, ShadingType, LevelFormat, PageNumber, Footer, LineRuleType, ImageRun } = D;


// ------------------------------------------------------------------ embedded figures (downsized PNG copies of the final figures)
const FIG_DIR = process.env.FIG_DIR || path.join(__dirname, "docxfigs");
function figImg(n) {
  const buf = fs.readFileSync(path.join(FIG_DIR, `Fig${n}.png`));
  const w0 = buf.readUInt32BE(16), h0 = buf.readUInt32BE(20);
  const maxW = 600, maxH = 640;                      // px at 96 dpi (about 6.25 in x 6.7 in)
  const k = Math.min(maxW / w0, maxH / h0);
  return new Paragraph({ alignment: AlignmentType.CENTER, keepNext: true, pageBreakBefore: n > 1,
    spacing: { before: 120, after: 120 },
    children: [new ImageRun({ type: "png", data: buf, transformation: { width: Math.round(w0 * k), height: Math.round(h0 * k) },
      altText: { title: `Figure ${n}`, description: `Figure ${n}`, name: `Fig${n}` } })] });
}

// ------------------------------------------------------------------ helpers
const PENDING = [];
let curSection = "Title page";
const plain = t => t.split("*").join("").replace(/^PENDING:?\s*/, "").trim();
function runs(text, base = {}) {
  const out = [];
  const parts = text.split(/(\*\*[^*]+\*\*|\*[^*]+\*|\[\[[^\]]+\]\])/g).filter(s => s !== "");
  for (const s of parts) {
    if (s.startsWith("**")) out.push(new TextRun({ ...base, text: s.slice(2, -2), bold: true }));
    else if (s.startsWith("*")) out.push(new TextRun({ ...base, text: s.slice(1, -1), italics: true }));
    else if (s.startsWith("[[")) { PENDING.push({ section: curSection, text: plain(s.slice(2, -2)) }); out.push(new TextRun({ ...base, text: "...", highlight: "yellow" })); }
    else out.push(new TextRun({ ...base, text: s }));
  }
  return out;
}
const P = (t, o = {}) => new Paragraph({ children: runs(t, o.run || {}), alignment: o.align || AlignmentType.JUSTIFIED,
  spacing: { after: 120, line: 360, lineRule: LineRuleType.AUTO }, indent: o.indent, ...(o.extra || {}) });
const H1 = t => { curSection = plain(t); return new Paragraph({ heading: HeadingLevel.HEADING_1, children: runs(t) }); };
const H2 = t => { curSection = plain(t); return new Paragraph({ heading: HeadingLevel.HEADING_2, children: runs(t) }); };
const H3 = t => new Paragraph({ heading: HeadingLevel.HEADING_3, children: runs(t) });
const B = t => new Paragraph({ numbering: { reference: "bullets", level: 0 }, children: runs(t), spacing: { after: 80, line: 320 } });
const N = (t, ref = "nums") => new Paragraph({ numbering: { reference: ref, level: 0 }, children: runs(t), spacing: { after: 80, line: 320 } });

const border = { style: BorderStyle.SINGLE, size: 4, color: "999999" };
const borders = { top: border, bottom: border, left: border, right: border };
function table(header, rows, widths, caption, note) {
  const W = widths.reduce((a, b) => a + b, 0);
  const cell = (t, i, head) => new TableCell({
    borders, width: { size: widths[i], type: WidthType.DXA },
    shading: head ? { fill: "E7EEF5", type: ShadingType.CLEAR } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [new Paragraph({ children: runs(String(t), { size: 18, bold: head }), alignment: i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER })]
  });
  const out = [new Paragraph({ children: runs(caption, { size: 20 }), spacing: { before: 200, after: 80 }, keepNext: true }),
    new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: widths,
      rows: [new TableRow({ tableHeader: true, children: header.map((h, i) => cell(h, i, true)) }),
             ...rows.map(r => new TableRow({ children: r.map((c, i) => cell(c, i, false)) }))] })];
  if (note) out.push(new Paragraph({ children: runs(note, { size: 16 }), spacing: { before: 60, after: 240 } }));
  return out;
}

// ------------------------------------------------------------------ content
const C = [];
C.push(new Paragraph({ children: runs("**Subclinical beak and feather disease virus infection in pet psittacine birds in Kerman, south-eastern Iran: prevalence, risk factors and full-genome-guided genotyping of partial Rep sequences**", { size: 30 }), alignment: AlignmentType.CENTER, spacing: { after: 240 } }));
C.push(P("[[PENDING: author names, affiliations, ORCID, corresponding author e-mail]]", { align: AlignmentType.CENTER }));
C.push(P("**Running title:** Subclinical BFDV in pet psittacines, Kerman, Iran", { align: AlignmentType.CENTER }));
C.push(P("**Draft v3:** all analyses completed and verified. Yellow ellipses (...) mark items that still need information from the authors (listed at the end of this file). Host names follow the final list supplied by the authors.", { run: { size: 18, color: "555555" }, align: AlignmentType.CENTER }));

C.push(H1("Abstract"));
C.push(P("Beak and feather disease virus (BFDV; *Circoviridae*) is the most important viral pathogen of captive and wild parrots. Subclinically infected birds are thought to sustain its circulation through the pet trade. We screened feathers from 136 clinically healthy pet psittacines presented to veterinary clinics in Kerman, south-eastern Iran (April 2024–March 2025), by PCR targeting ORF1 (Rep). BFDV DNA was detected in 43 birds (31.6%; 95% CI 24.4–39.8%).\n".replace("\n", "") +
  " Positivity did not differ between males and females (30.9% vs 32.4%; Fisher's exact p = 1.00) or among host groups (Monte Carlo exact p = 0.29). In Firth-penalised logistic regression no host group or sex effect reached significance. Among positive birds, the semi-quantitative PCR band intensity differed among host groups (p = 0.007), mainly because both PCR-positive cockatiels showed strong amplification. Ten amplicons were sequenced (651 nt; GenBank PX634653–PX634662). They shared 94.9–99.8% nucleotide identity, encoded an intact Rep reading frame and conserved the rolling-circle (motifs II–III) and SF3-helicase (Walker A/B, motif C) motifs. Variation was concentrated at third codon positions (46/60 variable sites).".replace("\n", "")));
C.push(P("To overcome the limited resolution of partial sequences, we placed the Kerman sequences on a maximum-likelihood backbone of 429 complete BFDV genomes, classified according to the global GI/GII scheme. All ten sequences belonged to genotype GI. Eight fell in sub-clade GI-c (cumulative likelihood weight ratio, LWR, 0.90–1.00), one (Pb23) in GI-f (LWR 1.00), and one (Pb64) could not be resolved between GI-d (LWR 0.56) and GI-f (0.36). The GI-c sequences were 97.7–98.9% identical to strains from Italy, Saudi Arabia, Poland and other Iranian provinces, and the GI-f/GI-d pair to Iranian budgerigar strains from 2014 and a Bangladeshi strain. No recombination breakpoint was detected, and Rep was under strong purifying selection (dN/dS = 0.14; no site under positive selection). Subclinical BFDV infection is common in pet parrots in Kerman. Viruses circulating locally belong to internationally distributed lineages, which supports pre-sale screening and quarantine in the Iranian pet-bird trade."));
C.push(P("**Keywords:** beak and feather disease virus; *Circovirus*; psittacine birds; subclinical infection; replication-associated protein; phylogenetic placement; Iran"));

C.push(H1("1. Introduction"));
C.push(P("Psittacine beak and feather disease (PBFD) is a chronic, often fatal disease of parrots and cockatoos. It was first characterised pathologically in Australian cockatoos (Pass and Perry, 1984). Its agent, beak and feather disease virus (BFDV), is a member of the genus *Circovirus* (family *Circoviridae*). The virus has a circular, single-stranded DNA genome of about 2 kb with two major open reading frames: ORF1 (V1), encoding the replication-associated protein (Rep), and ORF2 (C1), encoding the capsid protein (Cap) (Bassami et al., 2001; Todd, 2000). BFDV targets rapidly dividing cells of the feather and beak epidermis and the lymphoid organs. The resulting feather dystrophy, beak deformities and immunosuppression predispose affected birds to secondary infections (Pass and Perry, 1984; Todd, 2000)."));
C.push(P("BFDV has been reported from captive and wild psittacines on every inhabited continent, and its range now extends beyond the native Australasian hosts (Fogell et al., 2016; Varsani et al., 2011). Phylogeographic analyses link this dispersal to the international trade in live birds (Harkins et al., 2014; Fogell et al., 2016). The virus is highly diverse and frequently recombines (Julian et al., 2013; Shah et al., 2023). Detection in non-psittacine birds raises concern about spill-over to wildlife (Amery-Gale et al., 2017). At the same time, endogenous circoviral elements in some avian genomes call for caution in interpreting PCR detections outside the Psittaciformes (Palacios-Martínez et al., 2026). Many infected parrots show no clinical signs yet shed virus in feather dust and droppings. Such subclinical carriers are therefore central to the persistence of BFDV in aviaries and pet shops (Hess et al., 2004; Fogell et al., 2016; Guan, 2026)."));
C.push(P("Using 454 complete genomes, Shah et al. (2023) proposed a two-clade framework (GI, sub-clades GI-a to GI-f; GII, sub-clades GII-a and GII-b) for global BFDV diversity. However, most surveillance studies, including all previous reports from Iran, have sequenced only short fragments of ORF1 or ORF2 generated by diagnostic PCR (Ypelaar et al., 1999; Dolatyabi et al., 2022). Trees built from such fragments often have weak support, and their placement within genome-based classifications is uncertain (Sarker et al., 2014). Phylogenetic placement of short query sequences on a fixed reference tree, widely used in microbial ecology (Barbera et al., 2019), offers a principled way to assign partial sequences to genome-defined lineages while quantifying the uncertainty of each assignment."));
C.push(P("Keeping pet parrots has become increasingly popular in Iran, and parrots are traded between provinces through breeders and pet shops. Molecular surveys from Tehran and several northern and central provinces have reported BFDV detection rates of about 18–58% (Haddadmarandi et al., 2018; Dolatyabi et al., 2022; Khosravi et al., 2024). Kerman province in the south-east has not been studied. Moreover, earlier Iranian surveys did not consistently distinguish clinically affected from clinically normal birds, and they analysed risk factors with tests poorly suited to sparse contingency tables."));
C.push(P("The present study had three aims:"));
C.push(N("to estimate the frequency of subclinical BFDV infection in pet psittacines in Kerman and its association with host group and sex, using exact and penalised methods appropriate for small samples;"));
C.push(N("to characterise the Rep coding region of local viruses; and"));
C.push(N("to assign the partial Kerman sequences to the global genome-based classification by placing them on a full-genome backbone."));

C.push(H1("2. Materials and methods"));
C.push(H2("2.1. Birds and sampling"));
C.push(P("Between April 2024 and March 2025, feathers were collected from 136 pet psittacine birds presented to veterinary clinics in Kerman city, Kerman province, Iran. Birds were mostly presented for reasons unrelated to PBFD, chiefly molecular sex determination. None showed feather or beak lesions or other signs suggestive of PBFD at examination. Four growing (blood) feathers were plucked from each bird into sterile tubes and stored at [[PENDING: storage temperature]] until DNA extraction. Host group and sex were recorded for each bird (68 males, 68 females) (Table 1). [[PENDING: how sex was determined (PCR sexing / owner record) and age data, if available]]."));
C.push(P("Host groups were *Pyrrhura molinae* (green-cheeked conure; n = 73), *Myiopsitta monachus* (monk parakeet; n = 21), *Agapornis* spp. (lovebirds; n = 20), *Psittacus erithacus* (grey parrot; n = 10), *Nymphicus hollandicus* (cockatiel; n = 4) and *Psittacula krameri* (rose-ringed parakeet; n = 3). The remaining five birds were single individuals of other species, recorded together as 'other' ([[PENDING: confirm — *Pionites* sp. or *Pionus* sp.; *Psittacula eupatria*; *Cacatua* sp.; *Amazona* sp.; *Aratinga solstitialis*]])."));
C.push(P("[[PENDING: ethics statement — institutional animal-ethics approval code, owner consent; feather sampling is minimally invasive]]."));

C.push(H2("2.2. DNA extraction and PCR"));
C.push(P("The calamus of each feather was cut into pieces of about 5 mm with a sterile blade. Genomic DNA was extracted with a commercial tissue DNA kit (Parstous, Mashhad, Iran) according to the manufacturer's instructions. A 717-bp fragment of ORF1 was amplified with the universal primers P2 (5′-AACCCTACAGACGGCGAG-3′) and P4 (5′-GTCACAGTCCTCCTTGTACC-3′) (Ypelaar et al., 1999) (Table 2). Reactions (25 µL) contained 4 µL template DNA, 15 µL PCR master mix ([[PENDING: PCR master mix brand]]), 1 µL of each primer ([[PENDING: primer concentration]]) and nuclease-free water."));
C.push(P("Cycling comprised initial denaturation at 95 °C for 4 min; 33 cycles of 95 °C for 30 s, 60 °C for 30 s and 72 °C for 90 s; and a final extension at 72 °C for 5 min. Each run included a no-template control and DNA from a confirmed BFDV-positive bird. Products were separated on 1% agarose gels at 120 V for 60 min and visualised under UV light. The intensity of the specific band was scored semi-quantitatively as none, low, moderate or high relative to the positive control ([[PENDING: confirm scoring criteria and scorer blinding]]). This score reflects amplicon yield and is not a measure of viral load."));

C.push(H2("2.3. Sequencing"));
C.push(P("Ten amplicons with strong bands were purified and sequenced bidirectionally/with the forward primer ([[PENDING: confirm — only forward (F) chromatograms are present in the data set]]) on an ABI platform ([[PENDING: company, city]]). Chromatograms were inspected and trimmed, and primer sequences were removed. Base-call quality was summarised as mean Phred score and the proportion of bases with Phred ≥ 20 (Ewing and Green, 1998) using Biopython (Cock et al., 2009). The 651-nt consensus sequences were deposited in GenBank under accession numbers PX634653–PX634662."));

C.push(H2("2.4. Reference data sets"));
C.push(P("Reference sequences comprised the 589 BFDV sequences used in the thesis tree, which were retrieved from GenBank (NCBI Entrez, queried without a registered contact e-mail address) together with collection date, country and host. Records ≥ 1,800 nt were treated as circular complete genomes. They were placed in the same orientation and rotated to start at the conserved origin nonanucleotide (TAGTATTAC). Shorter records were treated as linear partial sequences. The region homologous to the Kerman amplicon was extracted from every reference by local alignment."));
C.push(P("Sub-clade labels (GI-a to GI-f, GII-a, GII-b) for the 453 genomes in the global classification of Shah et al. (2023) were transcribed from the branch colours of their supplementary full-genome tree. Labels were then propagated to unlabelled genomes on our own full-genome tree: a genome inherited the label of the largest clade whose labelled members were ≥ 90% from a single sub-clade. This procedure recovered only GI-d, GI-e and GI-f as clean clades. Genomes outside these clades kept the transcribed label (GI-a, GI-b or GI-c) when one existed, and 86 of the 429 backbone genomes remained unassigned."));

C.push(H2("2.5. Alignment and phylogenetic analyses"));
C.push(P("Amplicon sequences were aligned with MAFFT v7 (Katoh and Standley, 2013) and trimmed to the window spanned by all Kerman sequences. Maximum-likelihood (ML) trees were inferred with IQ-TREE 3 (version 3.1.4; Minh et al., 2020). The substitution model was selected by ModelFinder under the Bayesian information criterion (Kalyaanamoorthy et al., 2017). Branch support was assessed with 1,000 ultrafast bootstrap replicates (UFBoot2; Hoang et al., 2018) and 1,000 SH-like approximate likelihood-ratio tests (Guindon et al., 2010)."));
C.push(P("For direct comparison with the thesis analysis and with Shah et al. (2023), a second tree was inferred under the TIM3+F+I+G4 model. The two topologies were compared by Robinson–Foulds distance. Trees were midpoint-rooted for display."));
C.push(P("A full-genome ML backbone was inferred in the same way from 429 complete genomes aligned with MAFFT (2,358 alignment columns; TIM3+F+I+R4 selected). The Kerman sequences were then assigned to sub-clades in three independent ways:"));
C.push(N("evolutionary placement with EPA-ng (Barbera et al., 2019) on the backbone topology, with GTR+G parameters re-estimated on the amplicon region with RAxML-NG v2.0.4 (Kozlov et al., 2019); assignments are reported as the cumulative likelihood weight ratio (LWR) per sub-clade;", "nums2"));
C.push(N("an IQ-TREE amplicon tree constrained to the backbone topology; and", "nums2"));
C.push(N("the unconstrained amplicon tree.", "nums2"));

C.push(H2("2.6. Genetic distances, diversity and haplotypes"));
C.push(P("Uncorrected pairwise p-distances were computed with pairwise deletion of gaps. For each Kerman sequence, the ten closest references and the mean identity to each sub-clade were tabulated. Nucleotide diversity (π), haplotype diversity (Hd) and the number of segregating sites were calculated for the Kerman sequences and for all Iranian sequences. A minimum-spanning haplotype network (Bandelt et al., 1999) was built for the Kerman sequences and their nearest references."));

C.push(H2("2.7. Rep motifs, recombination and selection"));
C.push(P("Translated sequences were screened for the conserved rolling-circle replication motifs and the superfamily 3 helicase motifs of circovirus Rep. A codon alignment was built for a subset comprising the Kerman sequences, all Iranian references and the five closest references to each Kerman sequence: proteins were aligned with MAFFT L-INS-i and the nucleotides threaded back onto the protein alignment."));
C.push(P("Recombination was screened with GARD (Kosakovsky Pond et al., 2006). Selection was analysed with SLAC and FEL (Kosakovsky Pond and Frost, 2005) and MEME (Murrell et al., 2012) in HyPhy 2.5.103 (Kosakovsky Pond et al., 2020), with significance set at p ≤ 0.1."));

C.push(H2("2.8. Phylogeny–trait association"));
C.push(P("Clustering of sequences by country and by host genus on the ML amplicon tree was tested with the parsimony score, the association index (Wang et al., 2001) and the maximum monophyletic clade size. Each statistic was compared with a null distribution from 1,000 random permutations of tip labels (Parker et al., 2008)."));

C.push(H2("2.9. Statistical analysis of survey data"));
C.push(P("Detection frequencies are reported with Wilson score 95% confidence intervals (Wilson, 1927). Association between positivity and sex was tested with Fisher's exact test. Association between positivity and host group was tested with a Monte Carlo exact test of independence: 100,000 permutations with both margins fixed, with the Pearson χ² as test statistic. The test was run with all host groups and again with groups of fewer than 10 birds combined. Effect sizes are given as Cramér's V."));
C.push(P("Host group and sex were analysed jointly with Firth-penalised logistic regression (Heinze and Schemper, 2002), using *P. molinae* and females as reference categories. Because the band-intensity score is defined only for positive birds, intensity was analysed separately among the 43 positives. We used the Monte Carlo exact test for host group and sex, the Kruskal–Wallis test (host group) and the Mann–Whitney U test (sex)."));
C.push(P("To explain the discrepancy with the original thesis analysis, we also reproduced the host-group × four-level band table (none/low/moderate/high) analysed there with an asymptotic χ² test. All analyses used Python 3 (NumPy, SciPy, statsmodels). Code, software versions and random seeds are provided in the Supplementary Material. α = 0.05."));

C.push(H1("3. Results"));
C.push(H2("3.1. Detection of BFDV DNA in clinically healthy birds"));
C.push(P("BFDV DNA was detected in 43 of 136 birds (31.6%; 95% CI 24.4–39.8%). Positivity was similar in females (22/68, 32.4%) and males (21/68, 30.9%) (odds ratio 0.93; Fisher's exact p = 1.00) (Figure 1B). Detection frequency ranged from 9.5% in *M. monachus* (2/21; 95% CI 2.7–28.9%) to 50% in *P. erithacus* (5/10) and *N. hollandicus* (2/4) (Table 1; Figure 1A). Most groups were small and their confidence intervals overlapped widely."));
C.push(P("Positivity was not associated with host group, either with all groups (χ² = 7.35, Monte Carlo p = 0.29, Cramér's V = 0.23) or with groups of fewer than 10 birds combined (χ² = 7.12, p = 0.13) (Table 3). In the penalised logistic model, no host group differed significantly from *P. molinae*. *M. monachus* showed the lowest odds of positivity (OR 0.25; 95% CI 0.06–1.03; p = 0.055), and sex had no effect (OR male vs female 0.89; 95% CI 0.42–1.88; p = 0.76) (Table S2)."));

C.push(H2("3.2. Semi-quantitative band intensity among positive birds"));
C.push(P("Of the 43 positive birds, 24 (55.8%) showed low, 16 (37.2%) moderate and 3 (7.0%) high band intensity (Table S2). Intensity was unrelated to sex (Monte Carlo p = 0.35; Mann–Whitney p = 0.10). It differed among host groups (χ² = 34.05, Monte Carlo p = 0.007, Cramér's V = 0.63). This was driven largely by the two positive cockatiels, both of which showed high intensity, while the rank-based Kruskal–Wallis test was not significant (p = 0.07). The data should therefore be read as a signal for further quantitative work rather than as evidence of host differences in viral load."));
C.push(P("Analysed as in the original thesis, the host-group × four-level band table gave χ² = 59.2 (asymptotic p < 0.001). However, 79% of its cells had expected counts below 5 (minimum 0.07), and the Monte Carlo p-value was 0.002. Because 'none' is one of the four categories, this test mixes infection status with intensity; it should not be interpreted as an effect of host on the severity of infection."));

C.push(H2("3.3. Sequence quality and Rep characteristics"));
C.push(P("Mean Phred scores of the trimmed forward reads ranged from 38.5 to 50.2, and 74.7–94.7% of bases had Phred ≥ 20 (Table S1; Figure S2). Eight of the ten reads covered ≥ 94% of the submitted sequence, and nine matched it at ≥ 97.5% identity. The exceptions were S27 (read covering 71.7% of the submitted sequence), Pb64 (82.9%) and S7, whose read matched the submitted sequence at only 90.4% identity [[PENDING: authors to confirm how the S7, S27 and Pb64 consensus sequences were obtained, e.g. a reverse read or manual editing]]."));
C.push(P("All ten sequences were 651 nt long, without ambiguous bases except one in PX634661 (Pb28). Each translated into an uninterrupted 217-amino-acid Rep fragment. Pairwise identity among the Kerman sequences was 94.9–99.8% at the nucleotide level (mean 96.8%) and 96.3–100% at the amino-acid level. Of 60 variable nucleotide sites, 46 (76.7%) were at third codon positions, and only 16 amino-acid positions varied."));
C.push(P("The rolling-circle motifs II (HLQG) and III (YCSK), the SF3 helicase Walker A (GPPGCGKS) and motif C (ITSN) were fully conserved in all ten sequences. Walker B was VVLDD in nine sequences and VILDD in PX634660 (Pb23) (Table 4). On identity, the sequences formed three groups: seven sequences (Pb24, Pb25, Pb28, Pb30, Pb34, S7, S27; 96.8–99.8%), a pair (Pb23–Pb64; 98.9%), and a single divergent sequence (Pb22; ≤ 96.3% to all others)."));

C.push(H2("3.4. Phylogenetic relationships and genotype assignment"));
C.push(P("The amplicon alignment contained 599 sequences and 651 sites; ModelFinder selected TIM2+I+R4 (log-likelihood −16,659). The tree inferred under TIM3+F+I+G4 had a lower likelihood (−16,725) and differed from it by a Robinson–Foulds distance of 308 (maximum 1,192), reflecting the short, weakly supported internal branches typical of a 651-nt fragment. Seven Kerman sequences (Pb24, Pb25, Pb28, Pb30, Pb34, S7, S27) formed a group of 11 tips with Saudi Arabian (MK803398, MK803403), Italian and Iranian (OP039998) strains, but this node was weakly supported (SH-aLRT 78.5, UFBoot 25). Within it, S7 and S27 (SH-aLRT/UFBoot 94.1/96) and Pb30 and Pb34 (96.5/100) formed well-supported pairs. Pb22 grouped with Iranian and Pakistani rose-ringed parakeet and budgerigar strains (KT764931, KT764932, OP039999, KX641239). Pb23 and Pb64 formed a pair (97.1/99) that was sister to MK387075 (Bangladesh) and nested among Iranian budgerigar strains from 2014 (KT764926–28, KT764934) (Figure 3; Figure S1)."));
C.push(P("The full-genome backbone contained 429 genomes (Figure S4). EPA-ng placed all ten Kerman sequences in genotype GI. Eight were assigned to sub-clade GI-c with cumulative LWR of 0.90–1.00 (Pb24, Pb30, Pb34, S7, S27: ≥ 0.999; Pb22: 0.99; Pb25: 0.90; Pb28: 0.91). Pb23 was placed in GI-f (LWR 1.00). Pb64 was split between GI-d (0.56), GI-f (0.36) and GI-e (0.02) and is therefore reported as unresolved GI-d/GI-f (Figure 4, Table 5). The constrained tree also placed Pb23 and Pb64 in GI-f. For the other eight sequences neither amplicon tree contained a clade with a uniform sub-clade label, so no label could be transferred. Placement of Pb22 rests on an edge with 105 backbone tips, of which 30 are GI-c and 75 are unlabelled, so its assignment to GI-c should be regarded as provisional."));
C.push(P("The nearest references of the GI-c sequences were KF723390 (Italy, grey parrot, 2011; Pb24, 98.9%), MK803398 and MK803403 (Saudi Arabia, grey parrot, 2017; 98.0–98.5%), JX221037 (Poland) and OP039998 (Iran, rose-ringed parakeet, 2019; S7, 98.2%). Pb22 was closest to OP039999 (Iran, budgerigar, 2019; 97.8%). Pb23 was 99.1% identical to MK387075 (Bangladesh, budgerigar, 2017) and 98.9% to the Iranian budgerigar strain KT764927; Pb64 was 98.2% identical to KT764927. The mean identity of each GI-c sequence was highest to GI-c (95.3–95.9%), whereas Pb23 and Pb64 were most similar to GI-e (96.9–97.0%), then to GI-d and GI-f, which illustrates why identity alone should not be used to name a sub-clade (Table S3). The ten Kerman sequences represented ten distinct haplotypes (Hd = 1.00; π = 0.031; 60 segregating sites). Thirty-one Iranian sequences (including ours) formed 29 haplotypes (Hd = 0.996; π = 0.042), and the 599 sequences formed 444 haplotypes (π = 0.055). Kerman haplotypes were scattered across the minimum-spanning network without a single Kerman cluster (Figure 5)."));
C.push(H2("3.5. Recombination, selection and phylogeny–trait association"));
C.push(P("GARD, run on the 41-sequence subset (Kerman and Iranian sequences plus their nearest references), found no evidence of a recombination breakpoint in the 651-nt fragment (101 models examined; no improvement in AICc). The global dN/dS estimated by SLAC was 0.14. FEL detected 65 of 217 codons under purifying selection (p ≤ 0.1) and none under diversifying selection, and MEME found no site with episodic positive selection (Figure S3). The phylogeny was strongly associated with country (association index 3.3 vs 52.4 expected by chance; parsimony score 85 vs 404; p = 0.001, the smallest value possible with 1,000 permutations) and with host genus among reference sequences (association index 15.3 vs 54.3; parsimony score 141 vs 438; p = 0.001; Table S4). Both associations were dominated by clusters of sequences submitted from single studies and should not be read as evidence of host restriction."));

C.push(H1("4. Discussion"));
C.push(P("To our knowledge, this is the first molecular survey of BFDV in pet parrots in south-eastern Iran and the first Iranian study restricted to birds without clinical signs. Almost one in three apparently healthy parrots carried detectable BFDV DNA in their feathers. This figure lies within the range reported for Iranian aviary and pet-shop populations (Haddadmarandi et al., 2018; Dolatyabi et al., 2022; Khosravi et al., 2024) and within the broad range of 3.5% to 45% reported elsewhere (Rahaus and Wolff, 2003; Hsu et al., 2006; Hakimuddin et al., 2016; Ko et al., 2024; Zhang et al., 2025). Direct comparison is limited by differences in sampling frame (clinically normal versus mixed or sick birds), sample type, PCR target and species composition."));
C.push(P("Detecting viral DNA in feathers of birds without lesions indicates infection or contamination with shed virus. It does not establish productive infection or predict progression to disease (Hess et al., 2004). Nevertheless, such birds are plausible sources of transmission in clinics, shops and breeding facilities."));
C.push(P("Neither sex nor host group was associated with positivity. The absence of a sex effect is consistent with the Iranian survey of Khosravi et al. (2024) and with a recent global meta-analysis (Zhang et al., 2025). The low frequency in monk parakeets (9.5%) contrasts with the high prevalence reported in invasive monk and rose-ringed parakeets in Spain (Morinha et al., 2020) [[PENDING: verify figures cited in the draft (33% and 37%) against the source]]. This may reflect differences in population origin and density, but our sample was small."));
C.push(P("The original thesis reported a highly significant association between host and 'severity' (P < 0.001). We show that this result arose from a sparse four-level table that combined infection status with band intensity. Analysed appropriately, host group was not associated with infection, and evidence for host differences in band intensity among positives was weak and driven by two cockatiels. PCR band intensity is at best a crude proxy for viral load. Quantitative PCR would be needed to test whether cockatiels or grey parrots shed more virus, as suggested for grey parrots elsewhere (Amery-Gale et al., 2017; Dolatyabi et al., 2022)."));
C.push(P("The Kerman sequences preserved all Rep motifs required for rolling-circle replication and ATP-dependent helicase activity. Most variation was synonymous, at third codon positions. This pattern agrees with the strong purifying selection generally reported for BFDV Rep, which is more conserved than Cap (Shah et al., 2023; Ko et al., 2024) In our data a global dN/dS of 0.14, 65 purifying sites and no positively selected site point in the same direction. The finding that all ten sequences differ from those already in GenBank reflects the high genetic diversity of BFDV. It should not be taken as evidence of unusual mutation rates in Kerman."));
C.push(P("Phylogenetic placement assigned all ten Kerman sequences to genotype GI. Eight fell in GI-c, a sub-clade that contains viruses from grey parrots and other hosts in Europe and the Middle East, and these were closest to strains from Italy, Saudi Arabia, Poland and elsewhere in Iran. Two sequences (Pb23, Pb64) belonged to GI-f or GI-d and were almost identical to budgerigar strains previously reported from Iran (2014) and Bangladesh. A mixture of unrelated lineages in one city, with several close matches to foreign strains, is what would be expected if local viruses were introduced repeatedly through the live-bird trade (Harkins et al., 2014; Fogell et al., 2016). The data do not allow us to date or direct these introductions, and the strong association between phylogeny and country in the reference set partly reflects clusters submitted by single studies. Placement on the backbone gave a probability-based assignment for fragments that the amplicon tree alone could not classify, and for Pb64 it showed that the assignment is ambiguous. The approach is limited by the quality of the sub-clade labels on the backbone (86 of 429 genomes unlabelled), by the short fragment, and by recombination elsewhere in the genome, which GARD cannot exclude for the 651-nt region studied here (Julian et al., 2013; Shah et al., 2023)."));
C.push(H2("4.1. Limitations"));
C.push(P("The study has several limitations:"));
C.push(B("Birds were a convenience sample of clinic referrals, mostly for sex determination, and several host groups were represented by fewer than ten birds. This limits generalisation and statistical power."));
C.push(B("A single conventional PCR on feathers may miss low-level infections, so the true frequency may be underestimated. Band intensity is not a validated quantitative measure."));
C.push(B("No follow-up was possible to determine whether positive birds developed disease or cleared infection."));
C.push(B("Only 651 nt of ORF1 (about one-third of the genome) were sequenced. Placement on a full-genome backbone mitigates but does not remove this limitation, particularly where recombination has created mosaic genomes."));
C.push(B("Host species were recorded by the owners or clinic staff and were not confirmed by molecular identification. The voucher numbers of the sequenced birds (Pb series) cannot be matched individually to the 136 surveyed birds [[PENDING: authors to confirm that all sequenced birds belong to the 136 surveyed birds; five sequenced birds are cockatiels, whereas the survey table has only four cockatiels, two of them positive]]."));
C.push(P("Complete-genome sequencing of local strains and longitudinal sampling of breeding facilities would address these gaps."));
C.push(H1("5. Conclusions"));
C.push(P("Subclinical BFDV infection is frequent among pet parrots in Kerman and is not confined to particular host groups or to either sex. The local viruses carry a conserved Rep and belong mainly to sub-clade GI-c (eight of ten sequences) and, in two cases, to GI-f or GI-d of the global BFDV phylogeny, together with strains circulating in Italy, Saudi Arabia, Poland, Bangladesh and elsewhere in Iran. Routine PCR screening of birds entering shops and breeding facilities, quarantine of new birds and improved traceability of the domestic and international trade are recommended. Phylogenetic placement on a full-genome backbone offers a practical way to genotype the short diagnostic sequences that dominate BFDV surveillance."));

C.push(H1("Declarations"));
C.push(P("**Ethics statement:** [[PENDING: ethics approval code and statement, as in Section 2.1]]. **Funding:** [[PENDING: funding sources]]. **Conflict of interest:** The authors declare no conflict of interest. **Author contributions:** [[PENDING: CRediT]]. **Data availability:** Sequences are available in GenBank (PX634653–PX634662). Survey data, alignments, trees, the Colab notebook and the analysis code are provided in a public code repository (https://github.com/aliukbiotech-afk/pbfd-bfdv-kerman-code; archived release DOI [[PENDING: Zenodo DOI; make the repository public before submission]]) and as Supplementary Material."));

// ------------------------------------------------------------------ tables
C.push(H1("Tables"));
C.push(...table(["Host group", "n", "Positive (%)", "95% CI", "Low", "Moderate", "High"], [
  ["*Pyrrhura molinae*", 73, "25 (34.2)", "24.4–45.7", 13, 11, 1],
  ["*Myiopsitta monachus*", 21, "2 (9.5)", "2.7–28.9", 1, 1, 0],
  ["*Agapornis* spp.", 20, "6 (30.0)", "14.5–51.9", 4, 2, 0],
  ["*Psittacus erithacus*", 10, "5 (50.0)", "23.7–76.3", 5, 0, 0],
  ["*Nymphicus hollandicus*", 4, "2 (50.0)", "15.0–85.0", 0, 0, 2],
  ["*Psittacula krameri*", 3, "1 (33.3)", "6.1–79.2", 0, 1, 0],
  ["Other species", 5, "2 (40.0)", "11.8–76.9", 1, 1, 0],
  ["**Total**", "**136**", "**43 (31.6)**", "**24.4–39.8**", "**24**", "**16**", "**3**"]],
  [3000, 700, 1300, 1300, 900, 1200, 960],
  "**Table 1.** BFDV PCR results and semi-quantitative band intensity by host group (clinically healthy pet psittacines, Kerman, 2024–2025).",
  "CI, Wilson score interval. Low/Moderate/High, band intensity among positives."));
C.push(...table(["Primer", "Sequence (5′–3′)", "Position in ORF1", "Product"], [
  ["P2 (forward)", "AACCCTACAGACGGCGAG", "182–199", "717 bp"],
  ["P4 (reverse)", "GTCACAGTCCTCCTTGTACC", "879–898", ""]],
  [1800, 3500, 2200, 1860], "**Table 2.** Primers used for BFDV detection (Ypelaar et al., 1999).", null));
C.push(...table(["Comparison", "Test", "Statistic", "p", "Effect size"], [
  ["Positivity × sex", "Fisher's exact", "OR = 0.93", "1.00", "—"],
  ["Positivity × host group (7 groups)", "Monte Carlo exact (10⁵)", "χ² = 7.35", "0.29", "V = 0.23"],
  ["Positivity × host group (n < 10 pooled)", "Monte Carlo exact", "χ² = 7.12", "0.13", "V = 0.23"],
  ["Intensity × host group (positives, n = 43)", "Monte Carlo exact", "χ² = 34.05", "0.007", "V = 0.63"],
  ["Intensity × host group (positives)", "Kruskal–Wallis", "H = 11.60", "0.07", "—"],
  ["Intensity × sex (positives)", "Monte Carlo exact / Mann–Whitney", "χ² = 2.81 / U = 290.5", "0.35 / 0.10", "V = 0.26"],
  ["Thesis replication: host × 4-level band", "χ² asymptotic / Monte Carlo", "χ² = 59.2", "< 0.001 / 0.002", "79% cells E < 5"]],
  [3000, 2100, 1700, 1300, 1260], "**Table 3.** Association tests for BFDV positivity and band intensity.",
  "V, Cramér's V; E, expected count. Firth logistic regression results are given in Table S2."));
C.push(...table(["Motif", "Consensus (this study)", "Position (aa in fragment)", "Variation"], [
  ["RCR motif II", "HLQG", "17–20", "none"],
  ["RCR motif III", "YCSK", "55–58", "none"],
  ["SF3 Walker A", "GPPGCGKS", "131–138", "none"],
  ["SF3 Walker B", "VVLDD", "166–170", "VILDD in PX634660 (Pb23)"],
  ["SF3 motif C", "ITSN", "207–210", "none"]],
  [2200, 2300, 2300, 2560], "**Table 4.** Conserved Rep motifs in the ten Kerman sequences.",
  "Positions refer to the 217-aa translated amplicon."));
C.push(...table(["Isolate (accession)", "Host", "EPA-ng sub-clade (LWR)", "Constrained tree", "Nearest reference (% identity)"], [
  ["Pb24 (PX634653)", "*N. hollandicus*", "GI-c (1.00)", "no label", "KF723390, Italy (98.9)"],
  ["Pb25 (PX634654)", "*N. hollandicus*", "GI-c (0.90)", "no label", "Table S3"],
  ["Pb30 (PX634655)", "*P. erithacus*", "GI-c (1.00)", "no label", "Table S3"],
  ["S7 (PX634656)", "*N. hollandicus*", "GI-c (1.00)", "no label", "OP039998, Iran (98.2)"],
  ["S27 (PX634657)", "*N. hollandicus*", "GI-c (1.00)", "no label", "Table S3"],
  ["Pb34 (PX634658)", "*Agapornis* spp.", "GI-c (1.00)", "no label", "Table S3"],
  ["Pb64 (PX634659)", "*P. molinae*", "GI-d (0.56) / GI-f (0.36)", "GI-f", "KT764927, Iran (98.2)"],
  ["Pb23 (PX634660)", "*P. erithacus*", "GI-f (1.00)", "GI-f", "MK387075, Bangladesh (99.1)"],
  ["Pb28 (PX634661)", "*N. hollandicus*", "GI-c (0.91)", "no label", "Table S3"],
  ["Pb22 (PX634662)", "*P. krameri*", "GI-c (0.99)", "no label", "OP039999, Iran (97.8)"]],
  [1700, 1700, 2100, 1500, 2360], "**Table 5.** Sub-clade assignment of the ten Kerman sequences.",
  "LWR, cumulative likelihood weight ratio of EPA-ng placements on the 429-genome backbone. No label: the constrained amplicon tree had no uniformly labelled clade around the sequence. Pb64 is unresolved between GI-d and GI-f."));

C.push(H1("Figure legends"));
C.push(figImg(1));
C.push(P("**Figure 1.** BFDV detection in clinically healthy pet psittacines, Kerman. (A) Percentage of PCR-positive birds by host group with Wilson 95% confidence intervals (n = birds tested; dashed line, overall frequency of 31.6%). (B) Percentage positive by sex, with positive/total in the bars. Host-group colours are the same in Figures 1 and 2 and Figure S2."));
C.push(figImg(2));
C.push(P("**Figure 2.** Who was tested and who was positive. All 136 birds are shown in three columns of equal height: host group (with positive/total), PCR result and sex. PCR-positive birds (43) occur in every host group and in both sexes (22 females, 21 males), and positivity was not associated with host group (Monte Carlo p = 0.29) or sex (Fisher p = 1.00). PCR band intensity of the positive birds is not shown here and is reported in the text and Table S2."));
C.push(figImg(3));
C.push(P("**Figure 3.** Maximum-likelihood tree of the ORF1 amplicon showing the Kerman sequences (orange), other Iranian strains (pink), their nearest references and representatives of each sub-clade. Tip labels give host species, country and year (isolate and accession numbers are in Tables S6 and S7). The coloured bar shows the sub-clade (Shah et al., 2023); its colours are identical in Figures 3 and 4 and Figures S1 and S4. Dots mark nodes with UFBoot ≥ 95% and SH-aLRT ≥ 80%. Model TIM2+I+R4; 599 sequences × 651 sites."));
C.push(figImg(4));
C.push(P("**Figure 4.** Genotype assignment of the ten Kerman Rep sequences by phylogenetic placement. *What is shown:* each bar is one Kerman sequence, placed with EPA-ng on a reference tree of 429 complete BFDV genomes whose sub-clades GI-a to GI-f follow Shah et al. (2023). Segments give the likelihood weight ratio (LWR), i.e. the probability that the sequence belongs to each sub-clade; hatched, weight on branches without a sub-clade label. *Why it is needed:* a 651-nt fragment cannot be classified reliably by sequence identity alone, so placement on a whole-genome tree is used both to name the sub-clade and to measure how certain that name is. *Take-home message:* eight sequences fall in GI-c (LWR 0.90–1.00), one (Pb23) in GI-f (1.00), and one (Pb64) is split between GI-d (0.56) and GI-f (0.36) and is therefore left unresolved. Bars are labelled by host species and ordered, from top to bottom, S27, Pb34, Pb24, Pb30, S7, Pb22, Pb28, Pb25, Pb64, Pb23 (Table S6)."));
C.push(figImg(5));
C.push(P("**Figure 5.** (A) Pairwise nucleotide identity among the Kerman sequences (orange, *) and other Iranian sequences, labelled by host species (genus abbreviated). (B) Minimum-spanning haplotype network of the same sequences and their nearest non-Iranian references. Each circle is one haplotype and numbers on links are nucleotide substitutions; colours show origin (orange, this study; pink, other Iranian strains; grey, nearest non-Iranian strains). Isolate and accession numbers are in Tables S6 and S7."));
C.push(P("**Supplementary figures.** **Figure S1.** Complete maximum-likelihood tree of the 599 amplicon sequences (host species, country, year; sub-clade bar as in Figure 3). **Figure S2.** Quality of the forward Sanger reads, ranked from best to worst: bars give the percentage of trimmed bases with Phred ≥ 20 and the mean Phred score; colours are the host groups of Figure 1. The two lowest-ranking reads (both cockatiels, S7 and S27) are among those with incomplete coverage or reduced identity to the submitted sequence (Table S1). **Figure S3.** Selection and conserved motifs in the 217-codon Rep fragment. (A) Site-wise significance (−log10 p) of FEL and MEME; all 65 codons with FEL p < 0.1 are under purifying selection and no codon shows positive or episodic selection; blue bands mark the two rolling-circle-replication motifs (II, III) and the three SF3-helicase motifs (Walker A, Walker B, C). (B) Residues at the five motifs in the ten Kerman isolates compared with their consensus: all isolates are identical at every motif except one (Pb23, Walker B VILDD instead of VVLDD, highlighted). Rows are labelled by host species and ordered Pb34; Pb24, Pb25, Pb28, S27, S7; Pb22; Pb23, Pb30; Pb64. **Figure S4.** Maximum-likelihood tree of the 429 full genomes used as the placement backbone (host species, country, year; sub-clade bar as in Figure 3). Tables S1–S7 and the analysis notebook are provided as supplementary files."));

// ------------------------------------------------------------------ references
C.push(H1("References"));
const refs = JSON.parse(fs.readFileSync(path.join(__dirname, "references.json"), "utf8"));
for (const r of refs) C.push(new Paragraph({ children: runs(r), alignment: AlignmentType.LEFT, indent: { left: 360, hanging: 360 }, spacing: { after: 80, line: 300 } }));

// ------------------------------------------------------------------ pending corrections
C.push(H1("Appendix (internal; remove before submission): items to complete"));
C.push(P("Each yellow '...' in the text corresponds to one line below, in order of appearance."));
PENDING.forEach((it, i) => C.push(P(`${i + 1}. **${it.section}:** ${it.text}`, { align: AlignmentType.LEFT })));
C.push(H1("Appendix (internal; remove before submission): notes on corrections"));
for (const t of [
  "**Host names (resolved).** Final names were supplied by the authors: cockatiel *Nymphicus hollandicus* (Pb24, Pb25, S7, S27, Pb28), African grey parrot *Psittacus erithacus* (Pb30, Pb23), lovebird *Agapornis* spp. (Pb34), green-cheeked parakeet *Pyrrhura molinae* (Pb64) and rose-ringed parakeet *Psittacula krameri* (Pb22). These agree with the GenBank sheet. The thesis called Pb34 a parrotlet (*Forpus*); the *Agapornis* label has been adopted for the isolate and, by inference, for the 20-bird survey group (thesis Table 4-1 says *Forpus*), so please confirm the survey group. The five 'other' birds remain unspecified.",
  "**Monk parakeet label.** The 21 rows labelled 'Green-cheeked parakeet (Pyrrhura molinae)' in dataset_figures.xlsx correspond to *Myiopsitta monachus* in the thesis (identical band counts) and have been treated as such. Please confirm.",
  "**Sequenced birds vs survey.** Five sequenced isolates are cockatiels (Pb24, Pb25, S7, S27, Pb28), but the survey contains only four cockatiels, two of them PCR-positive (S7, S27). Pb vouchers do not match the S sample codes. Please confirm that all sequenced birds belong to the 136 surveyed birds.",
  "**Sequencing.** Only forward chromatograms exist. For S27 the read covers 72% of the submitted sequence, for Pb64 83%, and for S7 read and sequence agree at 90%. Please confirm how the consensus sequences were obtained (reverse reads? re-sequencing?).",
  "**Collection dates.** The GenBank submission table lists collection dates 1–10 May 2024, one per isolate per day; confirm real dates.",
  "**Literature figures.** Prevalence figures cited from the literature were carried over from the first draft and must be checked against the original papers: Dolatyabi 58%; Khosravi 56.2%, 65.7%, 50.5%; Haddadmarandi 18.2%; Morinha 33%/37%; country prevalences.",
  "**Draft text removed or corrected.** The draft said 'nine' sequenced samples (there are ten; Pb24 was missing from the thesis tree). 'Viral load' has been replaced by 'band intensity'. Claims that partial-gene and whole-genome analyses are consistent, that all sequences are novel and therefore indicate significant mutation, and that mutations 'originated from Africa, Asia and Eurasia' have been removed. Duplicate references (Hsu 2006a/b; Ko 2024a/b) were merged; 'Pilar et al.' was corrected to Ypelaar et al. (1999); the undated PubMed entry was replaced by Haddadmarandi et al. (2018). Off-topic references (Kim et al. 2022, polyomavirus; Todd et al. 2008, pigeon circovirus; Speer and van Zeeland 2017, a book review) were dropped.",
  "**Statistics.** SPSS χ² tests were replaced by exact, Monte Carlo and penalised methods (Section 2.9); all numbers in Section 3.1–3.2 come from analysis/results/run1/stats (Colab run, independently recomputed; see VERIFICATION.md)."
]) C.push(B(t));

// ------------------------------------------------------------------ document
const doc = new Document({
  creator: "PBFD Kerman project", title: "PBFD Kerman manuscript v2",
  styles: {
    default: { document: { run: { font: "Times New Roman", size: 24 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: "Times New Roman" }, paragraph: { spacing: { before: 360, after: 160 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, font: "Times New Roman" }, paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, italics: true, font: "Times New Roman" }, paragraph: { spacing: { before: 200, after: 100 }, outlineLevel: 2 } }]
  },
  numbering: { config: [
    { reference: "bullets", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] },
    { reference: "nums", levels: [{ level: 0, format: LevelFormat.LOWER_LETTER, text: "(%1)", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] },
    { reference: "nums2", levels: [{ level: 0, format: LevelFormat.LOWER_LETTER, text: "(%1)", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] }] },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } }, lineNumbers: { countBy: 1, restart: "newPage" } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: [PageNumber.CURRENT] })] })] }) },
    children: C
  }]
});
Packer.toBuffer(doc).then(b => { const out = path.join(__dirname, "..", "..", "manuscript_PBFD_v2.docx"); fs.writeFileSync(out, b); console.log("wrote", out); });
