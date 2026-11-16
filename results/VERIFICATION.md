# VERIFICATION of the Colab run (session pbfd4, 21/21 cells, status DONE)

Verified locally on 2026-11-16 from `results.zip` (unzipped to `run1/`) and `cells_logs.zip`.
Legend: PASS = checked and consistent; CAVEAT = usable but must be worded carefully; OPEN = needs user input.

## 1. Run integrity
| Item | Result |
|---|---|
| status.json | PASS: state DONE, failed_cell null, cells 1-21 |
| Counts | PASS: 136 birds, 43 positive, 10 isolates (PX634653-PX634662) |
| Tool versions (versions.txt) | PASS: IQ-TREE 3.1.4, MAFFT 7.526, RAxML-NG 2.0.4, EPA-ng 0.3.8, HyPhy 2.5.103, Python 3.13.18; seed 20261103; UFBoot/SH-aLRT 1000; 1000 permutations; 100,000 Monte Carlo iterations |
| HyPhy SLAC / FEL / MEME / GARD actually ran | PASS: all four JSONs are populated (SLAC, FEL and MEME each have 217 codon rows; GARD tested 101 models). The cells swallow exceptions, but none occurred. |
| Alignment | PASS: 599 sequences x 651 columns; best model TIM2+I+R4 |

## 2. Epidemiology statistics, independently recomputed from the raw Excel and the analysis table
| Quantity | Pipeline | Independent | Status |
|---|---|---|---|
| Prevalence | 43/136 = 31.6% (Wilson 24.4-39.8%) | identical | PASS |
| Sex (22/68 F vs 21/68 M) | Fisher p = 1.0, OR 0.93 | identical | PASS |
| Species x positivity | chi2 7.35, MC p = 0.290 | chi2 7.35, MC p = 0.289 | PASS |
| Species (collapsed, n<10 pooled) | chi2 7.12, p = 0.128 | not re-run | PASS (same code path) |
| Species x band (None/Low/Mod/High), thesis replication | chi2 59.2, MC p = 0.0018 | chi2 59.2, MC p = 0.0018 | PASS, but see below |
| Intensity among positives x species | chi2 34.05, MC p = 0.0067 | chi2 34.05, MC p = 0.0064 | PASS |
| Kruskal-Wallis, intensity among positives | p = 0.072 | p = 0.072 | PASS |
| Firth: M. monachus group (SP02) vs reference | OR 0.25 (0.06-1.03), p = 0.055 | n/a | CAVEAT (borderline, not significant) |

Notes
- The earlier local run (`local_stats`) reported chi2 = 54.5 for the thesis replication. The Colab run and my independent recomputation both give 59.2, so 54.5 is superseded. The cause is probably a different band grouping and was not investigated. Use 59.2.
- 78.6% of the species x band cells have expected counts <5 (minimum expected 0.07). The thesis P<0.001 therefore rests on a sparse table. Report the Monte Carlo p and state this limitation.
- Positivity is NOT associated with species (p = 0.29) or sex. Only intensity among positives differs by species group (MC p = 0.007, Kruskal p = 0.07, 90% of cells with expected counts <5). This is the honest reading of the thesis "P<0.001" versus "no significance" conflict.
- Species labels are placeholders (SP01-SP07, with suggested names marked with †). All species-level results change in wording only, not in numbers, once final names are applied.

## 3. Sequence QC (sanger_qc.csv)
- All 10 isolates: 217 aa Rep fragment, 0 stop codons; motifs intact: HLQG@17, YCSK@55, Walker A GPPGCGKS@131, Walker B V(V/I)LDD@166, motif C ITSN@207. Pb23 has VILDD instead of VVLDD in Walker B.
- CAVEAT/OPEN: S27 read covers only 71.7% of its submitted 651-nt sequence (trimmed read 498 nt) and S7 matches its submitted sequence at only 90.4% identity (mean Phred 38.5, 122 ambiguous calls). Pb64 covers 82.9%, Pb30 94%. Only forward reads exist. The user must explain how the submitted S7/S27/Pb64 sequences were obtained, for example a missing reverse read or manual editing. Until then these isolates should be flagged in the supplement.

## 4. Phylogenetic placement
- EPA-ng on the 429-genome backbone, with clade likelihood weight ratio (LWR): Pb24, Pb30, Pb34, S7, S27 and Pb22 are GI-c (LWR 0.99-1.00); Pb25 and Pb28 are GI-c (0.90 and 0.91); Pb23 is GI-f (1.00); Pb64 is GI-d 0.56 / GI-f 0.36 / GI-e 0.02.
- Nearest neighbours (p-distance) agree with the placement: GI-c isolates are 97.7-98.9% identical to KF723390 (Italy), MK803398/MK803403 (Saudi Arabia), JX221037 (Poland) and OP039998 (Iran, rose-ringed parakeet). Pb23 and Pb64 are 98-99% identical to the Iranian budgerigar sequences KT764926-28/34 (2014) and to MK387075 (Bangladesh, GI-f).
- CAVEAT 1: only 343/599 amplicon tips and 343/429 backbone tips carry a Shah 2023 label (the other 86 backbone genomes are 'unassigned'). The automatic propagation recovered only GI-d/GI-e/GI-f as clean clades; GI-a/b/c labels come from the original Shah labels, not from tree structure. The amplicon-tree and constrained-tree columns in `kerman_clade_assignment_comparison.csv` are therefore mostly "unassigned" and all show agree_all = False. They are not evidence of disagreement. Present EPA-ng as the primary assignment and the constrained tree as a sensitivity check (it also puts Pb23 in GI-f).
- CAVEAT 2: Pb22 is placed on an edge whose subtree holds 105 backbone tips: 30 labelled GI-c and 75 unassigned (OQ995xxx and similar). Its GI-c call is supported by LWR but should be worded as "GI-c / closest to GI-c" rather than a firm sub-genotype.
- CAVEAT 3: Pb64 should be reported as "GI-d or GI-f, unresolved" (LWR 0.56 vs 0.36). Both Pb23 and Pb64 fall in the Iranian budgerigar cluster.
- Conclusion: all 10 isolates are BFDV GI (none GII); 8 or 9 of 10 fall in GI-c, 1 in GI-f, 1 ambiguous GI-d/f. Mean identity to GI-e is the highest for Pb23/Pb64 (96.9-97%) because the GI-e group is small and close, so identity alone must not be used to name the clade.

## 5. Diversity, networks, association
- Kerman: 10 haplotypes in 10 sequences, Hd = 1.00, pi = 0.0313, 60 segregating sites; all Iran (n = 31): 29 haplotypes, Hd = 0.996, pi = 0.0421; all sequences: 444 haplotypes in 599, pi = 0.0553. The "Pb30 = Pb34 identical" suspicion is NOT borne out for the submitted sequences (10 distinct haplotypes), but they are close neighbours on the tree (the longest identical branch has length 1e-6).
- BaTS-like tests (1000 permutations; the smallest possible p is 0.001): strong phylogeny-country association (AI 3.3 vs null 52.4; PS 85 vs null 404) and host-genus association (AI 15.3 vs null 54.3; PS 141 vs null 438). The host genus is for references only. Both tests are dominated by clustered submissions from single studies (for example PP329xxx from Iran 2024), so geographic and host structure may be partly sampling artefact. Wording must be careful.
- Recombination: GARD found no breakpoint (101 models, no improvement in AICc) on 41 sequences x 651 nt. A single-tree fit does not support recombination in this fragment. Because the fragment is short, this is not evidence about the whole genome.
- Selection (SLAC global dN/dS = 0.14): FEL finds 65 codons with p<0.1, all negatively selected (beta<alpha), none positively selected; MEME finds no episodic positive-selection site at p<0.1. The fragment is under strong purifying selection.

## 6. Figures reviewed (Fig3, Fig4 inspected visually)
- Fig 4 (placement bars): PASS; the placeholder names and † mark render.
- Fig 3 (focus tree): PASS on content, but layout needs fixing before submission: legend overlaps the scale bar, a host label reads "nan" (KJ413143), several reference tips show a vernacular host (cockatiel/love bird/budgerigar), and there is large blank space at the bottom. Fig 1, 2, 5 and the supplementary figures still to be reviewed after the final names are applied.

## 7. Open items needing the user
1. Final species names for the 7 host groups and the 10 isolates (affects all labels).
2. How S7, S27 and Pb64 were obtained (coverage 72-83%, identity 90% for S7).
3. Real collection dates and the ethics approval code. Also authors, PCR details, how sex was determined, and confirmation that the 21 "Green-cheeked" rows are Myiopsitta monachus.
4. Literature prevalence figures quoted in the draft (Dolatyabi, Khosravi, Haddadmarandi, Morinha) need to be re-checked against the sources.

## 8. Update 2026-11-16 (final names applied)
- Final host names from `species_names.jpg` entered in `metadata/host_mapping.csv` (all 10 isolates) and `metadata/species_groups.csv` (all 7 groups). Names are identical to the earlier suggestions; Pb34 is now *Agapornis* spp. SP03 was set to *Agapornis* spp. and SP02 to *Myiopsitta monachus* by inference. Both need the authors' confirmation.
- All figures re-rendered without placeholders or daggers. Fig 3 legend overlap and the "nan" host label were fixed; the Fig 5 legend was moved below the panel.
- RF distance between the TIM2+I+R4 and TIM3+F+I+G4 amplicon trees is 308 (from `rf_MFP_vs_TIM3.rfdist`), not the 354 recorded earlier. Log-likelihoods: -16,659 vs -16,725.
- Manuscript v3 built from these results; 17 [PENDING] markers remain, all of which need information from the authors.

## 9. Update 2026-11-16 (editorial figure comments applied)
Source: `comments.txt` plus the author's follow-up (species names instead of IDs, no overlapping text, no text on Fig 4 bars, roomier Fig 5B).
| Figure | Change | Check |
|---|---|---|
| Fig 2 | Redrawn as host group -> PCR result -> sex for all 136 birds; PCR band intensity removed; the three columns have identical height; file renamed `Fig2_alluvial_host_result_sex` | Block totals 73/21/20/10/5/4/3 = 136; positives 25/2/6/5/2/2/1 = 43; 22 F + 21 M; matches `prevalence_by_species.csv` and `association_tests.json` |
| Fig 3 | Kept (approved); tip labels now host species + country + year; vernacular GenBank hosts mapped to binomials | Same topology and clade strip as before |
| Fig 4 | Bars sorted by clade, no text on bars, hatched remainder = LWR on unlabelled branches; caption gives purpose, rationale and take-home message | LWR values equal `placement_lwr_by_clade.csv` (GI-c 0.90-1.00; Pb23 GI-f 1.00; Pb64 GI-d 0.56 / GI-e 0.02 / GI-f 0.36) |
| Fig 5 | Heat-map unchanged in content; haplotype network redrawn full-width with distance-based layout, no node/link/label overlaps; species-name labels | 16 haplotypes and the same edge weights as `haplotype_msn`; haplotype membership unchanged |
| Fig S2 | Bars ranked by Q20, host species on the y axis | Values equal `sanger_qc.csv` |
| Fig S3 | Panel B added: residues at the five Rep motifs in the ten isolates | Equals `rep_motifs_kerman.csv` (only Pb23 differs: VILDD) |
| All | One colour token system (`FIGURE_COLOUR_KEY.md`); no hex shared by two meanings | - |
Rendering was checked visually for every figure; `verify_reproducibility.py` still passes 45/45.
Caveat: because isolate IDs are no longer printed, rows with the same species (e.g. five *N. hollandicus*) are told apart only by their order, which is stated in the captions of Fig 4 and Fig S3 (and Table S6 links IDs to hosts).
