# BFDV in pet psittacines, Kerman (Iran): analysis code and results

Code, inputs and results for the manuscript *Subclinical beak and feather disease virus infection in pet
psittacine birds in Kerman, south-eastern Iran: prevalence, risk factors and full-genome-guided genotyping of
partial Rep sequences*. Ten partial ORF1/Rep sequences are in GenBank (PX634653-PX634662).

## What is here
| Path | Content |
|---|---|
| `data/` | The 10 Kerman sequences (FASTA), raw Sanger chromatograms (`ab1/`), and the 598-sequence alignment used in the thesis |
| `metadata/` | Survey table (136 birds), host and species name tables (final names in `FINAL_*` columns), Shah et al. 2023 sub-clade labels, reference metadata |
| `scripts/` | Python modules (`pbfd_stats.py`, `pbfd_phylo.py`, `pbfd_figures.py`, `pbfd_names.py`) and entry points (`prepare_inputs.py`, `make_figures.py`, `make_supplement.py`) |
| `PBFD_analysis.ipynb`, `build_notebook.py` | Google Colab notebook that runs the whole analysis (21 code cells), and the script that generates it |
| `headless/`, `make_headless.py` | The same 21 cells as plain scripts for unattended runs (used with the Google Colab CLI) |
| `results/` | Everything the Colab run produced: trees, alignments, placement, HyPhy JSON, tables, figures, statistics, `versions.txt` |
| `manuscript/` | `build_manuscript.js` and `build_cover_letter.js` (docx-js) and `references.json` |
| `verify/` | Independent verification script and its latest report |

## Reproducing
1. **Statistics, figures, supplementary tables (minutes, any computer).**
   ```
   pip install -r requirements.txt
   python scripts/make_figures.py          # re-renders every figure from results/ and metadata/
   python scripts/make_supplement.py       # writes Supplementary_Tables.xlsx
   python verify/verify_reproducibility.py # 45 independent checks, writes verify/verification_report.txt
   ```
2. **Full phylogenetic analysis (about 4 h on a 2-CPU Colab VM).** Open `PBFD_analysis.ipynb` in Google Colab,
   upload `PBFD_colab_inputs.zip` when asked and run all cells. It installs the tools below, queries NCBI for the
   reference records (the notebook asks for your own e-mail address, which NCBI requires; the unattended run reported in the paper queried NCBI without one),
   and writes `results.zip`. The random seed is fixed (20261103). Re-running gives the same tree topologies only
   if the same tool versions are used; bootstrap values and ModelFinder output can differ slightly between versions.
3. **Unattended run with the Colab CLI:** `headless/_launcher.py` and `_start_blocking.py` run cells 1-21 in order
   (see comments in the files).

Tool versions of the reported run (`results/versions.txt`): IQ-TREE 3.1.4, MAFFT 7.526, RAxML-NG 2.0.4,
EPA-ng 0.3.8, HyPhy 2.5.103, Python 3.13.18, Biopython 1.88, NumPy 2.4.6, pandas 2.2.3, SciPy 1.16.3, statsmodels 0.15.0.

## Notes on the data
- Original thesis files (Excel sheets, PDF chapters, MEGA trees) are not redistributed. `scripts/prepare_inputs.py`
  documents how the inputs in `data/` and `metadata/` were derived from them.
- Host names follow the list supplied by the authors. The 20-bird survey group is recorded as *Agapornis* spp. and
  the 21-bird group as *Myiopsitta monachus*; both need confirmation by the authors.
- The S7, S27 and Pb64 submitted sequences are only partly supported by the forward Sanger reads (`results/tables/sanger_qc.csv`).
- Statistical conventions: Wilson 95% CIs; Monte Carlo exact tests with 100,000 permutations, margins fixed; Firth
  penalised logistic regression; HyPhy significance p <= 0.1; 1,000 permutations for phylogeny-trait tests.

## Verification summary
`verify/verification_report.txt` lists each check (survey statistics recomputed from the raw Excel values, Firth
regression re-implemented independently, alignment rows compared with the submitted FASTA, Rep translation and motifs,
diversity indices, placement re-derived from the `.jplace` file, HyPhy outputs, all 589 reference host/country values
compared with the raw GenBank records, and the numbers printed in the manuscript). The heavy tools (IQ-TREE, EPA-ng, HyPhy)
are not re-run by this script; their saved outputs are checked for consistency instead.

## License
Code: MIT (see `LICENSE`; the copyright holder has to be confirmed by the authors). Sequence data: see the GenBank records.
