# BFDV in pet psittacines, Kerman (Iran): analysis code and results

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aliukbiotech-afk/pbfd-bfdv-kerman-code/blob/main/PBFD_analysis.ipynb)

Analysis code, inputs and results for the study *Subclinical beak and feather disease virus infection in pet
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
| `verify/` | Independent verification script (writes `verify/verification_report.txt`) |
| `FIGURE_COLOUR_KEY.md` | Colour tokens and labelling rules shared by all figures |

## Reproducing
1. **Statistics, figures, supplementary tables (minutes, any computer).**
   ```
   pip install -r requirements.txt
   python scripts/make_figures.py          # re-renders every figure from results/ and metadata/
   python scripts/make_supplement.py       # writes Supplementary Tables S1-S7 (Supplementary_Tables.xlsx)
   python verify/verify_reproducibility.py # 45 independent checks, writes verify/verification_report.txt
   ```
2. **Full phylogenetic analysis (about 4 h on a 2-CPU Colab VM).** Click the *Open in Colab* badge above and run all
   cells. The notebook clones this repository for its inputs, so nothing has to be uploaded. It installs the tools below, queries NCBI for the
   reference records (the notebook asks for your own e-mail address, which NCBI requires; the unattended run reported in the paper queried NCBI without one),
   and writes `results.zip`. The random seed is fixed (20261103). Re-running gives the same tree topologies only
   if the same tool versions are used; bootstrap values and ModelFinder output can differ slightly between versions.
3. **Unattended run with the Colab CLI:** `colab exec -f headless/_start_blocking.py` clones the repository and runs
   cells 1-21 in order through `headless/_launcher.py` (see comments in the files).

Tool versions of the reported run (`results/versions.txt`): IQ-TREE 3.1.4, MAFFT 7.526, RAxML-NG 2.0.4,
EPA-ng 0.3.8, HyPhy 2.5.103, Python 3.13.18, Biopython 1.88, NumPy 2.4.6, pandas 2.2.3, SciPy 1.16.3, statsmodels 0.15.0.

## Figures
`scripts/make_figures.py` writes Fig 1–5 and Fig S1–S4 (TIFF 600 dpi, PDF, PNG) into `results/figures/`. Figure labels use host species
instead of isolate IDs; the ID-to-host table is Supplementary Table S6 and the 589 references are in Table S7. The haplotype-network
layout (Fig 5B) is computed deterministically (fixed seed) and takes about a minute.

## How to cite
Cite the manuscript and the archived release of this repository: (DOI to be added after the Zenodo release).

## License
Code: MIT (see `LICENSE`; the copyright holder has to be confirmed by the authors). Sequence data: see the GenBank records.
