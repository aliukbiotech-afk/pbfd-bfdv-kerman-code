# Figure colour key and labelling rules

One colour = one meaning, identical in every main and supplementary figure. The tokens live in
`scripts/pbfd_figures.py` (`CLADE_COL`, `KERMAN_COL`/`IRAN_COL`/`OTHER_COL`, `HOST_COL`, `RESULT_COL`, `SEX_COL`).

| Meaning | Token | Hex | Used in |
|---|---|---|---|
| Sub-clade GI-a | `CLADE_COL["GI-a"]` | `#0072B2` | Fig 3, 4; Fig S1, S4 |
| GI-b | | `#56B4E9` | same |
| GI-c | | `#009E73` | same |
| GI-d | | `#D9B300` | same |
| GI-e | | `#332288` | same |
| GI-f | | `#882255` | same |
| GII-a / GII-b / unassigned | | `#8C8C8C` / `#4D4D4D` / `#DDDDDD` | Fig S1, S4 |
| Sequence from this study (Kerman) | `KERMAN_COL` | `#D55E00` | Fig 3, 5; Fig S1, S4 |
| Other Iranian strains | `IRAN_COL` | `#CC79A7` | Fig 3, 5; Fig S1, S4 |
| Nearest non-Iranian strains | `OTHER_COL` | `#707070` | Fig 5 |
| Host group (7 groups) | `HOST_COL` | `#264653 #2A9D8F #8AB17D #E9C46A #F4A261 #E76F51 #9A9A9A` | Fig 1, 2; Fig S2 |
| PCR result | `RESULT_COL` | positive `#2B2B2B`, negative `#D0D0D0` | Fig 2 |
| Sex | `SEX_COL` | female `#5B6770`, male `#A9B4BC` | Fig 1, 2 |

Clade, origin and host colours were chosen so that no hex value is shared between two meanings (earlier drafts used
vermilion for both "this study" and GI-f, and pink for both "other Iranian" and GI-d).

## Labels
* Specimens are labelled by **host species** (binomial, or common name + binomial for host groups).
  Vernacular host names found in GenBank records (budgerigar, cockatiel, love bird, ...) are converted to binomials.
* Isolate IDs and accessions are not printed in the figures; they are in Supplementary Table S6 (isolates) and S7 (references).
* Fig 5 abbreviates the genus (e.g. *N. hollandicus*) because of space.
