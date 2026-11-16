"""Single source of truth for host / species labels.

Until FINAL_* columns are filled in the metadata CSVs, every label is rendered as
    <placeholder> (suggested: <name>?)†
Once FINAL_* is filled, the final name is used everywhere (figures, tables, docx).
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

DAGGER = "†"


def _blank(x) -> bool:
    return x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() == ""


def load_species_groups(meta_dir: Path) -> pd.DataFrame:
    return pd.read_csv(Path(meta_dir) / "species_groups.csv", dtype=str, keep_default_na=False)


def load_host_mapping(meta_dir: Path) -> pd.DataFrame:
    return pd.read_csv(Path(meta_dir) / "host_mapping.csv", dtype=str, keep_default_na=False)


def species_label(group_id: str, groups: pd.DataFrame, style: str = "plot") -> str:
    row = groups.loc[groups.group_id == group_id]
    if row.empty:
        return str(group_id)
    r = row.iloc[0]
    if not _blank(r.FINAL_SPECIES):
        return r.FINAL_SPECIES if _blank(r.FINAL_COMMON) or style == "sci" else f"{r.FINAL_COMMON} ({r.FINAL_SPECIES})"
    if style == "long":
        return f"{group_id} (suggested: {r.suggested_species}; confidence {r.confidence}){DAGGER}"
    if r.confidence.startswith("high"):
        return f"{r.suggested_short}{DAGGER}"          # confident suggestion, flagged until confirmed
    return f"{group_id}: {r.suggested_short}?{DAGGER}"


def host_label(voucher: str, hosts: pd.DataFrame, style: str = "plot") -> str:
    row = hosts.loc[hosts.voucher.str.lower() == str(voucher).lower()]
    if row.empty:
        return str(voucher)
    r = row.iloc[0]
    if not _blank(r.FINAL_HOST):
        return r.FINAL_HOST
    if style == "long":
        sug = r.suggested_host if not _blank(r.suggested_host) else "unknown"
        return f"{r.placeholder} (suggested: {sug}; {r.confidence}){DAGGER}"
    return f"{r.placeholder}: {r.suggested_short}?{DAGGER}"


def kerman_tip_label(tip: str, hosts: pd.DataFrame) -> str:
    """'PX634656_Kerman_S7' -> '<host> · Iran-Kerman'"""
    parts = tip.split("_Kerman_")
    if len(parts) != 2:
        return tip
    acc, v = parts
    return f"{host_label(v, hosts)} · Iran-Kerman"


def pending_names(meta_dir: Path) -> dict:
    g, h = load_species_groups(meta_dir), load_host_mapping(meta_dir)
    return dict(species_groups_pending=int(g.FINAL_SPECIES.map(_blank).sum()),
                isolate_hosts_pending=int(h.FINAL_HOST.map(_blank).sum()))


def species_legend_table(meta_dir: Path) -> pd.DataFrame:
    """Footnote table for figures: placeholder -> full suggestion and evidence."""
    g = load_species_groups(meta_dir)
    return g[["group_id", "host_excel", "suggested_species", "confidence", "FINAL_SPECIES"]]
