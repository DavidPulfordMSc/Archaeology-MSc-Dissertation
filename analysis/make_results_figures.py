#!/usr/bin/env python3
"""
make_results_figures.py

Generate the descriptive figures and tables reported in the dissertation Results section.

This script is deliberately limited to descriptive summaries. It does not attempt to
pool effect sizes or infer causal directionality.

Inputs (preferred):
- study_level_n33.csv (study-level, locus-extractable subset)
- gene_level_summary.csv (gene-level extraction table)

Outputs:
- PNG and SVG figures
- CSV tables supporting figures

Author: David Pulford
Year: 2026
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import List, Dict, Optional

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from wordcloud import WordCloud


# -----------------------------
# Configuration
# -----------------------------

FIG_DPI = 300
FIG_FORMATS = ("png", "svg")


# -----------------------------
# Utilities
# -----------------------------

def ensure_dir(path: str) -> str:
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def save_fig(fig: plt.Figure, out_dir: str, name: str) -> None:
    for ext in FIG_FORMATS:
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), dpi=FIG_DPI, bbox_inches="tight")


def read_table(path: str) -> pd.DataFrame:
    """
    Read CSV or Excel into a dataframe.
    """
    lower = path.lower()
    if lower.endswith(".csv"):
        return pd.read_csv(path)
    if lower.endswith(".xlsx") or lower.endswith(".xls"):
        return pd.read_excel(path)
    raise ValueError(f"Unsupported file type: {path}")


def clean_text_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip()


def drop_accidental_header_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove rows where a cell contains the column name itself.
    This occasionally happens when copying from spreadsheets.
    """
    for c in df.columns:
        if df[c].dtype == object:
            df = df[df[c].astype(str).str.strip().str.lower() != str(c).strip().lower()]
    return df


# -----------------------------
# Normalisation helpers
# -----------------------------

DOMAIN_MAP: Dict[str, str] = {
    "cognitive": "Cognition",
    "cognition": "Cognition",
    "psychiatric": "Psychiatric",
    "neurodevelopment": "Neurodevelopmental",
    "neurodevelopmental": "Neurodevelopmental",
    "circadian": "Circadian",
    "chronotype": "Circadian",
    "sleep": "Circadian",
    "addiction": "Addiction",
    "smoking": "Addiction",
    "pain": "Pain",
    "social": "Social",
    "language": "Cognition",
    "brain": "Brain structure/imaging",
    "neuroimaging": "Brain structure/imaging",
    "cranial": "Brain structure/imaging",
    "none": "None/Context only",
    "mixed": "Mixed/multiple",
    "multiple": "Mixed/multiple",
}

def normalise_domains(raw: str) -> List[str]:
    if pd.isna(raw) or str(raw).strip() == "":
        return []
    parts = re.split(r"[;,/|]+", str(raw))
    out: List[str] = []
    for p in parts:
        t = p.strip().lower()
        if not t:
            continue
        mapped: Optional[str] = None
        for k, v in DOMAIN_MAP.items():
            if k in t:
                mapped = v
                break
        out.append(mapped if mapped else "Other")
    return sorted(set(out))


def normalise_study_type(raw: str) -> str:
    if pd.isna(raw):
        return "Not reported"
    t = str(raw).lower()

    if any(x in t for x in ["gwas", "phewas", "uk biobank"]):
        return "GWAS/PheWAS association"
    if any(x in t for x in ["heritab", "ldsc", "partition"]):
        return "Heritability/enrichment"
    if any(x in t for x in ["eqtl", "gtex", "expression", "transcript"]):
        return "Gene expression / eQTL"
    if "splic" in t:
        return "Splicing / isoform"
    if any(x in t for x in ["organoid", "crispr", "functional", "assay", "reporter"]):
        return "Functional validation"
    if any(x in t for x in ["selection", "adaptive", "introgression", "population"]):
        return "Population genetics / introgression"
    if any(x in t for x in ["review", "framework"]):
        return "Review / synthesis"
    if any(x in t for x in ["archaeolog", "petri", "tar"]):
        return "Archaeological context"

    return "Other / mixed"


def normalise_loci_level(raw: str) -> str:
    if pd.isna(raw):
        return "Not reported"
    t = str(raw).lower()

    if any(x in t for x in ["genome-wide", "burden", "enrichment", "score"]):
        if any(x in t for x in ["rs", "snp", "variant", "haplotype", "region", "gene"]):
            return "Mixed (genome-wide + loci)"
        return "Genome-wide only"
    if any(x in t for x in ["rs", "snp", "variant"]):
        return "Variant-level (SNP/rsID)"
    if any(x in t for x in ["haplotype", "region", "tract", "segment"]):
        return "Region/haplotype"
    if any(x in t for x in ["gene", "locus"]):
        return "Gene-level"
    return "Mixed/unclear"


# -----------------------------
# Plot helpers
# -----------------------------

def bar_plot(categories: List[str], counts: List[int], title: str, ylabel: str, out_path_base: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(categories, counts)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    save_fig(fig, os.path.dirname(out_path_base), os.path.basename(out_path_base))
    plt.close(fig)


def wordcloud_from_frequencies(freq: Dict[str, int], title: str, out_path_base: str) -> None:
    wc = WordCloud(width=1600, height=800, background_color="white").generate_from_frequencies(freq)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.imshow(wc, interpolation="bilinear")
    ax.axis("off")
    ax.set_title(title)
    save_fig(fig, os.path.dirname(out_path_base), os.path.basename(out_path_base))
    plt.close(fig)


# -----------------------------
# Main figure generation
# -----------------------------

def make_study_figures(study_df: pd.DataFrame, out_dir: str) -> None:
    # Study type distribution
    st = study_df.copy()
    st["Study_type_norm"] = st["Study_type"].apply(normalise_study_type)
    counts = st["Study_type_norm"].value_counts()
    bar_plot(
        counts.index.tolist(),
        counts.values.tolist(),
        "Studies by study type (normalised)",
        "Number of studies",
        os.path.join(out_dir, "fig03_study_types")
    )
    counts.reset_index().rename(columns={"index": "Study_type", "Study_type_norm": "Count"}).to_csv(
        os.path.join(out_dir, "fig03_study_types_counts.csv"), index=False
    )

    # Trait domain distribution (each domain counted once per study)
    td = study_df[["Study_ID", "Phenotype_domains"]].copy()
    td["domain_list"] = td["Phenotype_domains"].apply(normalise_domains)
    td_long = td.explode("domain_list").dropna(subset=["domain_list"])
    td_long = td_long.drop_duplicates(["Study_ID", "domain_list"])
    dcounts = td_long["domain_list"].value_counts()

    bar_plot(
        dcounts.index.tolist(),
        dcounts.values.tolist(),
        "Trait domain distribution (normalised)",
        "Number of studies (domain counted once per study)",
        os.path.join(out_dir, "fig04_trait_domains")
    )
    dcounts.reset_index().rename(columns={"index": "Phenotype_domain", "domain_list": "Count"}).to_csv(
        os.path.join(out_dir, "fig04_trait_domains_counts.csv"), index=False
    )

    # Loci resolution
    lr = study_df.copy()
    lr["Loci_level_norm"] = lr["Loci_level"].apply(normalise_loci_level)
    lcounts = lr["Loci_level_norm"].value_counts()

    bar_plot(
        lcounts.index.tolist(),
        lcounts.values.tolist(),
        "Resolution of introgression evidence (normalised)",
        "Number of studies",
        os.path.join(out_dir, "fig05_loci_resolution")
    )
    lcounts.reset_index().rename(columns={"index": "Loci_level_norm", "Loci_level_norm": "Count"}).to_csv(
        os.path.join(out_dir, "fig05_loci_resolution_counts.csv"), index=False
    )


def make_gene_figures(gene_df: pd.DataFrame, out_dir: str) -> None:
    g = gene_df.copy()

    # Standardise basic fields
    for c in ["Study_ID", "Gene", "Phenotype_domain", "Evidence_type"]:
        g[c] = clean_text_series(g[c])

    g = drop_accidental_header_rows(g)
    g = g[g["Study_ID"].astype(str).str.strip() != ""]
    g = g[g["Gene"].astype(str).str.strip() != ""]
    g["Gene"] = g["Gene"].str.upper()

    # Table 1: gene recurrence (count distinct studies)
    gene_counts = (
        g.drop_duplicates(["Study_ID", "Gene"])
         .groupby("Gene")["Study_ID"]
         .nunique()
         .sort_values(ascending=False)
    )
    gene_counts.rename("n_studies").reset_index().to_csv(
        os.path.join(out_dir, "table01_gene_recurrence.csv"), index=False
    )

    # Figure 6: top genes bar chart
    top_n = 20
    top = gene_counts.head(top_n).reset_index()
    top.columns = ["Gene", "n_studies"]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(top["Gene"], top["n_studies"])
    ax.set_title(f"Top {top_n} genes (counted by distinct studies)")
    ax.set_ylabel("Number of studies")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    save_fig(fig, out_dir, "fig06_top_genes")
    plt.close(fig)

    # Figure 7: gene word cloud (study-weighted)
    wordcloud_from_frequencies(
        gene_counts.to_dict(),
        "Gene word cloud (weighted by number of distinct studies)",
        os.path.join(out_dir, "fig07_wordcloud_genes")
    )

    # Figure 13: evidence types (study-weighted)
    ecounts = (
        g.drop_duplicates(["Study_ID", "Evidence_type"])
         .groupby("Evidence_type")["Study_ID"]
         .nunique()
         .sort_values(ascending=False)
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(ecounts.index.tolist(), ecounts.values.tolist())
    ax.set_title("Evidence types among gene-naming studies")
    ax.set_ylabel("Number of studies")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    save_fig(fig, out_dir, "fig13_evidence_types")
    plt.close(fig)
    ecounts.rename("n_studies").reset_index().to_csv(
        os.path.join(out_dir, "fig13_evidence_types_counts.csv"), index=False
    )

    # Figure 8: phenotype-domain word cloud within gene-naming subset (study-weighted)
    domain_counts = (
        g.drop_duplicates(["Study_ID", "Phenotype_domain"])
         .groupby("Phenotype_domain")["Study_ID"]
         .nunique()
         .sort_values(ascending=False)
    )
    wordcloud_from_frequencies(
        domain_counts.to_dict(),
        "Phenotype-domain word cloud (weighted by number of studies)",
        os.path.join(out_dir, "fig08_wordcloud_domains")
    )

    # Figures 9 and 10: domain-specific gene word clouds
    def domain_gene_wordcloud(domain_label: str, out_name: str) -> None:
        sub = g[g["Phenotype_domain"].str.strip().str.lower() == domain_label.lower()].copy()
        if sub.empty:
            return
        freq = (
            sub.drop_duplicates(["Study_ID", "Gene"])
               .groupby("Gene")["Study_ID"]
               .nunique()
               .sort_values(ascending=False)
        )
        if len(freq) < 3:
            return
        wordcloud_from_frequencies(
            freq.to_dict(),
            f"Gene word cloud within domain: {domain_label} (study-weighted)",
            os.path.join(out_dir, out_name)
        )

    domain_gene_wordcloud("Neurodevelopmental", "fig09_wordcloud_genes_neurodevelopmental")
    domain_gene_wordcloud("Mixed/multiple", "fig10_wordcloud_genes_mixed_multiple")

    # Figure 11: gene x domain heatmap (top genes)
    top_heatmap = gene_counts.head(25).index.tolist()
    sub = g[g["Gene"].isin(top_heatmap)].copy()
    co = (
        sub.drop_duplicates(["Study_ID", "Gene", "Phenotype_domain"])
           .groupby(["Gene", "Phenotype_domain"])["Study_ID"]
           .nunique()
           .reset_index(name="n_studies")
    )
    pivot = co.pivot(index="Gene", columns="Phenotype_domain", values="n_studies").fillna(0)
    pivot.to_csv(os.path.join(out_dir, "fig11_gene_domain_matrix.csv"))

    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Gene x phenotype domain (counts by studies)")
    fig.colorbar(im, ax=ax, label="Number of studies")
    plt.tight_layout()
    save_fig(fig, out_dir, "fig11_gene_domain_heatmap")
    plt.close(fig)

    # Figure 12: gene-domain network (genes with >= 2 studies)
    min_studies = 2
    keep = gene_counts[gene_counts >= min_studies].index.tolist()
    subn = g[g["Gene"].isin(keep)].copy()

    edges = (
        subn.drop_duplicates(["Study_ID", "Gene", "Phenotype_domain"])
            .groupby(["Gene", "Phenotype_domain"])["Study_ID"]
            .nunique()
            .reset_index(name="weight")
    )

    G = nx.Graph()
    for gene in keep:
        G.add_node(gene, node_type="gene")
    domains = sorted(subn["Phenotype_domain"].dropna().unique())
    for dom in domains:
        G.add_node(dom, node_type="domain")

    for _, row in edges.iterrows():
        if pd.notna(row["Phenotype_domain"]):
            G.add_edge(row["Gene"], row["Phenotype_domain"], weight=int(row["weight"]))

    pos = nx.spring_layout(G, seed=42, k=0.8)

    fig, ax = plt.subplots(figsize=(14, 10))
    gene_nodes = [n for n in G.nodes if G.nodes[n].get("node_type") == "gene"]
    dom_nodes = [n for n in G.nodes if G.nodes[n].get("node_type") == "domain"]

    nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, node_size=600, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=dom_nodes, node_size=800, ax=ax)
    nx.draw_networkx_edges(G, pos, width=1, alpha=0.6, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

    ax.set_title(f"Gene-phenotype domain network (genes with >={min_studies} studies)")
    ax.axis("off")
    plt.tight_layout()
    save_fig(fig, out_dir, "fig12_gene_domain_network")
    plt.close(fig)


# -----------------------------
# Entrypoint
# -----------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce dissertation Results figures (descriptive).")
    parser.add_argument("--study-data", required=True, help="Path to study-level dataset (CSV or Excel).")
    parser.add_argument("--gene-data", required=True, help="Path to gene-level dataset (CSV or Excel).")
    parser.add_argument("--output", default="outputs/figures", help="Directory for output figures and tables.")
    args = parser.parse_args()

    out_dir = ensure_dir(args.output)

    study_df = read_table(args.study_data)
    gene_df = read_table(args.gene_data)

    # Basic guard: keep only populated studies
    if "Study_ID" in study_df.columns:
        study_df["Study_ID"] = clean_text_series(study_df["Study_ID"])
        study_df = study_df[study_df["Study_ID"] != ""].copy()

    make_study_figures(study_df, out_dir)
    make_gene_figures(gene_df, out_dir)

    print(f"Done. Outputs written to: {out_dir}")


if __name__ == "__main__":
    main()
