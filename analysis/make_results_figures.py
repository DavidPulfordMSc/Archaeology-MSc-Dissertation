import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from wordcloud import WordCloud

# --- DISSERTATION CONFIG ---
# Paths for the datasets used in the Results section
STUDY_DATA = "study_level_n33.csv"  # The 33-study locus-resolved subset
GENE_DATA = "gene_level_summary.csv" 
OUT_DIR = "outputs/figures"

# Ensure output folders exist
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)
if not os.path.exists(OUT_DIR + "/tables"):
    os.makedirs(OUT_DIR + "/tables")

# Domain order as defined in the dissertation text (Fig 11)
# Includes specific categories: Circadian, Pain, Psychiatric, etc.
DOMAIN_ORDER = [
    "Circadian", "Mixed / multiple", "Neurodevelopmental", 
    "Pain", "Psychiatric", "Social"
]

# Map for cleaning the Phenotype_domains column
DOMAIN_MAP = {
    "cognitive": "Cognition", "psychiatric": "Psychiatric",
    "neurodevelopment": "Neurodevelopmental", "circadian": "Circadian",
    "sleep": "Circadian", "pain": "Pain", "social": "Social",
    "mixed": "Mixed/multiple", "multiple": "Mixed/multiple"
}

# --- FUNCTIONS ---

def get_study_type(text):
    """Categorises study design for Fig 03 (Composition of Evidence)"""
    t = str(text).lower()
    if any(x in t for x in ["gwas", "phewas"]): return "GWAS/PheWAS association"
    if any(x in t for x in ["eqtl", "expression"]): return "Gene expression / eQTL"
    if any(x in t for x in ["functional", "assay", "organoid"]): return "Functional validation"
    if any(x in t for x in ["selection", "introgression"]): return "Population genetics"
    return "Other / mixed"

def get_resolution(text):
    """Categorises biological resolution for Fig 05"""
    t = str(text).lower()
    if "variant" in t or "snp" in t: return "Variant-level (SNP/rsID)"
    if "gene" in t: return "Gene-level"
    if "region" in t or "haplotype" in t: return "Region/haplotype"
    return "Genome-wide only"

# --- MAIN ANALYSIS ---

print("Reading dissertation datasets...")
# Loading the core 33 studies identified in the PRISMA flow (Fig 2)
df_study = pd.read_csv(STUDY_DATA)
df_gene = pd.read_csv(GENE_DATA)

# 1. Evidence Base Composition (Fig 03)
print("Plotting Figure 3: Study Types...")
df_study['Type_Norm'] = df_study['Study_type'].apply(get_study_type)
type_counts = df_study['Type_Norm'].value_counts()
plt.figure(figsize=(10, 6))
type_counts.plot(kind='bar', color='skyblue')
plt.title("Figure 3: Composition of the Evidence Base (n=33)")
plt.ylabel("Number of Studies")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig03_study_types.png", dpi=300)

# 2. Trait Domain Distribution (Fig 04)
# Note: One study can cover multiple domains
print("Plotting Figure 4: Trait Domains...")
# Split domains and explode to get counts per domain
df_study['domain_list'] = df_study['Phenotype_domains'].str.split(';')
domains_exploded = df_study.explode('domain_list')
domain_counts = domains_exploded['domain_list'].str.strip().value_counts()
plt.figure(figsize=(10, 6))
domain_counts.plot(kind='bar', color='salmon')
plt.title("Figure 4: Trait-Domain Distribution")
plt.ylabel("Number of Studies")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig04_trait_domains.png", dpi=300)

# 3. Biological Resolution (Fig 05)
print("Plotting Figure 5: Biological Resolution...")
df_study['Res_Norm'] = df_study['Loci_level'].apply(get_resolution)
res_counts = df_study['Res_Norm'].value_counts()
plt.figure(figsize=(10, 6))
res_counts.plot(kind='bar', color='lightgreen')
plt.title("Figure 5: Resolution of Introgression Evidence")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig05_loci_resolution.png", dpi=300)

# 4. Gene Recurrence Analysis (Table 1 & Fig 06)
print("Analysing gene recurrence...")
# Clean gene names (ensure uppercase for OAS1, TLR, etc.)
df_gene['Gene'] = df_gene['Gene'].astype(str).str.strip().str.upper()
# Count unique studies per gene
gene_stats = df_gene.drop_duplicates(['Study_ID', 'Gene']).groupby('Gene')['Study_ID'].nunique().sort_values(ascending=False)
gene_stats.to_csv(f"{OUT_DIR}/tables/table01_gene_recurrence.csv")

# 5. Top Genes Word Cloud (Fig 07)
print("Generating Fig 7: Gene Word Cloud...")
word_freq = gene_stats.to_dict()
wc = WordCloud(width=1200, height=600, background_color="white").generate_from_frequencies(word_freq)
plt.figure(figsize=(10, 5))
plt.imshow(wc, interpolation='bilinear')
plt.axis("off")
plt.title("Figure 7: Gene Emphasis (Weighted by Study Count)")
plt.savefig(f"{OUT_DIR}/fig07_gene_frequency.png")

# 6. Gene-Domain Heatmap (Fig 11)
print("Generating Fig 11: Heatmap...")
top_genes = gene_stats.head(20).index.tolist()
subset = df_gene[df_gene['Gene'].isin(top_genes)]
# Pivot to get Gene vs Domain counts
heatmap_data = subset.drop_duplicates(['Study_ID', 'Gene', 'Phenotype_domain'])
pivot = heatmap_data.groupby(['Gene', 'Phenotype_domain'])['Study_ID'].nunique().unstack(fill_value=0)
# Reorder columns to match DOMAIN_ORDER
pivot = pivot.reindex(columns=DOMAIN_ORDER, fill_value=0)

plt.figure(figsize=(12, 8))
plt.imshow(pivot.values, cmap="YlOrRd", aspect="auto")
plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=45, ha="right")
plt.yticks(range(len(pivot.index)), pivot.index)
plt.colorbar(label="Number of Studies")
plt.title("Figure 11: Gene x Phenotype Domain Heatmap")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig11_gene_domain_heatmap.png")

# 7. Gene-Domain Network (Fig 12)
print("Generating Fig 12: Network Diagram...")
G = nx.Graph()
# Only use genes appearing in 2+ studies for clarity
reliable_genes = gene_stats[gene_stats >= 2].index.tolist()
net_data = df_gene[df_gene['Gene'].isin(reliable_genes)]

for _, row in net_data.iterrows():
    if pd.notna(row['Phenotype_domain']):
        G.add_edge(row['Gene'], row['Phenotype_domain'])

plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_size=800, node_color='lightgrey', font_size=8)
plt.title("Figure 12: Gene-Phenotype Domain Relationships")
plt.savefig(f"{OUT_DIR}/fig12_gene_domain_network.png")

print(f"All figures generated in {OUT_DIR}")
