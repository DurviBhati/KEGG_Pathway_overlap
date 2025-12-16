#!/usr/bin/env python
# coding: utf-8

# In[43]:


#!pip install venn


# In[32]:


##============================================
## Importing important libraries 
##============================================

import pandas as pd
import requests
import io
import itertools
import matplotlib.pyplot as plt
from venn import venn
from itertools import combinations


# In[68]:


##==========================================
## Fetching data from KEGG 
##==========================================

def load_kegg_data(url):
    """
    Fetching KEGG REST data and load it as a Dataframe
    """
    print(f"Fetching data from {url}...")
    response = requests.get(url)
    if response.status_code == 200:
        return pd.read_csv(io.StringIO(response.text), sep='\t', header=None)
    else:
        raise Exception(f"Failed to fetch data from {url}")

# Load Pathways
pathway = load_kegg_data("https://rest.kegg.jp/list/pathway/hsa")
pathway.columns = ['PATHWAY_ID', 'PATHWAY_NAME']

# Load Gene-Pathway Link
gene_pathway = load_kegg_data("https://rest.kegg.jp/link/pathway/hsa")
gene_pathway.columns = ['GENE_ID', 'PATHWAY_ID']
# Removing 'path:' prefix so it matches pathway list
gene_pathway["PATHWAY_ID"] = gene_pathway["PATHWAY_ID"].str.replace("path:", "", regex=False)


# Load Gene Info 
gene_info = load_kegg_data("https://rest.kegg.jp/list/hsa")
gene_info.columns = ['GENE_ID', 'GENE_INFO','TYPE','TYPE_DESCRIPTION']

print("-"*100)
print(pathway.head())
print("-"*100)
print(gene_pathway.head())
print("-"*100)
print(gene_info.head())
print("-"*100)


# In[54]:


# ==========================================
#  Data Cleaning and Parsing
# ==========================================
print("Parsing gene information...")

# The gene list from KEGG comes as "Symbol; Description" in one column. 
# We need to split it.

def extract_gene_symbol(row):
    """
    Extracts gene symbol and gene name from KEGG TYPE_DESCRIPTION.
    
    TYPE_DESCRIPTION format:
    SYMBOL1, SYMBOL2, ... ; Full gene name
    """
    desc = row["TYPE_DESCRIPTION"]

    if pd.isna(desc):
        return pd.Series(["NA", "NA"])

    # Split symbol(s) and description
    if ";" in desc:
        symbol_part, gene_name = desc.split(";", 1)
    else:
        symbol_part = desc
        gene_name = "NA"

    # First symbol is the primary gene symbol
    gene_symbol = symbol_part.split(",")[0].strip()

    return pd.Series([gene_symbol, gene_name.strip()])

gene_info[["GENE_SYMBOL", "GENE_NAME"]] = gene_info.apply(
    extract_gene_symbol, axis=1
)

genes = gene_info[
    ["GENE_ID", "GENE_SYMBOL", "GENE_NAME", "TYPE"]
]

print(genes.head())
print("Unique gene symbols:", genes["GENE_SYMBOL"].nunique())

# Drop rows where gene symbol could not be extracted
genes = genes.dropna(subset=["GENE_SYMBOL"])

# Very rare edge cases
genes = genes[~genes["GENE_SYMBOL"].isin(["-", "NA"])]

print("After cleaning:", genes.shape)


# In[55]:


##==========================================
##  Merging all data
##==========================================
print("Merging datasets...")

# Merge all data into one master table
merged = (
    gene_pathway
    .merge(pathways, on="PATHWAY_ID", how="left")
    .merge(genes, on="GENE_ID", how="inner")
)

print(merged.head())


# In[56]:


##==========================================
##  Computing pathway and Gene Set
##==========================================

pathway_genes = (
    merged
    .groupby(["PATHWAY_ID", "PATHWAY_NAME"])["GENE_SYMBOL"]
    .apply(set)
    .reset_index()
)

print(pathway_genes.head())


# In[57]:


##==========================================
##  Computing crosstalk analysis 
##==========================================

crosstalk_results = []

# Compare every pair
for(path_id_1, path_name_1, genes1), (path_id_2, path_name_2, genes2) in combinations(pathway_genes.values, 2):
    overlap = genes1 & genes2
    if overlap:
        crosstalk_results.append({
            'PATHWAY_ID1': path_id_1,
            'PATHWAY_NAME1': path_name_1,
            'PATHWAY_ID2': path_id_2,
            'PATHWAY_NAME2': path_name_2,
            'NUMBER_OF_OVERLAPPING_GENES': len(overlap),
            'LIST_OF_OVERLAPPING_GENES': ";".join(sorted(overlap))
        })

## Save Results
crosstalk_df = pd.DataFrame(crosstalk_results)
crosstalk_df = crosstalk_df.sort_values(
    "NUMBER_OF_OVERLAPPING_GENES", ascending=False
)
crosstalk_df.to_csv("KEGG_crosstalk.csv", index=False)
print("Saved KEGG_crosstalk.csv")


# In[58]:


##==========================================
## Gene Centrality Analysis
##==========================================
# Count number of pathways per gene
gene_rank = (
    merged.groupby("GENE_SYMBOL")["PATHWAY_ID"]
    .nunique()
    .sort_values(ascending=False)
    .reset_index(name="NUM_PATHWAYS")
)
gene_rank.to_csv("gene_pathway_rank.csv", index=False)
print("Saved gene_pathway_rank.csv")


# In[59]:


##==========================================
## Shared Pathways for Top genes
##==========================================
top_genes = gene_rank["GENE_SYMBOL"].head(4).tolist()

# Retrieve sets of pathways for these genes
gene_to_pathways = (
    merged.groupby("GENE_SYMBOL")["PATHWAY_ID"]
    .apply(set)
)

# Find intersection
common = set.intersection(*[gene_to_pathways[g] for g in top_genes])
if not common:
    top_genes = top_genes[:3]
    common = set.intersection(*[gene_to_pathways[g] for g in top_genes])

print("Top genes:", top_genes)
#print("Common pathways:", common)


# In[60]:


##==========================================
## Venn Diagram
##==========================================

venn_data = {g: gene_to_pathways[g] for g in top_genes}
venn(venn_data)
plt.title("Pathway overlap among top genes")
plt.show()


# ## New Feature added -> JACCARD INDEX 
# 
# There are some KEGG pathways which are very larges and overlap with many other pathways because they contain many genes, this overlap doesn't always mean that they are functionally related. So, the question I would ask in this scenario is :
# 
# **Are pathway overlaps biologically meaningful or are they just driven by pathway size?**
# 
# An interesting solution which has been widely is **JACCARD SIMILARITY**
# The Jaccard index normalizes overlap by pathway size. This can be biologically interpreted as :
# 
# 1. Hight Overlap + High Jaccard mean -> there is true functional coupling between the genes.
# 2. High Overlap + Low Jaccard mean -> the overalp is size driven artifact.
# 3. Low Overlap + High Jaccard mean -> Small but tightly connected pathways.
# 

# In[70]:


##==========================================
## Jaccard-normalized pathway crosstalk
##==========================================

jaccard_results = []
for (path_id_1, path_name_1, genes1), (path_id_2, path_name_2, genes2) in combinations(pathway_genes.values, 2):
    intersection = genes1 & genes2
    union = genes1 | genes2

    if intersection:
        jaccard_results.append({
            "PATHWAY_ID1": path_id_1,
            "PATHWAY_NAME1": path_name_1,
            "PATHWAY_ID2": path_id_2,
            "PATHWAY_NAME2": path_name_2,
            "NUM_OVERLAPPING_GENES": len(intersection),
            "JACCARD_SCORE": len(intersection) / len(union),
            "OVERLAPPING_GENES": ";".join(sorted(intersection))
        })

## Creating the results table
jaccard_df = (
    pd.DataFrame(jaccard_results)
    .sort_values(
        by=["JACCARD_SCORE", "NUM_OVERLAPPING_GENES"],
        ascending=[False, False]
    )
)

print(jaccard_df.head(5))


# ##### In the table we see the top pathways which are fucntionally similar. Each row shows 2 pathways that share the largest proportion of genes.
# ##### For eg: Hypertophic Cardiomypathy and Dialated Cardiomypathy are clinically related diseases,
# ##### So we see a strong overlap which is expected biologically.
# ##### The table also shows numebr of overlapping genes .i.e. 89 genes in this case.
# ##### The jaccard score 0.774 shows that 77% of all genes in these pathways are same further inferring that there is tight functional coupling between the pathways 
# 

# In[74]:


##==============================================
## Visualizing through a netwrok lens
##==============================================
## Nodes -> KEGG pathways
## Edges-> Strong crosstalk/ High Jaccard Similarity)
## Dense Clusters(Hubs) -> Pathways involved in the same biological process


import networkx as nx
import matplotlib.pyplot as plt

## Filter pathway pairs by Jaccard threshold
## Keeps only strong, biologically meaningful connections
JACCARD_THRESHOLD = 0.5

network_df = jaccard_df[
    jaccard_df["JACCARD_SCORE"] >= JACCARD_THRESHOLD
]

## Build pathway crosstalk network
G = nx.Graph()

# Add edges
for _, row in network_df.iterrows():
    G.add_edge(
        row["PATHWAY_NAME1"],
        row["PATHWAY_NAME2"],
        weight=row["JACCARD_SCORE"]
    )

## Visualize the network

plt.figure(figsize=(12, 12))

# Layout: spring layout for biological networks
pos = nx.spring_layout(G, seed=42, k=0.4)

# Node sizes proportional to degree to make sure hub pathways stand out
node_sizes = [G.degree(n) * 300 for n in G.nodes()]

# Edge widths proportional to Jaccard similarity
edge_widths = [G[u][v]["weight"] * 4 for u, v in G.edges()]

# Draw nodes
nx.draw_networkx_nodes(G,pos,node_size=node_sizes,alpha=0.85)

# Draw edges
nx.draw_networkx_edges(G,pos,width=edge_widths,alpha=0.6)

# Draw labels
nx.draw_networkx_labels(G,pos,font_size=8)

plt.title(
    "KEGG Pathway Crosstalk Network\n"
    f"(Edges: Jaccard ≥ {JACCARD_THRESHOLD})"
)
plt.axis("off")
plt.show()

