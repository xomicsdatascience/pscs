import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from flask import current_app
import os
import pathlib

common_directory = './common/'

def analyze(filename: str,
            save_dir: str,
            exp_name: str,
            image_format: str = '.png') -> str:
    """
    Runs analysis.
    Parameters
    ----------
    filename : str
        File to analyze.
    save_dir : str
        Directory where to store data.
    exp_name : str
        Name of the experiment; used for naming output files.
    image_format : str
        List of image formats to save the figures as.
    Returns
    -------
    str
        Top-level directory where results are stored.
    """
    sc._settings.ScanpyConfig.figdir = pathlib.Path(save_dir)
    mt_thresh = 100
    save_suffix = f"_{exp_name}{image_format}"

    # First load data
    if filename.endswith('.tsv'):
        sep = '\t'
    else:
        sep = ','
    data = pd.read_csv(filename, sep=sep, header=0, index_col=0)
    ann_data = ad.AnnData(data, dtype=float)
    # Highest expressed genes
    sc.pl.highest_expr_genes(ann_data, n_top=30, save=save_suffix)

    # Filter out outlier cells
    sc.pp.filter_cells(ann_data, min_genes=0)
    sc.pp.filter_genes(ann_data, min_cells=0)

    # Identify MT genes
    f = open(os.path.join(common_directory, 'mt_list.txt'), 'r')
    mt_list = f.read().split(',')
    f.close()

    ann_data.var['mt'] = [c in mt_list for c in ann_data.var_names]
    sc.pp.calculate_qc_metrics(ann_data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(ann_data,
                 ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, save=save_suffix)

    # Scatter plots
    sc.pl.scatter(ann_data, x='total_counts', y='pct_counts_mt', title="Percent MT vs. Total Proteins", save=save_suffix)
    sc.pl.scatter(ann_data, x='total_counts', y='n_genes_by_counts', title='Number of genes expressed vs. protein count', save=save_suffix)

    # MT quantity thresholding; removing everything above threshold
    ann_data = ann_data[ann_data.obs.pct_counts_mt <= mt_thresh, :]

    # Normalize quantities
    sc.pp.normalize_total(ann_data, target_sum=1e7)

    # Apply log
    sc.pp.log1p(ann_data)
    sc.pp.highly_variable_genes(ann_data, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(ann_data, save=save_suffix)

    ann_data.raw = ann_data

    # UMAP!
    sc.pp.neighbors(ann_data, n_neighbors=12, n_pcs=20, use_rep='X')
    sc.tl.umap(ann_data)

    sc.pl.umap(ann_data, add_outline=True, save=save_suffix)
    sc.tl.leiden(ann_data)
    sc.pl.umap(ann_data, color=['leiden'], add_outline=True, save=save_suffix)

    # Automatically-identified markers
    max_genes = ann_data.var.shape[0]
    sc.tl.rank_genes_groups(ann_data, 'leiden', n_genes=min([max_genes, 20]), sharey=False, fontsize=16, save=save_suffix)

    plot_count = 4  # Number of markers to plot
    table_count = min([50,max_genes])  # number to store in table

    gene_rankings = list(ann_data.uns['rank_genes_groups']['names'])
    gene_rankings = [list(r) for r in gene_rankings]

    genes_to_plot = []
    for group_idx in range(len(gene_rankings[0])):
        for gene_idx in range(plot_count):
            genes_to_plot.append(gene_rankings[gene_idx][group_idx])

    sc.pl.umap(ann_data, color=['leiden','leiden','leiden','leiden', *genes_to_plot], add_outline=True, save=save_suffix)

    return