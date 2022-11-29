import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from flask import current_app, g, session
import os
import pathlib
import plotly.express as px
from pscs.db import get_db
import uuid
import hashlib
common_directory = './common/'
import matplotlib
matplotlib.use('Agg')


def analyze(filename: str,
            save_dir: str,
            exp_name: str,
            image_format: str = '.svg'):
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
    annotated_data
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
    out_path = os.path.join(save_dir, f'highest_expr_genes{save_suffix}')
    register_result(out_path, 'graph', title='Highest-expressed genes', description="Ranking of the most-expressed genes")

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
    out_path = os.path.join(save_dir, f"violin{save_suffix}")
    register_result(out_path, 'graph', title="Violin plot of gene counts and quantities", description="Violin plots showing the distribution of the number of genes expressed, the total protein counts, and the fraction of proteins from the mitochondria.")
    # Scatter plots
    sc.pl.scatter(ann_data, x='total_counts', y='pct_counts_mt', title="Percent MT vs. Total Proteins", save=f"_percentmt{save_suffix}")
    out_path = os.path.join(save_dir, f'scatter_percentmt{save_suffix}')
    register_result(out_path, 'graph', 'Fraction MT', description='Fraction of MT vs. total proteins')
    sc.pl.scatter(ann_data, x='total_counts', y='n_genes_by_counts', title='Number of genes expressed vs. protein count', save=f"_genes{save_suffix}")
    out_path = os.path.join(save_dir, f'scatter_genes{save_suffix}')
    register_result(out_path, 'graph', 'Genes expressed and total protein count', description='The plot shows the number of genes expressed related to the total protein count.')

    # MT quantity thresholding; removing everything above threshold
    ann_data = ann_data[ann_data.obs.pct_counts_mt <= mt_thresh, :]

    # Normalize quantities
    sc.pp.normalize_total(ann_data, target_sum=1e7)

    # Apply log
    sc.pp.log1p(ann_data)
    sc.pp.highly_variable_genes(ann_data, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(ann_data, save=f"_log{save_suffix}")
    out_path = os.path.join(save_dir, f'filter_genes_dispersion_log{save_suffix}')
    register_result(out_path, 'graph', 'Highly-variable genes', description='The graph shows which genes/proteins are highly variable across the dataset')
    ann_data.raw = ann_data

    # UMAP!
    sc.pp.neighbors(ann_data, n_neighbors=12, n_pcs=20, use_rep='X')
    sc.tl.umap(ann_data)

    sc.pl.umap(ann_data, add_outline=True, save=save_suffix)
    sc.tl.leiden(ann_data)
    sc.pl.umap(ann_data, color=['leiden'], add_outline=True, save=save_suffix)
    out_path = os.path.join(save_dir, f'umap{save_suffix}')
    register_result(out_path, 'graph', 'UMAP', description='The UMAP shows a two-dimensional reduction of the proteome while attempting to preserve the high-dimensional structure.')

    # Automatically-identified markers
    max_genes = ann_data.var.shape[0]
    sc.tl.rank_genes_groups(ann_data, 'leiden', n_genes=min([max_genes, 20]), sharey=False, fontsize=16, save=save_suffix)
    # out_path = os.path.join(save_dir, f'rank_genes_groups{save_suffix}')
    # register_result(out_path, 'graph', 'Rank Genes by Groups', description='The plots show which genes best differentiate the clusters from one another.')

    plot_count = 4  # Number of markers to plot
    table_count = min([50,max_genes])  # number to store in table

    gene_rankings = list(ann_data.uns['rank_genes_groups']['names'])
    gene_rankings = [list(r) for r in gene_rankings]

    genes_to_plot = []
    for group_idx in range(len(gene_rankings[0])):
        for gene_idx in range(plot_count):
            genes_to_plot.append(gene_rankings[gene_idx][group_idx])

    sc.pl.umap(ann_data, color=['leiden','leiden','leiden','leiden', *genes_to_plot], add_outline=True, save=f"_markers_{save_suffix}")
    out_path = os.path.join(save_dir, f"umap_markers_{save_suffix}")
    register_result(out_path, 'graph', title='Automatically-identified markers', description="Marker genes overlaid on the UMAP, demonstrating which genes are useful for demarking each cluster.")

    umap_points = ann_data.obsm['X_umap']
    px.scatter(umap_points[:,0], umap_points[:,1])

    return ann_data

def register_result(filename: str,
                    result_type: str,
                    title: str,
                    description: str):
    """
    SQL wrapper for registering a result into the results table.
    Parameters
    ----------
    filename : str
        Name of the file
    result_type : str
        Type of result (e.g., 'graph', 'table')
    title : str
        Title of the result, appropriate as a header.
    description : str
        Full description of the result.

    Returns
    -------
    None
    """
    db = get_db()
    id_result = get_hash(filename)
    id_user = g.user['id_user']
    id_project = session['CURRENT_PROJECT']
    id_analysis = 0  # TODO

    db.execute('INSERT INTO results (id_result, id_project, id_analysis, file_path, result_type, title, description) VALUES (?,?,?,?,?,?,?)', (id_result, id_project, id_analysis, filename, result_type, title, description))
    db.commit()
    return

def get_hash(filename: str):
    f = open(filename, 'rb')
    d = f.read(1024)
    s = hashlib.sha256()
    while d:
        s.update(d)
        d = f.read(1024)
    return s.hexdigest()