from pscs.analysis.pipeline.base import PipelineNode
import pandas as pd
import numpy as np
import pathlib
import scanpy as sc
from anndata import AnnData


class FilterCellValues(PipelineNode):
    def __init__(self,
                 minimum_unique_protein: int = None,
                 maximum_unique_protein: int = None,
                 minimum_total_protein: float = None,
                 maximum_total_protein: float = None
                 ):
        """
        Excludes cells (rows) from the dataset based on proteins/genes values. The "unique" protein parameters examine the
         number of non-zero columns. The "total" protein parameters examine the sum across all columns.
        Parameters
        ----------
        minimum_unique_protein : int
            Minimum number of unique proteins for cells to be retained (equality is retained).
        maximum_unique_protein : int
            Maximum number of unique proteins for cells to be retained (equality is retained).
        minimum_total_protein : float
            Minimum value of the sum of all proteins in a cell for a cell to be retained (equality is retained).
        maximum_total_protein : float
            Maximum value of the sum of all proteins in a cell for a cell to be retained (equality is retained).
        """
        super().__init__()
        self.minimum_unique_protein = minimum_unique_protein
        self.maximum_unique_protein = maximum_unique_protein
        self.minimum_total_protein = minimum_total_protein
        self.maximum_total_protein = maximum_total_protein
        self.effect = ["~obs['obs_keep']"]
        return

    def run(self):
        ann_data = self.previous[0].result
        # get obs keep if it already exists
        if 'obs_keep' in ann_data.obs.keys():
            obs_keep = ann_data.obs['obs_keep']
        else:
            obs_keep = pd.Series([True]*ann_data.shape[0], index=ann_data.obs_names)
        obs_keep_list = []
        # filter_cells only accepts one arg at a time; go through each one
        if self.minimum_unique_protein is not None:
            ok, _ = sc.pp.filter_cells(ann_data, min_genes=self.minimum_unique_protein, inplace=False)
            obs_keep_list.append(ok)
        if self.maximum_unique_protein is not None:
            ok, _ = sc.pp.filter_cells(ann_data, max_genes=self.maximum_unique_protein, inplace=False)
            obs_keep_list.append(ok)
        if self.minimum_total_protein is not None:
            ok, _ = sc.pp.filter_cells(ann_data, min_counts=self.minimum_total_protein, inplace=False)
            obs_keep_list.append(ok)
        if self.maximum_total_protein is not None:
            ok, _ = sc.pp.filter_cells(ann_data, max_counts=self.maximum_total_protein, inplace=False)
            obs_keep_list.append(ok)
        # Iterate through the kept observations
        for ok in obs_keep_list:
            obs_keep = obs_keep & ok
        ann_data.obs['obs_keep'] = obs_keep  # only needed if 'obs_keep' wasn't in ann_data.obs

        self._terminate(ann_data)
        return

class FilterGeneValues(PipelineNode):
    def __init__(self,
                 minimum_cell_presence: int = None,
                 maximum_cell_presence: int = None,
                 minimum_cell_amount: float = None,
                 maximum_cell_amount: float = None):
        """
        Excludes proteins/genes (columns) from the dataset based on the number of cells that contain them and the
        quantities in which they are detected. The "cell_presence" parameters examine how many cells have non-zero
        amounts of a particular protein/gene. The "cell_amount" parameters examine the total amount of a protein/gene
        is detected across all cells.
        Parameters
        ----------
        minimum_cell_presence : int
            Minimum number of cells in which a protein/gene is found. Removes rare proteins/genes.
        maximum_cell_presence : int
            Maximum number of cell in which a protein/gene is found. Removes common proteins/genes.
        minimum_cell_amount : float
            Minimum total quantity across all cells for a protein/gene to be considered valid. Removes proteins/genes
            that are found in overall low quantities.
        maximum_cell_amount : float
            Maximum total quantity across all cells for a protein/gene to be considered valid. Removes proteins/genes
            that are found in overall high quantities.
        """
        super().__init__()
        self.minimum_cell_presence = minimum_cell_presence
        self.maximum_cell_presence = maximum_cell_presence
        self.minimum_cell_amount = minimum_cell_amount
        self.maximum_cell_amount = maximum_cell_amount
        self.effect = ["~var['var_keep']"]
        return

    def run(self):
        # filter_genes only accepts on arg at a time
        ann_data = self.previous[0].result
        if "var_keep" in ann_data.var.keys():
            var_keep = ann_data.var["var_keep"]
        else:
            var_keep = pd.Series([True]*ann_data.n_vars, index=ann_data.var.index)
        var_keep_list = []
        if self.minimum_cell_presence is not None:
            vk, _ = sc.pp.filter_genes(ann_data, min_cells=self.minimum_cell_presence, inplace=False)
            var_keep_list.append(vk)
        if self.maximum_cell_presence is not None:
            vk, _ = sc.pp.filter_genes(ann_data, max_cells=self.maximum_cell_presence, inplace=False)
            var_keep_list.append(vk)
        if self.minimum_cell_amount is not None:
            vk, _ = sc.pp.filter_genes(ann_data, min_counts=self.minimum_cell_amount, inplace=False)
            var_keep_list.append(vk)
        if self.maximum_cell_amount is not None:
            vk, _ = sc.pp.filter_genes(ann_data, max_counts=self.maximum_cell_amount, inplace=False)
            var_keep_list.append(vk)
        # Iterate through kept variables
        for vk in var_keep_list:
            var_keep = var_keep & vk
        ann_data.var["var_keep"] = var_keep

        self._terminate(ann_data)
        return
