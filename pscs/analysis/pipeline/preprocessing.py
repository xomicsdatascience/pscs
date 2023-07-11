from pscs.analysis.pipeline.base import PipelineNode
import pandas as pd
import numpy as np
import pathlib
import scanpy as sc
from scanpy import preprocessing as pp
from anndata import AnnData
from typing import Collection, Optional, Literal, Union, Sequence
import math

class CalculateQCMetrics(PipelineNode):
    def __init__(self,
                 expr_type: str = "counts",
                 var_type: str = "genes",
                 qc_vars: Collection[str] = (),
                 percent_top: Optional[Collection[int]] = (50, 100, 200, 500),
                 layer: Optional[str] = None,
                 use_raw: bool = False,
                 log1p: bool = True):
        # Assign all arguments to self.parameters
        super().__init__()
        self.parameters = {}
        for param, value in vars().items():
            if param == "self":
                continue
            self.parameters[param] = value
        self.parameters["inplace"] = True
        self.effect = ["total_{var_type}_by_{expr_type}", "total_{expr_type}",
                       "pct_{expr_type}_in_top_{percent_top}_{var_type}",
                       "total_{expr_type}_{qc_var}", "pct_{expr_type}_{qc_var}"]
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.calculate_qc_metrics(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class FilterCells(PipelineNode):
    def __init__(self,
                 min_counts: Optional[int] = None,
                 min_genes: Optional[int] = None,
                 max_counts: Optional[int] = None,
                 max_genes: Optional[int] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        vars_without_self["inplace"] = True
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        # filter_cells only accepts one arg at a time; go through each one
        for par_key, par_val in self.parameters.items():
            par = {par_key: par_val}
            if par_val is not None:
                _ = sc.pp.filter_cells(ann_data, **par)
        self._terminate(ann_data)
        return
#
# class FilterCells(PipelineNode):
#     def __init__(self,
#                  minimum_unique_protein: int = None,
#                  maximum_unique_protein: int = None,
#                  minimum_total_protein: float = None,
#                  maximum_total_protein: float = None
#                  ):
#         """
#         Excludes cells (rows) from the dataset based on proteins/genes values. The "unique" protein parameters examine the
#          number of non-zero columns. The "total" protein parameters examine the sum across all columns.
#         Parameters
#         ----------
#         minimum_unique_protein : int
#             Minimum number of unique proteins for cells to be retained (equality is retained).
#         maximum_unique_protein : int
#             Maximum number of unique proteins for cells to be retained (equality is retained).
#         minimum_total_protein : float
#             Minimum value of the sum of all proteins in a cell for a cell to be retained (equality is retained).
#         maximum_total_protein : float
#             Maximum value of the sum of all proteins in a cell for a cell to be retained (equality is retained).
#         """
#         super().__init__()
#         vars_without_self = vars()
#         del vars_without_self["self"]
#         self.store_vars_as_parameters(**vars_without_self)
#         # self.minimum_unique_protein = None if minimum_unique_protein is None else int(minimum_unique_protein)
#         # self.maximum_unique_protein = None if maximum_unique_protein is None else int(maximum_unique_protein)
#         # self.minimum_total_protein = None if minimum_total_protein is None else float(minimum_total_protein)
#         # self.maximum_total_protein = None if maximum_total_protein is None else float(maximum_total_protein)
#         self.effect = ["~obs['obs_keep']"]
#         return
#
#     def run(self):
#         ann_data = self._previous[0].result
#         # get obs keep if it already exists
#         if 'obs_keep' in ann_data.obs.keys():
#             obs_keep = ann_data.obs['obs_keep']
#         else:
#             obs_keep = pd.Series([True]*ann_data.shape[0], index=ann_data.obs_names)
#         obs_keep_list = []
#         # filter_cells only accepts one arg at a time; go through each one
#         if self.minimum_unique_protein is not None:
#             ok, _ = sc.pp.filter_cells(ann_data, min_genes=float(self.minimum_unique_protein), inplace=False)
#             obs_keep_list.append(ok)
#         if self.maximum_unique_protein is not None:
#             ok, _ = sc.pp.filter_cells(ann_data, max_genes=float(self.maximum_unique_protein), inplace=False)
#             obs_keep_list.append(ok)
#         if self.minimum_total_protein is not None:
#             ok, _ = sc.pp.filter_cells(ann_data, min_counts=float(self.minimum_total_protein), inplace=False)
#             obs_keep_list.append(ok)
#         if self.maximum_total_protein is not None:
#             ok, _ = sc.pp.filter_cells(ann_data, max_counts=float(self.maximum_total_protein), inplace=False)
#             obs_keep_list.append(ok)
#         # Iterate through the kept observations
#         for ok in obs_keep_list:
#             obs_keep = obs_keep & ok
#         ann_data.obs['obs_keep'] = obs_keep  # only needed if 'obs_keep' wasn't in ann_data.obs
#
#         self._terminate(ann_data)
#         return


class FilterGenes(PipelineNode):
    def __init__(self,
                 min_counts: Optional[int] = None,
                 min_cells: Optional[int] = None,
                 max_counts: Optional[int] = None,
                 max_cells: Optional[int] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        vars_without_self["inplace"] = True
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        # filter_genes only accepts one arg at a time; go through each one
        for par_key, par_val in self.parameters.items():
            par = {par_key: par_val}
            if par_val is not None:
                _ = sc.pp.filter_genes(ann_data, **par)
        self._terminate(ann_data)
        return


# class FilterGenes(PipelineNode):
#     def __init__(self,
#                  minimum_cell_presence: int = None,
#                  maximum_cell_presence: int = None,
#                  minimum_cell_amount: float = None,
#                  maximum_cell_amount: float = None):
#         """
#         Excludes proteins/genes (columns) from the dataset based on the number of cells that contain them and the
#         quantities in which they are detected. The "cell_presence" parameters examine how many cells have non-zero
#         amounts of a particular protein/gene. The "cell_amount" parameters examine the total amount of a protein/gene
#         is detected across all cells.
#         Parameters
#         ----------
#         minimum_cell_presence : int
#             Minimum number of cells in which a protein/gene is found. Removes rare proteins/genes.
#         maximum_cell_presence : int
#             Maximum number of cell in which a protein/gene is found. Removes common proteins/genes.
#         minimum_cell_amount : float
#             Minimum total quantity across all cells for a protein/gene to be considered valid. Removes proteins/genes
#             that are found in overall low quantities.
#         maximum_cell_amount : float
#             Maximum total quantity across all cells for a protein/gene to be considered valid. Removes proteins/genes
#             that are found in overall high quantities.
#         """
#         super().__init__()
#         self.minimum_cell_presence = None if minimum_cell_presence is None else int(minimum_cell_presence)
#         self.maximum_cell_presence = None if maximum_cell_presence is None else int(maximum_cell_presence)
#         self.minimum_cell_amount = None if minimum_cell_amount is None else float(minimum_cell_amount)
#         self.maximum_cell_amount = None if maximum_cell_amount is None else float(maximum_cell_amount)
#         self.effect = ["~var['var_keep']"]
#         return
#
#     def run(self):
#         # filter_genes only accepts on arg at a time
#         ann_data = self._previous[0].result
#         if "var_keep" in ann_data.var.keys():
#             var_keep = ann_data.var["var_keep"]
#         else:
#             var_keep = pd.Series([True]*ann_data.n_vars, index=ann_data.var.index)
#         var_keep_list = []
#         if self.minimum_cell_presence is not None:
#             vk, _ = sc.pp.filter_genes(ann_data, min_cells=float(self.minimum_cell_presence), inplace=False)
#             var_keep_list.append(vk)
#         if self.maximum_cell_presence is not None:
#             vk, _ = sc.pp.filter_genes(ann_data, max_cells=float(self.maximum_cell_presence), inplace=False)
#             var_keep_list.append(vk)
#         if self.minimum_cell_amount is not None:
#             vk, _ = sc.pp.filter_genes(ann_data, min_counts=float(self.minimum_cell_amount), inplace=False)
#             var_keep_list.append(vk)
#         if self.maximum_cell_amount is not None:
#             vk, _ = sc.pp.filter_genes(ann_data, max_counts=float(self.maximum_cell_amount), inplace=False)
#             var_keep_list.append(vk)
#         # Iterate through kept variables
#         for vk in var_keep_list:
#             var_keep = var_keep & vk
#         ann_data.var["var_keep"] = var_keep
#
#         self._terminate(ann_data)
#         return


class HighlyVariableGenes(PipelineNode):
    def __init__(self,
                 layer: Optional[str] = None,
                 n_top_genes: Optional[int] = None,
                 min_mean: Optional[float] = 0.0125,
                 max_mea: Optional[float] = 3,
                 min_disp: Optional[float] = 0.5,
                 max_disp: Optional[float] = str(math.inf),
                 span: Optional[float] = 0.3,
                 n_bins: int = 20,
                 flavor: Literal['seurat', 'cell_ranger', 'seurat_v3'] = 'seurat',
                 subset: bool = False,
                 batch_key: Optional[str] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        self.parameters["inplace"] = True
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.highly_variable_genes(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Log1p(PipelineNode):
    def __init__(self,
                 base: Optional[float] = None,
                 chunked: Optional[bool] = None,
                 chunk_size: Optional[int] = None,
                 layer: Optional[str] = None,
                 obsm: Optional[str] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.log1p(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PCA(PipelineNode):
    def __init__(self,
                 n_comps: Optional[int] = None,
                 zero_center: Optional[bool] = True,
                 svd_solver: str = "arpack",
                 random_state: Union[None, int] = 0,
                 use_highly_variable: Optional[bool] = None,
                 dtype: str = 'float32',
                 chunked: bool = False,
                 chunk_size: Optional[int] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.pca(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class NormalizeTotal(PipelineNode):
    def __init__(self,
                 target_sum: Optional[float] = None,
                 exclude_highly_expressed: bool = False,
                 max_fraction: float = 0.05,
                 key_added: Optional[str] = None,
                 layer: Optional[str] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        vars_without_self["inplace"] = True
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.normalize_total(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RegressOut(PipelineNode):
    def __init__(self,
                 keys: Union[str, Sequence[str]] = None,
                 n_jobs: Optional[int] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        vars_without_self["inplace"] = True
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.regress_out(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Scale(PipelineNode):
    def __init__(self,
                 zero_center: bool = True,
                 max_value: Optional[float] = None,
                 layer: Optional[str] = None,
                 obsm: Optional[str] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.regress_out(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Subsample(PipelineNode):
    def __init__(self,
                 fraction: Optional[float] = None,
                 n_obs: Optional[int] = None,
                 random_state: Union[None, int] = 0):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.subsample(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DownsampleCounts(PipelineNode):
    def __init__(self,
                 counts_per_cell: Union[int, Collection[int], None] = None,
                 total_counts: Optional[int] = None,
                 random_state: Union[None, int] = 0,
                 replace: bool = False):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.downsample_counts(ann_data, **self.parameters)
        self._terminate(ann_data)
        return

