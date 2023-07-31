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
    important_parameters = ["expr_type", "var_type"]

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
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.calculate_qc_metrics(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class FilterCells(PipelineNode):
    important_parameters = None

    def __init__(self,
                 min_counts: Optional[int] = None,
                 min_genes: Optional[int] = None,
                 max_counts: Optional[int] = None,
                 max_genes: Optional[int] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
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


class FilterGenes(PipelineNode):
    important_parameters = None

    def __init__(self,
                 min_counts: Optional[int] = None,
                 min_cells: Optional[int] = None,
                 max_counts: Optional[int] = None,
                 max_cells: Optional[int] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
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


class HighlyVariableGenes(PipelineNode):
    important_parameters = ["n_top_genes"]

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
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.highly_variable_genes(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Log1p(PipelineNode):
    important_parameters = None

    def __init__(self,
                 base: Optional[float] = None,
                 chunked: Optional[bool] = None,
                 chunk_size: Optional[int] = None,
                 layer: Optional[str] = None,
                 obsm: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.log1p(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PCA(PipelineNode):
    important_parameters = ["n_comps", "zero_center"]

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
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.pca(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class NormalizeTotal(PipelineNode):
    important_parameters = ["target_sum", "exclude_highly_expressed"]

    def __init__(self,
                 target_sum: Optional[float] = None,
                 exclude_highly_expressed: bool = False,
                 max_fraction: float = 0.05,
                 key_added: Optional[str] = None,
                 layer: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.normalize_total(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RegressOut(PipelineNode):
    important_parameters = ["keys"]

    def __init__(self,
                 keys: Union[str, Sequence[str]] = None,
                 n_jobs: Optional[int] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.regress_out(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Scale(PipelineNode):
    important_parameters = ["zero_center"]

    def __init__(self,
                 zero_center: bool = True,
                 max_value: Optional[float] = None,
                 layer: Optional[str] = None,
                 obsm: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.regress_out(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Subsample(PipelineNode):
    important_parameters = None

    def __init__(self,
                 fraction: Optional[float] = None,
                 n_obs: Optional[int] = None,
                 random_state: Union[None, int] = 0):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.subsample(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DownsampleCounts(PipelineNode):
    important_parameters = ["counts_per_cell", "total_counts"]

    def __init__(self,
                 counts_per_cell: Union[int, Collection[int], None] = None,
                 total_counts: Optional[int] = None,
                 random_state: Union[None, int] = 0,
                 replace: bool = False):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.downsample_counts(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Combat(PipelineNode):
    important_parameters = ["covariates"]

    def __init__(self,
                 key: str = "batch",
                 covariates: Optional[Collection[str]] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.combat(ann_data, **self.parameters, inplace=True)
        self._terminate(ann_data)
        return


class Neighbors(PipelineNode):
    important_parameters = ["n_neighbors", "n_pcs"]

    def __init__(self,
                 n_neighbors: int = 15,
                 n_pcs: Optional[int] = None,
                 use_rep: Optional[str] = None,
                 knn: bool = True,
                 random_state: Union[None, int] = 0,
                 method: Optional[Literal["umap", "gauss", "rapids"]] = "umap",
                 metric: Literal['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan', 'braycurtis', 'canberra',
                                 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis',
                                 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener',
                                 'sokalsneath', 'sqeuclidean', 'yule'] = "euclidean",
                 key_added: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.neighbors(ann_data, **self.parameters)
        self._terminate(ann_data)
        return
