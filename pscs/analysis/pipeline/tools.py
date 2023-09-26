from pscs.analysis.pipeline.base import PipelineNode
import pandas as pd
import numpy as np
import pathlib
import scanpy as sc
from scanpy import preprocessing as pp
from scanpy import tools as tl
from anndata import AnnData
from typing import Collection, Optional, Literal, Union, Sequence, Tuple
import math


class PCA(PipelineNode):
    important_parameters = ["n_comps", "zero_center"]

    def __init__(self,
                 n_comps: Optional[int] = None,
                 zero_center: Optional[bool] = True,
                 svd_solver: Literal["arpack", "randomized", "auto", "lobpcg"] = "arpack",
                 random_state: Union[None, int] = 0,
                 use_highly_variable: Optional[bool] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        tl.pca(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class TSNE(PipelineNode):
    important_parameters = ["n_pcs", "perplexity"]

    def __init__(self,
                 n_pcs: Optional[int] = None,
                 use_rep: Optional[str] = None,
                 perplexity: Union[float, int] = 30,
                 metric: str = "euclidean",
                 early_exaggeration: Union[float, int] = 12,
                 learning_rate: Union[float, int] = 1000,
                 random_state: Union[None, int] = 0):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        tl.tsne(**self.parameters)
        self._terminate(ann_data)
        return


class UMAP(PipelineNode):
    important_parameters = None

    def __init__(self,
                 min_dist: float = 0.5,
                 spread: float = 1.0,
                 n_components: int = 2,
                 maxiter: Optional[int] = None,
                 alpha: float = 1.0,
                 gamma: float = 1.0,
                 negative_sample_rate: int = 5,
                 init_pos: Union[Literal['paga','spectral', 'random'], None] = "spectral",
                 random_state: Union[None, int] = 0,
                 a: Optional[float] = None,
                 b: Optional[float] = None,
                 method: Literal["umap", "rapids"] = "umap",
                 neighbors_key: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        tl.umap(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DrawGraph(PipelineNode):
    important_parameters = None

    def __init__(self,
                 layout: Literal['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa'] = "fa",
                 root: Optional[int] = None,
                 random_state: Union[None, int] = 0,
                 key_added_ext: Optional[str] = None,
                 init_pos: Union[str, bool, None] = None,
                 neighbors_key: Optional[str] = None,
                 obsp: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        tl.umap(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DiffMap(PipelineNode):
    important_parameters = ["n_comps"]

    def __init__(self,
                 n_comps: int = 15,
                 neighbors_key: Optional[str] = None,
                 random_state: Union[None, int] = 0):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        tl.diffmap(ann_data, **self.parameters)
        return


class EmbeddingDensity(PipelineNode):
    important_parameters = None

    def __init__(self,
                 basis: str = "umap",
                 groupby: Optional[str] = None,
                 key_added: Optional[str] = None,
                 components: Union[str, Sequence[str], None] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        tl.embedding_density(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Leiden(PipelineNode):
    important_parameters = ["resolution"]

    def __init__(self,
                 resolution: float = 1,
                 random_state: Union[None, int] = 0,
                 restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
                 key_added: str = "leiden",
                 directed: bool = True,
                 use_weights: bool = True,
                 n_iterations: int = 1,
                 neighbors_key: Optional[str] = None,
                 obsp: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        tl.leiden(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Louvain(PipelineNode):
    important_parameters = ["resolution"]

    def __init__(self,
                 resolution: Optional[float] = None,
                 random_state: Union[None, int] = 0,
                 restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
                 key_added: str = "louvain",
                 flavor: Literal["vtraag", "igraph", "rapids"] = "vtraag",
                 directed: bool = True,
                 use_weights: bool = False,
                 neighbors_key: Optional[str] = None,
                 obsp: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())

    def run(self):
        ann_data = self._previous[0].result
        tl.louvain(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Dendrogram(PipelineNode):
    important_parameters = ["groupby", "n_pcs"]

    def __init__(self,
                 groupy: Collection[str],
                 n_pcs: Optional[int] = None,
                 use_rep: Optional[str] = None,
                 var_names: Optional[Collection[str]] = None,
                 use_raw: Optional[bool] = None,
                 cor_method: Literal['pearson', 'kendall', 'spearman'] = 'pearson',
                 linkage_method: str = "complete",
                 optimal_ordering: bool = False,
                 key_added: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)

    def run(self):
        ann_data = self._previous[0].result
        tl.dendrogram(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DPT(PipelineNode):
    important_parameters = ["n_dcs", "n_branchings", "min_group_size"]

    def __init__(self,
                 n_dcs: int = 10,
                 n_branchings: int = 0,
                 min_group_size: float = 0.01,
                 allow_kendall_tau_shift: bool = True,
                 neighbors_key: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())

    def run(self):
        ann_data = self._previous[0].result
        tl.dpt(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PAGA(PipelineNode):
    important_parameters = ["groups"]

    def __init__(self,
                 groups: Optional[str] = None,
                 use_rna_velocity: bool = False,
                 model: Literal["v1.2", "v1.0"] = "v1.2",
                 neigbors_key: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        tl.paga(ann_data, **self.parameters)
        self._terminate(ann_data)
        return
