from pscs.analysis.pipeline.base import PipelineNode
import pandas as pd
import numpy as np
import pathlib
import scanpy as sc
from anndata import AnnData
from typing import Optional, Literal, Union

class Neighbors(PipelineNode):
    def __init__(self,
                 n_neighbors: int = 15,
                 n_pcs: Optional[int] = None,
                 use_rep: Optional[str] = None,
                 knn: bool = True,
                 random_state: Union[None, int] = 0,
                 method: Optional[Literal["umap", "gauss", "rapids"]] = "umap",
                 key_added: Optional[str] = None):
        """
        Computes the distances between the num_neighbours nearest neighbours.
        Parameters
        ----------
        num_neighbours : int
            Number of neighbours for which to return the distances.
        num_pcs : int
            Number of principal components to consider.

        """
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc.pp.neighbors(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class UMAP(PipelineNode):
    def __init__(self,
                 min_dist: float = 0.5,
                 spread: float = 1.0,
                 n_components: int = 2,
                 maxiter: Optional[int] = None,
                 alpha: float = 1.0,
                 gamma: float = 1.0,
                 negative_sample_rate: int = 5,
                 neighbors_key: Optional[str] = None):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc.tl.umap(ann_data, **self.parameters)
        self._terminate(ann_data)
        return

