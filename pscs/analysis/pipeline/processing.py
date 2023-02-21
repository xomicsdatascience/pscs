from pscs.analysis.pipeline.base import PipelineNode
import pandas as pd
import numpy as np
import pathlib
import scanpy as sc
from anndata import AnnData


class Neighbors(PipelineNode):
    def __init__(self,
                n_neighbors: int = 12,
                n_pcs: int = 20):
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
        self.n_neighbors = n_neighbors
        self.n_pcs = n_pcs
        self.effect = ["+.uns['neighbors']", "+.obsp['distances']", "+.obsp['connectivities']"]
        return

    def run(self):
        ann_data = self._previous[0].result
        sc.pp.neighbors(ann_data, n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)
        self._terminate(ann_data)
        return


class UMAP(PipelineNode):
    def __init__(self):
        super().__init__()
        self.requirements = ["+.uns['neighbors']", "+.obsp['distances']", "+.obsp['connectivities']"]
        self.effect = ["+.obsm['X_umap']"]
        return

    def run(self):
        ann_data = self._previous[0].result
        sc.tl.umap(ann_data)
        self._terminate(ann_data)
        return
