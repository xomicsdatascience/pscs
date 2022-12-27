from base import PipelineNode
import pandas as pd
import numpy as np
import pathlib
import scanpy as sc
from anndata import AnnData


class Neighbours(PipelineNode):
    def __int__(self,
                num_neighbours: int = 12,
                num_pcs: int = 20):
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
        self.num_neighbours = num_neighbours
        self.num_pcs = num_pcs
        return

    def run(self):
        raise NotImplementedError
        # ann_data = self.previous[0].result
        return


class UMAP(PipelineNode):
    def __init__(self):
        super().__init__()
        return

    def run(self):
        raise NotImplementedError
        return