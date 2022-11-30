import numpy as np
import anndata
from pscs.analysis.pipeline.base import InputNode
from copy import deepcopy
import pandas as pd
import os


# This file contains abstract classes for nodes.
class CSVLoadingNode(InputNode):
    def __init__(self,
                 path: str = None,
                 index_col: list = None):
        """
        Loads a .csv and stores the data in an AnnData object.
        Parameters
        ----------
        csv_path : str
            Path to the .csv to load
        """
        super().__init__()
        path = str(path)
        if path.endswith('.tsv'):
            self.sep = '\t'
        else:
            self.sep = ','

        self.path = path
        self.index_col = index_col
        self.effect = ['+X', '+obs', '+var']  # creates the data object
        return

    def run(self):
        # Load csv
        data = pd.read_csv(self.path, sep=self.sep, index_col=self.index_col)
        self._terminate(result=anndata.AnnData(data))
        return

    def _validate(self, suppress=False):
        """
        Verifies that input parameters are valid; returns True if valid, and raises an exception if not. If
        suppress=True, returns False instead of raising an exception if parameters are invalid.
        Returns
        -------
        bool
            True if parameters are valid; False if parameters are invalid and suppress=True.
        """
        exception_str = []
        if not self.path.endswith('.tsv') or not self.path.endswith('.csv'):
            exception_str.append(f"File {os.path.basename(path)} must be either .csv or .tsv")
        if len(exception_str) > 0:
            raise ValueError(exception_str)
        return
