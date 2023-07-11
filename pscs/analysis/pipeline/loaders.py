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
                 index_col: str = None):
        """
        Loads a .csv and stores the data in an AnnData object.
        Parameters
        ----------
        path : str
            Path to the .csv to load
        index_col : [str, int]
            Variable indicating the column to use as an index. If string, looks for the column with the header.
            If int, takes the column at that position.
        """
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        if isinstance(path, str) and path.endswith('.tsv'):
            self.parameters["sep"] = '\t'
        else:
            self.parameters["sep"] = ','
        return

    def run(self):
        # Load csv
        data = pd.read_csv(self.parameters["path"], sep=self.parameters["sep"], index_col=self.parameters["index_col"])
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
        if not self.parameters["save"].endswith('.tsv') or not self.parameters["save"].endswith('.csv'):
            exception_str.append(f"File {os.path.basename(self.parameters['save'])} must be either .csv or .tsv")
        if len(exception_str) > 0:
            raise ValueError(exception_str)
        return
