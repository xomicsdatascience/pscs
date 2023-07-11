from pscs.analysis.pipeline.base import OutputNode
import matplotlib
from matplotlib import pyplot as plt
import scanpy as sc
from scanpy import plotting as pl
from typing import Optional, Sequence, Literal, Collection, Union
import os
from pathlib import Path


class Scatter(OutputNode):
    def __init__(self,
                 x: Optional[str] = None,
                 y: Optional[str] = None,
                 color: Union[str, Collection[str], None] = None,
                 layers: Union[str, Collection[str], None] = None,
                 save: Union[str, bool] = "scatter.png"):
        super().__init__()
        vars_without_self = vars()
        save = "_" + str(save)
        if not save.endswith(".png") and not save.endswith(".svg") and not save.endswith(".jpeg"):
            save += ".png"
        vars_without_self["show"] = False
        vars_without_self["save"] = save
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        tmp_parameters = self.parameters.copy()
        tmp_parameters["save"] = os.path.basename(tmp_parameters["save"])
        pl.scatter(ann_data, **tmp_parameters)
        return


class HighestExpressedGenes(OutputNode):
    def __init__(self,
                 n_top: int = 30,
                 save: str = "highest_expressed.png",
                 log: bool = False):
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        vars_without_self["show"] = False
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = os.path.dirname(self.parameters["save"])
        tmp_parameters = self.parameters.copy()
        tmp_parameters["save"] = os.path.basename(tmp_parameters["save"])
        pl.highest_expr_genes(ann_data, **tmp_parameters)
        return


class UMAPPlot(OutputNode):
    def __init__(self,
                 add_outline: Optional[bool] = None,
                 color: Union[str, Sequence[str], None] = None,
                 save: str = "UMAP.png"):
        """
        Creates the UMAP plot based on the ann_data.obsm['X_umap'] matrix.
        Parameters
        ----------
        add_outline : bool
            Whether to add outlines to clusters on the plot.
        color : list
            List of column names to use to apply colors.
        save : str
            Name to use for the output file.
        """
        super().__init__()
        vars_without_self = vars()
        del vars_without_self["self"]
        self.store_vars_as_parameters(**vars_without_self)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = os.path.dirname(self.parameters["save"])
        tmp_parameters = self.parameters.copy()
        tmp_parameters["save"] = os.path.basename(tmp_parameters["save"])
        sc.pl.umap(ann_data, **tmp_parameters)
        return
