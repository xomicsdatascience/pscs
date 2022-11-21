from .pipeline import PipelineNode, OutputNode
from anndata import AnnData
import pandas as pd
import numpy as np
import scanpy as sc
import pathlib
from matplotlib import pyplot as plt
import json
import plotly

class Plotter(OutputNode):
    def __init__(self,
                 plot_fn: callable,
                 plot_fn_kwargs: dict,
                 save_path: str):
        """
        Used to plot data according to the input plotting function.
        Parameters
        ----------
        plot_fn : callable
            Function with signature plot_fn(AnnData, **plot_fn_kwargs)
        plot_fn_kwargs : dict
            Dictionary containing the arguments to pass to the plotting function
        """
        super().__init__()
        self.plot_fn = plot_fn
        self.plot_fn_kwargs = plot_fn_kwargs
        self.save_path = save_path
        return

    def run(self,
            ann_data: AnnData):
        self.plot_fn(ann_data, **self.plot_fn_kwargs)
        plt.savefig(self.save_path)
        return


class ScanpyPlotter(OutputNode):
    def __init__(self,
                 plot_fn: callable,
                 plot_fn_kwargs: dict,
                 save_dir: str):
        """
        Same as Plotter, but with consideration for fixing Scanpy's weird file saving behaviour.
        Parameters
        ----------
        plot_fn : callable
            Function with signature plot_fn(AnnData, **plot_fn_kwargs)
        plot_fn_kwargs : dict
            Dictionary containing the arguments to pass to the plotting function
        save_dir : str
            Path to directory where to save data.
        """
        super().__init__()
        self.plot_fn = plot_fn
        self.plot_fn_kwargs = plot_fn_kwargs
        self.save_dir = save_dir
        return

    def run(self,
            ann_data: AnnData):
        sc._settings.ScanpyConfig.figdir = pathlib.Path(self.save_dir)
        sc._settings.ScanpyConfig.autoshow = False
        self.plot_fn(ann_data, **self.plot_fn_kwargs)
        return


class InteractivePlot(OutputNode):
    def __init__(self,
                 plot_fn: callable,
                 ann_data_parameter: str,
                 save_dir: str,
                 plot_fn_kwargs: dict = None
                 ):
        """
        Interactive node intended to interact with the JavaScript of the website. User input is likely to require
        re-running the node multiple times.
        Parameters
        ----------
        plot_fn : callable
            Function to call for plotting. The expected signature is plot_fn(pandas.DataFrame, **kwargs)
        ann_data_parameter : str
            Parameter from the AnnData object that should be used to plot (e.g. obsm, layer, etc)
        save_dir : str
            Directory where to save the output.
        plot_fn_kwargs : dict
            The keyword arguments to use for the plotting function.
        """
        self.plot_fn = plot_fn
        self.ann_data_parameter = ann_data_parameter
        self.plot_fn_kwargs = plot_fn_kwargs
        self.save_dir = save_dir
        return

    def run(self, ann_data: AnnData) -> str:
        """
        Executes
        Parameters
        ----------
        ann_data : AnnData
            Annotated data to be plotted.
        Returns
        -------
        str
            JSON object as string representing the plot, typically passed to a webpage's JavaScript for display.
        """
        data = get_dataframe_from_anndata(ann_data=ann_data,
                                          data_parameter=self.ann_data_parameter)
        fig = self.plot_fn(data, **self.plot_fn_kwargs)
        self.result = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        self.complete = True
        return

    def rerun(self, ann_data: AnnData,
              plot_fn_kwargs: dict = None):
        self.reset()
        if plot_fn_kwargs is None:
            return self.run(ann_data=ann_data)
        else:
            self.plot_fn_kwargs = plot_fn_kwargs
            return self.run(ann_data=ann_data)


def get_dataframe_from_anndata(ann_data: AnnData,
                               data_parameter: str):
    """
    Extracts the relevant array from the AnnData object and returns it. This is useful for allowing the original data
    array, 'X' to be addressed the same way as the different layers.
    Parameters
    ----------
    ann_data : AnnData
        AnnData object from which to extract the Pandas DataFrame
    data_parameter : str
        Parameter of AnnData to convert to DataFrame, e.g. 'X', "layer['raw']"
    Returns
    -------
    pd.DataFrame
        DataFrame of the specified AnnData parameter
    """
    if data_parameter == 'X':
        return ann_data.to_df()
    elif data_parameter.startswith('layer'):
        # Parse things
        data_param = data_parameter[len('layer')+2:-2]  # Remove square brackets, quotes
        return ann_data.to_df(layer=data_param)
    else:
        raise ValueError(f"Data parameter, {data_parameter} not recognized.")
