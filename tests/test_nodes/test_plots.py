# These tests check whether the nodes run when provided with appropriate input.
# Given that the hash on any given plot might change for entirely cosmetic reasons, we only check that the files
# are created.
import unittest
import os
import math
from os.path import join
from pscs.analysis.pipeline.base import Pipeline
from pscs.analysis.pipeline.loaders import CSVLoadingNode
from pscs.analysis.pipeline import plotting
import tempfile
import scanpy as sc
import pathlib
test_data_path = os.path.join(os.path.dirname(__file__), "..", "test_data")
sc_figdir = os.path.join("/tmp", "scanpytest")
sc._settings.ScanpyConfig.figdir = pathlib.Path(sc_figdir)


class PlottingTests(unittest.TestCase):
    def test_Scatter(self):
        data_path = os.path.join(test_data_path, "metadata", "sample_csv.csv")
        inp = CSVLoadingNode(path=data_path)
        pl = plotting.Scatter(x="Edbrg6tpHfoR32tJGGS3",
                              y="0yZNMrBo5t1tNczQsSCf",
                              save=".png")
        inp.connect_to_output(pl)
        nodes = {"inp": inp, "pl": pl}
        pipe = Pipeline(nodes)
        pipe.run()
        out = join(sc_figdir, "scatter.png")
        self.assertTrue(os.path.exists(out))
        self.assertTrue(os.stat(out).st_size > 0)
        return

    def test_HeatMap(self):
        data_path = join(test_data_path, "metadata", "sample_csv.csv")
        inp = CSVLoadingNode(path=data_path)
        pl = plotting.HeatMap(var_names="Edbrg6tpHfoR32tJGGS3",
                              groupby="test",
                              save=".png")
        inp.connect_to_output(pl)
        inp.run()
        n = inp.result.shape[0]
        labels = ["0"]*math.floor(n/2) + ["1"]*math.ceil(n/2)
        inp.result.obs["test"] = labels
        pl.run()
        out = join(sc_figdir, "heatmap.png")
        print(sc_figdir)
        self.assertTrue(os.path.exists(out))
        self.assertTrue(os.stat(out).st_size > 0)
        return

    def test_DotPlot(self):
        data_path = join(test_data_path, "metadata", "sample_csv.csv")
        inp = CSVLoadingNode(path=data_path)
        pl = plotting.DotPlot(var_names="Edbrg6tpHfoR32tJGGS3",
                              groupby="test",
                              save=".png")
        inp.connect_to_output(pl)
        inp.run()
        n = inp.result.shape[0]
        labels = ["0"] * math.floor(n / 2) + ["1"] * math.ceil(n / 2)
        inp.result.obs["test"] = labels
        pl.run()
        out = join(sc_figdir, "dotplot_.png")
        self.assertTrue(os.path.exists(out))
        self.assertTrue(os.stat(out).st_size > 0)
        return

