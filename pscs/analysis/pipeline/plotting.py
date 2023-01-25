from pscs.analysis.pipeline.base import OutputNode
import matplotlib
from matplotlib import pyplot as plt
import scanpy as sc
# matplotlib.use('Agg')

class ScatterPlot(OutputNode):
    def __init__(self,
                 x_gene: str = None,
                 y_gene: str = None,
                 output_name: str = None,
                 color: str = None):
        super().__init__()
        self.x = x_gene
        self.y = y_gene
        self.output_name = output_name
        self.color = color
        return

    def run(self):
        ann_data = self._previous[0].result
        plt.scatter(ann_data[:, self.x].X, ann_data[:, self.y].X)
        plt.xlabel(self.x)
        plt.ylabel(self.y)
        plt.savefig(self.output_name)  ## TODO: make sure this is sanitized before reaching this
        # sc.pl.scatter(ann_data, x=self.x, y=self.y, color=self.color)
        return

class HighestExpressed(OutputNode):
    def __init__(self):
        super().__init__()
        return

    def run(self):
        raise NotImplementedError
        return

class Violin(OutputNode):
    def __init__(self):
        super().__init__()
        return

    def run(self):
        raise NotImplementedError
        return