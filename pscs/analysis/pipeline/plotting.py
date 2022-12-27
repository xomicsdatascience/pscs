from pscs.analysis.pipeline.base import OutputNode
import matplotlib
from matplotlib import pyplot as plt
import scanpy as sc
matplotlib.use('Agg')

class ScatterPlot(OutputNode):
    def __init__(self,
                 x: str = None,
                 y: str = None,
                 color: str = None):
        super().__init__()
        self.x = x
        self.y = y
        self.color = color
        return

    def run(self):
        ann_data = self.previous[0].result
        sc.pl.scatter(ann_data, x=self.x, y=self.y, color=self.color)
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