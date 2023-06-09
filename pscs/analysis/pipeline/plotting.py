from pscs.analysis.pipeline.base import OutputNode
import matplotlib
from matplotlib import pyplot as plt
import scanpy as sc
# matplotlib.use('Agg')


class ScatterPlot(OutputNode):
    def __init__(self,
                 x_gene: str = None,
                 y_gene: str = None,
                 output_name: str = "scatter.png"):
        super().__init__()
        self.x = x_gene
        self.y = y_gene
        self.output_name = output_name
        return

    def run(self):
        ann_data = self._previous[0].result
        plt.figure()
        plt.scatter(ann_data[:, self.x].X, ann_data[:, self.y].X)
        plt.xlabel(self.x)
        plt.ylabel(self.y)
        plt.savefig(self.output_name)
        plt.close()
        return


class HighestExpressedGenes(OutputNode):
    def __init__(self,
                 n_top: int = 30,
                 output_name: str = "highest_expressed.png"):
        """
        Fraction of counts assigned to each gene over all cells.
        Computes, for each gene, the fraction of counts assigned to that gene within a cell. The n_top genes with the highest mean fraction over all cells are plotted as boxplots.
        Parameters
        ----------
        n_top : int
            Number of top-expressed genes to plot. Default: 30.
        log : bool
            Whether to plot the x-axis in log scale. Default: False.
        output_name : str
            Name of the output file.
        """
        super().__init__()
        self.n_top = n_top
        self.output_name = output_name
        return

    def run(self):
        ann_data = self._previous[0].result
        sc.pl.highest_expr_genes(ann_data, n_top=self.n_top, show=False)
        plt.savefig(self.output_name)
        return


class UMAPPlot(OutputNode):
    def __init__(self,
                 add_outline: bool = None,
                 color: list = None,
                 output_name: str = "UMAP.png"):
        """
        Creates the UMAP plot based on the ann_data.obsm['X_umap'] matrix.
        Parameters
        ----------
        add_outline : bool
            Whether to add outlines to clusters on the plot.
        color : list
            List of column names to use to apply colors.
        output_name : str
            Name to use for the output file.
        """
        super().__init__()
        self.add_outline = add_outline
        self.output_name = output_name
        self.requirements = ["+.obsm['X_umap']"]
        self.color = color
        return

    def run(self):
        ann_data = self._previous[0].result
        sc.pl.umap(ann_data, add_outline=self.add_outline)
        plt.savefig(self.output_name)
        return
