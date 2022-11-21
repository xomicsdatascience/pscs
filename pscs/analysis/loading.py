from .pipeline import PipelineNode
from anndata import AnnData
import pandas as pd
import numpy as np

class CSVLoader(PipelineNode):
    def __init__(self,
                 file_path: str,
                 sep: str = ',',
                 index_col: str = 0,
                 dtype=np.float32):
        super().__init__()
        self.file_path = file_path
        self.sep = sep
        self.index_col = index_col
        self.dtype = dtype
        self.effect = ['X', 'obs', 'var']
        return

    def run(self) -> AnnData:
        if self.complete:
            return
        data = pd.read_csv(self.file_path, sep=self.sep, index_col=self.index_col)
        self.result = AnnData(data, dtype=self.dtype)
        self.complete = True
        return

    def connect_to_input(self, node):
        raise ValueError('This node can''t receive input.')