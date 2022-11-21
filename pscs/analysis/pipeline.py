import numpy as np
from abc import ABC, abstractmethod
import anndata
from pscs.exceptions import PipeLineException, PreviousNodesNotRun, NodeRequirementsNotMet
from copy import deepcopy
import pandas as pd
import os

class Pipeline:
    """
    Class for storing, designing, and running related processes. Useful for keeping track of results and pipelines run.
    """
    def __init__(self, segment_list: list = None):
        self.pipeline = segment_list
        return

    def run(self):
        # identify top-level nodes, get data from those to pass along
        input_nodes = []
        for p in self.pipeline:
            if isinstance(p, InputNode):
                input_nodes.append(p)


        ann_data = pipes[0].run()
        for p in pipes[1:]:
            p.run(ann_data)
        return

    def reset(self):
        for p in self.pipeline:
            p.reset()
        return

class PipelineNode(ABC):
    """
    Class for describing an individual pipeline segment, with all necessary classes indicated.
    """
    def __init__(self):
        self.effect = []  # properties added to annotated data
        self.next = []  # list of nodes that follow this one
        self.previous = []  # list of nodes that leads to this node
        self.requirements = []  # list of effects that are expected to be completed before reaching this node
        self.inplace = False  # Whether the effect is done in the same data structure as the input (i.e., no copy)
        self.complete = False  # Whether node has been run
        self.node_type = None  # input, output, internal, pipeline
        self.result = None  # stored output
        return

    @abstractmethod
    def run(self,
            ann_dat: anndata.AnnData):
        """
        Method for executing this node's effects.
        Returns
        -------
        anndata.AnnData

        """
        return

    def _determine_self_type(self):
        """
        Determines which type of node this is.
        Returns
        -------
        None
        """
        len_next = len(self.next)
        len_prev = len(self.previous)
        if len_next == 0 and len_prev != 0:
            self.node_type = 'output'
        elif len_next != 0 and len_prev != 0:
            self.node_type = 'internal'
        elif len_next != 0 and len_prev == 0:
            self.node_type = 'input'
        else:
            self.node_type = 'pipeline'
        return

    @property
    def cumulative_effect(self) -> list:
        """
        Gets the effect of all preceding nodes, appends its own effect, and returns the list
        Returns
        -------
        list
            List of the cumulative effect of the pipeline
        """
        cumul = []
        for prev in self.previous:
            cumul += prev.cumulative_effect
        return cumul + self.effect

    @property
    def cumulative_requirements(self) -> list:
        """
        Gets the requirements of all preceding nodes, appends its own requirements, and returns the list. This is most
        useful when compared with a node's cumulative effect.
        Returns
        -------
        list
            List of all requirements needed up to this point
        """
        requirements = []
        for prev in self.previous:
            requirements += prev.cumulative_requirements
        return requirements + self.requirements

    @property
    def result(self):
        # Check if we have multiple outputs
        if len(self.next) > 1:
            # Return a copy of the data to prevent cross-contamination
            return deepcopy(self._result)
        else:
            return self._result

    @result.setter
    def result(self, value):
        self._result = value
        return

    def connect_to_output(self, node):
        """
        Connects the output of this node to the input of the specified node.
        Parameters
        ----------
        node : PipelineNode
            Node whose input should be connected.
        Returns
        -------
        None
        """
        self.next.append(node)
        node.previous.append(self)
        return

    def connect_to_input(self, node):
        """
        Connects the input of this node to the output of the specified node.
        Parameters
        ----------
        node : PipelineNode
            Node whose output should be connected

        Returns
        -------
        None
        """
        self.previous.append(node)
        node.next.append(node)
        return

    def reset(self):
        """Resets the output value of the node."""
        self.result = None
        self.complete = False
        return

    def validate_inputs(self) -> bool:
        """
        Checks that the inputs have been run and that they satisfy this node's requirements.
        Returns
        -------
        bool
            True if input is valid; raises exception otherwise

        Raises
        ------
        PreviousNodesNotRun
            If the nodes leading to this node have not been run and don't hold an output.
        NodeRequirementsNotMet
            If the effect of all nodes leading to this node do not produce the required effect.
        """
        # Check that input nodes have been run
        for inp in self.previous:
            if not inp.complete:
                raise PreviousNodesNotRun()
        # Check that cumulative effects of inputs meet this node's requirements
        cumul_effect = set()
        for inp in self.previous:
            cumul_effect.update(inp.cumulative_effect())
        req_set = set(self.requirements)
        unmet_reqs = req_set.difference(cumul_effect)
        if len(unmet_req) > 0:
            # Not all requirements have been met
            raise NodeRequirementsNotMet(unmet_reqs=list(unmet_reqs), reqs=list(req_set))
        return True

    def _terminate(self,
                   result=None):
        """
        Runs commands to finish the run.
        Parameters
        ----------
        result
            Result of the node to store.

        Returns
        -------
        None
        """
        self.result = result
        self.complete = True
        return


class InputNode(PipelineNode):
    """
    Input node prototype; serves to indicate that the node loads data from disk.
    """
    def connect_to_input(self, node):
        raise ValueError(f"This node doesn't receive input.")


class OutputNode(PipelineNode):
    """
    Output node prototype; serves to indicate that the node produces a file to disk.
    """
    def connect_to_output(self, node):
        raise ValueError(f"This node doesn't have an output to be received.")

class DataIntegrationNode(PipelineNode):
    """
    Combines
    """
    def run(self):
        # Combine results of inputs
        self.validate_inputs()
        # Result at 0 taken to be main
        base_data = self.previous[0].result
        base_effect = set(self.previous[0].cumulative_effect)

        # Merge others
        for inp in self.previous[1:]:
            inp_effect = set(inp.cumulative_effect)
            # Get the difference
            diff_effect = inp_effect.difference(base_effect)

def parse_node_effect(eff: str) -> tuple:
    """
    Examines the input string and identifies which property of AnnData was modified.
    Parameters
    ----------
    eff : str
        String containing the effect

    Returns
    -------
    tuple
        Tuple containing
    """
    return ()


class CSVLoadingNode(InputNode):
    def __init__(self,
                 path: str,
                 index_col: list = None):
        """
        Loads a .csv and stores the data in an AnnData object.
        Parameters
        ----------
        csv_path : str
            Path to the .csv to load
        """
        super().__init__()
        if path.endswith('.csv'):
            self.sep = ','
        elif path.endswith('.tsv'):
            self.sep = '\t'
        else:
            raise ValueError(f"File {os.path.basename(path)} must be either .csv or .tsv")
        self.path = path
        self.index_col = index_col
        self.effect = ['+X', '+obs', '+var']  # creates the data object
        return

    def run(self):
        # Load csv
        data = pd.read_csv(self.path, sep=self.sep, index_col=self.index_col)
        self._terminate(result=anndata.AnnData(data))
        return