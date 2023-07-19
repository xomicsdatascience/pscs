from abc import ABC, abstractmethod
from pscs.exceptions import PreviousNodesNotRun, NodeRequirementsNotMet
from copy import deepcopy
from ast import literal_eval


class PipelineNode(ABC):
    """
    Class for describing an individual pipeline segment, with all necessary classes indicated.
    """
    def __init__(self):
        self.num_inputs = 1
        self.num_outputs = 1
        self.parameters = {}
        self._effect = []  # properties added to annotated data
        self._next = []  # list of nodes that follow this one
        self._previous = []  # list of nodes that lead to this node
        self._requirements = []  # list of effects that are expected to be completed before reaching this node
        self._inplace = False  # Whether the _effect is done in the same data structure as the input (i.e., no copy)
        self._result = None  # stored output
        return

    @abstractmethod
    def run(self):
        """
        Method for executing this node's effects.
        """
        return

    def run_pipeline(self):
        if self.is_ready:
            self.run()
            for next_node in self._next:
                next_node.run_pipeline()
        return

    def store_vars_as_parameters(self, **kwargs):
        """Stores the key-value pairs in self.parameters. Excludes "self"""
        for param, value in kwargs.items():
            if param == "self" or param == "return" or param.startswith("_"):
                continue
            if value is None or value == "null":
                self.parameters[param] = None
            else:
                self.parameters[param] = value

    @property
    def cumulative_effect(self) -> list:
        """
        Gets the _effect of all preceding nodes, appends its own _effect, and returns the list
        Returns
        -------
        list
            List of the cumulative _effect of the pipeline
        """
        cumul = []
        for prev in self._previous:
            cumul += prev.cumulative_effect
        return cumul + self._effect

    @property
    def cumulative_requirements(self) -> list:
        """
        Gets the _requirements of all preceding nodes, appends its own _requirements, and returns the list. This is most
        useful when compared with a node's cumulative _effect.
        Returns
        -------
        list
            List of all _requirements needed up to this point
        """
        requirements = []
        for prev in self._previous:
            requirements += prev.cumulative_requirements
        return requirements + self._requirements

    @property
    def result(self):
        # Check if we have multiple outputs
        if len(self._next) > 1:
            # Return a copy of the data to prevent cross-contamination
            return deepcopy(self._result)
        else:
            return self._result

    @result.setter
    def result(self, value):
        self._result = value
        return

    @property
    def is_complete(self) -> bool:
        return self._result is not None

    @property
    def is_ready(self) -> bool:
        """Checks whether all inputs to this node have a result."""
        for p in self._previous:
            if not p.is_complete:
                return False
        return True

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
        self._next.append(node)
        node._previous.append(self)
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
        self._previous.append(node)
        node._next.append(node)
        return

    def reset(self):
        """Resets the output value of the node."""
        self.result = None
        return

    def validate_inputs(self) -> bool:
        """
        Checks that the inputs have been run and that they satisfy this node's _requirements.
        Returns
        -------
        bool
            True if input is valid; raises exception otherwise

        Raises
        ------
        PreviousNodesNotRun
            If the nodes leading to this node have not been run and don't hold an output.
        NodeRequirementsNotMet
            If the _effect of all nodes leading to this node do not produce the required _effect.
        """
        # Check that input nodes have been run
        for inp in self._previous:
            if not inp.is_complete:
                raise PreviousNodesNotRun()
        # Check that cumulative effects of inputs meet this node's _requirements
        cumul_effect = set()
        for inp in self._previous:
            cumul_effect.update(inp.cumulative_effect())
        req_set = set(self._requirements)
        unmet_reqs = req_set.difference(cumul_effect)
        if len(unmet_reqs) > 0:
            # Not all _requirements have been met
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
        return


class InputNode(PipelineNode):
    """
    Input node prototype; serves to indicate that the node loads data from disk.
    """
    def __init__(self):
        super().__init__()
        self.num_inputs = 0

    def connect_to_input(self, node):
        raise ValueError(f"This node doesn't receive input.")


class OutputNode(PipelineNode):
    """
    Output node prototype; serves to indicate that the node produces a file to disk.
    """
    def __init__(self):
        super().__init__()
        self.num_outputs = 0

    def connect_to_output(self, node):
        raise ValueError(f"This node doesn't have an output to be received.")


class Pipeline:
    def __init__(self, nodes: dict = None):
        """
        Class for storing, designing, and running related processes. Useful for keeping track of results and pipelines run.
        Parameters
        ----------
        nodes : dict
            Nodes indexed by their ID. Default: None.
        """
        self.pipeline = nodes

    def run(self):
        # Get how many are ready, how many are done
        ready_list = []
        for _, node in self.pipeline.items():
            if node.is_ready:
                ready_list.append(node)
        while len(ready_list) > 0:
            run_list = ready_list
            ready_list = []
            for node in run_list:
                node.run()
                for next_node in node._next:
                    if next_node.is_ready:
                        ready_list.append(next_node)
        return

    def reset(self):
        for p in self.pipeline:
            p.reset()
        return
