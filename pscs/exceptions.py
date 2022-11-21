# This file is for custom exceptions specific to the project

class PipeLineException(Exception):
    pass


class PreviousNodesNotRun(PipeLineException):
    def __init__(self, err: str = None):
        """
        Exception used to indicate that a particular node can't run.
        Parameters
        ----------
        err : str
            String giving additional context to the exception.
        """
        msg = "The previous nodes in the pipeline have not been run"
        if err is None:
            msg += "."
        else:
            msg += f": {err}"
        super().__init__(msg)
        return


class NodeRequirementsNotMet(PipeLineException):
    def __init__(self,
                 unmet_reqs: list = None,
                 reqs: list = None):
        """
        Exception used to indicate that a node's requirements are not met; previous nodes do not produce the expected
        modifications to the AnnData object.
        Parameters
        ----------
        unmet_reqs : list
            List of requirements that have not been met.
        reqs : list
            List of requirements for the node.
        """
        msg = "The requirements of this node have not been met"
        if unmet_reqs is None or reqs is None:
            msg += "."
        else:
            msg += f": {unmet_reqs} not in {reqs}."
        super().__init__(msg)
        return
