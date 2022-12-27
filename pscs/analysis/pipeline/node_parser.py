from os.path import join, basename, dirname
# from base import PipelineNode, InputNode, OutputNode
from pscs.analysis.pipeline.base import PipelineNode
from pscs.analysis.pipeline.base import InputNode, OutputNode
import os
import json
from importlib import import_module
import inspect
import re
exclude_files = set(["__init__.py", "base.py", "exceptions.py", basename(__file__)])
reserved_params = set(['num_inputs', 'num_outputs', 'effect'])  # These are parameters that should be public but unused

# This file parses available nodes and creates a JSON object describing coarse node properties
# for use by the pipeline designer.


def without_leading_underscore(d: dict) -> list:
    """
    Parses the keys of a dict and returns a list of those without leading underscores.
    Parameters
    ----------
    d : dict
        Dict whose key should be parsed

    Returns
    -------
    list
        List of keys without leading underscores.
    """
    without_underscore = []
    for k in d.keys():
        if not k.startswith('_'):
            without_underscore.append(k)
    return without_underscore


def get_node_parameters(node: callable) -> dict:
    """
    Extracts the relevant node parameters from the input node.
    Parameters
    ----------
    node : callable
        Reference to a PipelineNode (or subclass) from which to extract the parameters.
    Returns
    -------
    dict
        Dictionary containing the relevant parameters.
    """
    d = dict()
    param_dict = inspect.signature(node).parameters
    params = parse_params(param_dict)

    # To support MIMO nodes in the designer, this part needs to be updated.
    # Check which type of node this is
    d['num_inputs'] = 1
    d['num_outputs'] = 1
    if issubclass(node, InputNode):
        d['num_inputs'] = 0
    elif issubclass(node, OutputNode):
        d['num_outputs'] = 0
    d['parameters'] = params
    return d


def parse_params(params_dict: dict) -> dict:
    """
    Parses the param dict and returns a dict holding only the annotation (type) and default value.
    Parameters
    ----------
    params_dict : dict
        Dict of parameters as retruend by inspect.signature(callable).parameters

    Returns
    -------
    dict
        Dict of tuples of the form (annotation, default_value)
    """
    params = {}
    for param_name, param_value in params_dict.items():
        annot = str(param_value.annotation)
        annot = re.search("'(.*)'", annot).group()[1:-1]  # look for the string between quotes, then omit the quotes
        params[param_name] = (annot, param_value.default)
    return params


def find_unique_name(d: dict, name: str) -> str:
    """
    Creates a string from "name" that is not present in "d" by appending "_X".
    Parameters
    ----------
    d : dict
        Dictionary into which "name" is trying to fit uniquely.
    name : str
        Base name

    Returns
    -------
    str
        New name that is not in d
    """
    dkeys = d.keys()
    id = -1
    newname = name
    while newname in dkeys:
        id += 1
        newname = f"{name}_{id}"
    return newname


def main():
    self_path = dirname(__file__)
    tmp_files = os.listdir(self_path)
    node_files = []
    # Go through the other files in this file's directory, grab all .py that aren't excluded
    for f in tmp_files:
        if f.endswith('.py') and f not in exclude_files:
            node_files.append(f[:-len('.py')])  # remove the .py to be able to load directly
    js_dict = {}
    for module in node_files:  # iterate through the pipeline files
        imp = import_module(module, package=__package__)
        node_classes = inspect.getmembers(imp, lambda mem: inspect.isclass(mem) and mem.__module__ == module)
        for node_tuple in node_classes:  # iterate through the classes in the pipeline file
            node_name = node_tuple[0]
            node_name = find_unique_name(js_dict, node_name)  # in case nodes are named the same
            node_params = get_node_parameters(node_tuple[1])
            node_params['module'] = module
            js_dict[node_name] = node_params

    # find static dir
    idx = self_path.rfind('pscs')
    static_json_path = join(self_path[:idx + len('pscs')], 'static', 'node_data.json')
    f = open(static_json_path, 'w')
    json.dump(js_dict, f, indent=1)
    f.close()
    return


if __name__ == "__main__":
    main()
