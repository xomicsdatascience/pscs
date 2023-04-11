from os.path import join, basename, dirname
from pscs.analysis.pipeline.base import PipelineNode
from pscs.analysis.pipeline.base import InputNode, OutputNode
from pscs.analysis.pipeline.base import Pipeline
from werkzeug.utils import secure_filename
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


def load_from_nodes(node_json: str) -> Pipeline:
    """
    Loads a pipeline and its parameters from a file. Intended to be paired with the pipeline export from the website.
    Parameters
    ----------
    node_json : str
        Path to the node JSON exported from the designer.

    Returns
    -------
    Pipeline
    """
    # Load data
    f = open(node_json, 'r')
    pipeline = json.load(f)
    node_data = pipeline['nodes']
    f.close()
    node_dict = {}
    src_dict = {}
    dst_dict = {}
    for node in node_data:
        node_module = node['module']
        node_name = node['procName']
        # Restrict imports to PSCS pipeline:
        module = import_module(f'pscs.analysis.pipeline.{node_module}', package=__package__)
        node_class = inspect.getmembers(module, lambda mem: inspect.isclass(mem) and mem.__name__ == node_name)[0][1]
        # Instantiate the class with specified parameters
        node_instance = node_class(**node['paramsValues'])
        node_instance.nodeId = node['nodeId']
        node_num, node_srcs, node_dsts = identify_connections(node)
        node_dict[node_num] = node_instance
        src_dict[node_num] = node_srcs
        dst_dict[node_num] = node_dsts
    connect_nodes(node_dict, src_dict)
    return node_dict


def initialize_pipeline(node_json: str, input_files: dict, output_dir: str):
    """
    Starts a pipeline by loading it from a file and setting the input files specified for the relevant nodes.
    Parameters
    ----------
    node_json : str
        Path to the node JSON exported from the designer.
    input_files : dict
        Dict keyed by the node id with values corresponding to the id_data of the file to fetch from the database.
    output_dir : str
        Leading path where files should be placed.
    Returns
    -------

    """
    pipeline_nodes = load_from_nodes(node_json)
    pipeline = Pipeline(pipeline_nodes)
    assign_inputs(pipeline, input_files=input_files)
    assign_outputs(pipeline, output_dir)
    return pipeline


def assign_outputs(pipeline: Pipeline, output_dir: str) -> None:
    """
    Adds "output_dir" as prefix to "output_name" attributes of pipeline nodes.
    Parameters
    ----------
    pipeline : Pipeline
        Pipeline whose outputs should be assigned.
    output_dir : str
        Directory where to put the outputs.

    Returns
    -------
    None
    """
    for node_num, node in pipeline.pipeline.items():
        if isinstance(node, OutputNode):
            node.output_name = join(output_dir, secure_filename(node.output_name))
    return

def assign_inputs(pipeline: Pipeline, input_files: dict, path_keyword: str = 'path') -> None:
    """
    Assigns the input files
    Parameters
    ----------
    pipeline : Pipelime
        Pipeline whose input nodes should be assigned inputs.
    input_files : dict
        Dict keyed by the node id with values corresponding to the id_data of the file to fetch from the database.
    path_keyword : str
        Keyword for the input path.

    Returns
    -------
    None
    """
    input_keys = input_files.keys()
    for node_num, node in pipeline.pipeline.items():
        if node.nodeId in input_keys:
            node.__dict__[path_keyword] = input_files[node.nodeId]
    return


def connect_nodes(node_dict: dict, src_dict: dict) -> None:
    """
    Connects the nodes in node_dict according to the sources in src_dict.
    Parameters
    ----------
    node_dict : dict
        Dictionary of node instances keyed by their ID.
    src_dict : dict
        Dictionary of lists keyed by a node ID. Each ID in the list represents a node that has the key node as a source.

    Returns
    -------
    None
    """
    for key_node, node in node_dict.items():
        for srcnode in src_dict[key_node]:
            node.connect_to_output(node_dict[srcnode])
    return

def identify_connections(node: dict) -> (str, list, list):
    """
    Identifies the node's ID and input & output connections.
    Parameters
    ----------
    node : dict
        Dictionary describing a single node and its connections, as exported from the pipeline designer.

    Returns
    -------
    str
        Node ID
    list
        List of node IDs that have this node as their source.
    list
        List of node IDs that have this node as their destination.
    """
    node_id_src = -1
    node_id_dst = -1
    srcs = []
    dsts = []
    for s in node["srcConnectors"]:
        ssplit = s.split('-')
        srcs.append(ssplit[2])
        node_id_src = ssplit[1]
    for d in node['dstConnectors']:
        dsplit = d.split('-')
        dsts.append(dsplit[1])
        node_id_dst = dsplit[2]
    if node_id_dst != -1:
        node_id = node_id_dst
    elif node_id_src != -1:
        node_id = node_id_src
    else:
        raise ValueError('Could not determine ID for node.')
    return node_id, srcs, dsts



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
