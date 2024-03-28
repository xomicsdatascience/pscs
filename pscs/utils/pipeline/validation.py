import json
from typing import Collection, Mapping
from collections import defaultdict as dd


def load_analysis(pipeline_json_file):
    with open(pipeline_json_file, 'r') as f:
        return json.load(f)["nodes"]


def validate_pipeline(nodelist: Collection[dict]) -> (bool, list):
    """
    Checks that every node in the list is valid.
    Parameters
    ----------
    nodelist : Collection[dict]
        List of nodes. Nodes are keyed with "nodeId" and "params".

    Returns
    -------
    bool:
        Whether the entire pipeline is valid.
    list:
        List of reasons why the
    list[dict]:
        List of invalid nodes with the reason that they are invalid.
    """

    invalid_nodes = dd(list)
    invalid_pipeline_reasons = []
    if not nodes_have_unique_ids(nodelist):
        invalid_pipeline_reasons.append("Nodes must have unique IDs")
    for node in nodelist:
        is_valid, invalid_reasons = validate_node(node)
        if not is_valid:
            invalid_nodes[node["nodeId"]].append(invalid_reasons)
    return len(invalid_nodes) == 0 and len(invalid_pipeline_reasons) == 0, invalid_pipeline_reasons, invalid_nodes


def validate_node(node: dict) -> (bool, list):
    validation_steps = [required_parameters_are_defined,
                        input_ports_receive_single_connection,
                        all_input_ports_receive_connection,
                        at_least_one_output_connected]
    invalid_reasons = []
    for v in validation_steps:
        if not v(node):
            invalid_reasons.append(v.__name__)
    return len(invalid_reasons) == 0, invalid_reasons


def nodes_have_unique_ids(nodelist: Collection[dict]):
    s = set()
    for node in nodelist:
        node_id = node["nodeId"]
        if node_id in s:
            return False
        s.add(node_id)
    return True


def required_parameters_are_defined(node: Mapping):
    """Verifies that all required parameters are defined."""
    if len(node["required_parameters"]) > 0:
        for param in node["required_parameters"]:
            if param not in node["paramsValues"] or node["paramsValues"][param] is None:
                return False
    return True


def input_ports_receive_single_connection(node: Mapping):
    """Verifies that input ports receive only one connection."""
    connected_ports = set()
    for dst_connector in node["dstConnectors"]:
        _, (_, dst_port) = parse_connector_id(dst_connector)
        if dst_port in connected_ports:
            return False
        connected_ports.add(dst_port)
    return True


def all_input_ports_receive_connection(node: Mapping):
    """Verifies that each input port has been connected."""
    connected_ports = set()
    for dst_connector in node["dstConnectors"]:
        _, (_, dst_port) = parse_connector_id(dst_connector)
        connected_ports.add(dst_port)
    return len(connected_ports) == node["num_inputs"]


def at_least_one_output_connected(node: Mapping):
    """Verifies that at least one of the node's outputs is connected."""
    return node["num_outputs"] == 0 or len(node["srcConnectors"]) > 0  # srcConnectors is the list of connections starting at this node


def parse_connector_id(connector_id: str) -> ((str, str), (str, str)):
    _, source, dst = connector_id.split("-")
    source_node, source_port = source.split(".")
    dst_node, dst_port = dst.split(".")
    return (source_node, source_port), (dst_node, dst_port)
