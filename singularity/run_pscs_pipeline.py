#!/bin/python3.11
from pscs_api import node_parser as panp
import pscs_api as pa
import json
import os
import argparse

from pscs_ops.io import AnnDataLoadingNode
from collections import defaultdict as dd

def main(node_json: dict,
         input_file_dict: dict):
    # Load and connect nodes; expected to have dangling nodes
    nodes = panp.load(node_json)
    partial_pipeline = pa.Pipeline(nodes)

    # For every node receiving an input from outside the container, create a loader
    input_mapping = dict
    ordered_input = dd(dict)
    for node_id, file_path in input_file_dict.items():
        loader_id = generate_unique_node_id(nodes)
        port_id = get_port_from_id(node_id)
        # ordered_input

        input_mapping[node_id] = loader_id
        loader = AnnDataLoadingNode(file_path)
        nodes[loader_id] = loader


def get_port_from_id(node_id) -> str:
    """Extracts the port ID from the input node id spec"""
    return node_id.split(".")[-1]


def generate_unique_node_id(nodes: dict) -> str:
    node_ids = set(nodes.keys())
    idx = 0
    s = f"temp_loader_{idx}"
    while s in node_ids:
        idx += 1
        s = f"temp_loader_{idx}"
    return s


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = 'PSCS Pipeline',
                                     description = 'Runs a PSCS pipeline')
    parser.add_argument('node_json')
    parser.add_argument('input_files_json')
    parser.add_argument('out_dir')
    args = parser.parse_args()
    
    f = open(args.input_files_json, 'r')
    input_files = json.load(f)
    f.close()
    
    os.makedirs(args.out_dir, exist_ok=True)
    pipeline = panp.initialize_pipeline(node_json=args.node_json,
                                               input_files=input_files,
                                               output_dir=args.out_dir)
    pipeline.run()
                                               
