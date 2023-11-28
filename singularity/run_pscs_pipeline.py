#!/bin/python3.11
from pscs_api import node_parser
import json
import os
import argparse

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
    pipeline = node_parser.initialize_pipeline(node_json=args.node_json,
                                               input_files=input_files,
                                               output_dir=args.out_dir)
    pipeline.run()
                                               
