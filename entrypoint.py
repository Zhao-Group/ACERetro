"""
This is called by MMLI-backend, relevant PR here: https://github.com/moleculemaker/mmli-backend/pull/56

It takes input and returns output via Minio.
"""

import argparse
import pickle
import numpy as np
import io
import sys
import joblib
import six
import json
import time
# Keep these lines before any other imports (need to override defaults)
sys.modules['sklearn.externals.joblib'] = joblib
sys.modules['sklearn.externals.six'] = six

from sfscore import SFScore
from pathway_search_standalone.scripts.search_utils import hybridSearch, build_graph_from_async
import networkx as nx
from minio import Minio
from minio.error import S3Error
import json
import os

def score_from_smi(smiles):
    PROJECT_DIRECTORY = '/app/ACERetro'
    sfscore_path = PROJECT_DIRECTORY + '/sfscore'
    sys.path.append(str(sfscore_path))

    sfscore_model = SFScore()
    sfscore_model.load()
    sfscore = sfscore_model.score_from_smi(smiles)
    print('SMILES:', smiles, ', SFScore:', sfscore)
    return sfscore

def get_chemoenzy_path_async(smiles, max_depth=10, chem_topk=10, max_num_templates=250, max_branching=15, time_lim=180):
    hybrid_search = hybridSearch()
    explored_rxns, explored_nodes, start_node = hybrid_search.get_chemoenzy_path_async(
        smiles, max_depth=max_depth, chem_topk=chem_topk, max_num_templates=max_num_templates,
        max_branching=max_branching, time_lim=time_lim
    )
    return explored_rxns, explored_nodes, start_node

def build_graph_from_async_wrapper(explored_rxns, explored_nodes, start_node):
    graph_json = build_graph_from_async(explored_rxns, explored_nodes, start_node)
    g = nx.node_link_graph(graph_json)

    graph_data = {}
    for gg in g.nodes():
        graph_data[gg] = g.nodes[gg]
        print(gg)
        print(g.nodes[gg])
    return graph_data

def add_smiles_svgs_to_json(json_graph):
    """
    Observation: 36 of 66 objects in the graph dictionary do not contain a smiles string. And it happens to be that the "whole first half" have it and "the 2nd half" don't have it.
    """
    from draw_chemical_svg import draw_chemical_svg
    errors = 0
    total_count = 0

    for index, data in json_graph['graph'].items():
        total_count += 1
        try: 
            smile = data['smiles']
            # print("Building SVG for smile: ", smile)
            svg = draw_chemical_svg(smile)
            json_graph['graph'][index]['smiles_svg'] = svg
        except Exception as e: 
            errors += 1
            # print("ERROR drawing svg for SMILE:", data.get('smiles', 'No SMILES provided'))
    # print("Total errors: ", errors)
    # print("Total count: ", total_count)
    return json_graph

def download_json_from_minio(bucket: str ='', path: str=''):
    client = Minio(
        os.environ['MINIO_URL'],  
        access_key=os.environ['MINIO_ACCESS_KEY'],
        secret_key=os.environ['MINIO_SECRET_ACCESS_KEY'],
        secure=True
    )

    try:
        # Get the object from Minio
        response = client.get_object(bucket, path)
        # Read the data from the response
        data = response.read()
        # Parse the JSON data
        return json.loads(data)
    except S3Error as e:
        print("Error occurred.", e)
        raise

def upload_json_to_minio(data, bucket: str ='', path: str=''):
    client = Minio(
        os.environ['MINIO_URL'],  
        access_key=os.environ['MINIO_ACCESS_KEY'],
        secret_key=os.environ['MINIO_SECRET_ACCESS_KEY'],
        secure=True
    )

    json_data = json.dumps(convert_to_json_serializable(data)).encode('utf-8')
    client.put_object(bucket, path, data=io.BytesIO(json_data), length=len(json_data), content_type='application/json')
    print(f"Successfully uploaded results to {bucket}{path}")

def numpy_to_python(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, dict):
        return {k: numpy_to_python(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [numpy_to_python(item) for item in obj]
    else:
        return obj

def convert_to_json_serializable(data):
    return json.loads(json.dumps(data, default=numpy_to_python))

def main():
    parser = argparse.ArgumentParser(description="Run various chemical analysis functions.")
    parser.add_argument('--job_id', required=True, help='Minio path to the JSON file containing input parameters')
    args = parser.parse_args()

    params = None
    try: 
        # Download the input JSON file from Minio
        params = download_json_from_minio(bucket='aceretro', path=f"/{args.job_id}/in/input.json")
    except Exception as e: 
        print(f"xxx Failed to download input from Minio with error: {e}")
        print("WARNING: FALLING BACK TO PLACEHOLDER DEFAULT SMILES STRING... `O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO`")
        params = {
            'smiles': 'O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO',
        }

    results = {}
    start_time = time.monotonic()
    if 'smiles' not in params:
        parser.error("xxx Error: The input JSON file must contain 'smiles' for score_from_smi")
    else: 
        results['sfscore'] = score_from_smi(params['smiles'])
        print(f"[1/3] Completed score_from_smi(). Runtime: {(time.monotonic() - start_time):.2f} seconds")
    
    start_time = time.monotonic()
    if 'smiles' not in params:
        parser.error("xxx Error: The input JSON file must contain 'smiles' for get_chemoenzy_path_async")
    else: 
        explored_rxns, explored_nodes, start_node = get_chemoenzy_path_async(
            params['smiles'], params.get('max_depth', 10), params.get('chem_topk', 10), 
            params.get('max_num_templates', 250), params.get('max_branching', 15), 
            params.get('time_lim', 20)
        )
        results['explored_rxns'] = explored_rxns
        results['explored_nodes'] = explored_nodes
        results['start_node'] = start_node
        print(f"[2/3] Completed get_chemoenzy_path_async(). Runtime: {(time.monotonic() - start_time):.2f} seconds")

    start_time = time.monotonic()
    if 'explored_rxns' not in results or 'explored_nodes' not in results or 'start_node' not in results:
        parser.error("xxx Error: The results dictionary must contain 'explored_rxns', 'explored_nodes', and 'start_node' for build_graph_from_async")
    else: 
        results['graph'] = build_graph_from_async_wrapper(results['explored_rxns'], results['explored_nodes'], results['start_node'])
        print(f"[3/3] build_graph_from_async(). Runtime: {(time.monotonic() - start_time):.2f} seconds")
    

    # Generate SVGs of the Smiles strings
    start_time = time.monotonic()
    results = add_smiles_svgs_to_json(results)
    print(f"[4/3] add_smiles_svgs_to_json(). Runtime: {(time.monotonic() - start_time):.2f} seconds")
    
    # Upload the results to Minio
    try:
        upload_json_to_minio(results, bucket='aceretro', path=f"/{args.job_id}/out/output.json")
    except Exception as e: 
        print(f"Failed to upload results to Minio with error: {e}")

if __name__ == "__main__":
    main()