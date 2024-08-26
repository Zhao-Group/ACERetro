"""
This is called by MMLI-backend, relevant PR here: https://github.com/moleculemaker/mmli-backend/pull/56

It takes input and returns output via Minio.
"""

import argparse
import io
import sys
import joblib
import six
import json
import time
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
    sys.modules['sklearn.externals.joblib'] = joblib
    sys.modules['sklearn.externals.six'] = six

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

def download_json_from_minio(bucket: str ='', path: str=''):
    # Initialize the Minio client
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
    # Initialize the Minio client
    client = Minio(
        os.environ['MINIO_URL'],  
        access_key=os.environ['MINIO_ACCESS_KEY'],
        secret_key=os.environ['MINIO_SECRET_ACCESS_KEY'],
        secure=True
    )

    try:
        # Convert data to JSON
        json_data = json.dumps(data).encode('utf-8')
        # Upload the JSON data to Minio
        client.put_object(bucket, path, data=io.BytesIO(json_data), length=len(json_data), content_type='application/json')
        print(f"Successfully uploaded results to {bucket}/{path}")
    except S3Error as e:
        print("Error occurred while uploading.", e)
        raise

def main():
    parser = argparse.ArgumentParser(description="Run various chemical analysis functions.")
    # parser.add_argument('function', choices=['score_from_smi', 'get_chemoenzy_path_async', 'build_graph_from_async'],
    #                     help='Function to run')
    parser.add_argument('--job_id', required=True, help='Minio path to the JSON file containing input parameters')
    args = parser.parse_args()

    # Download the JSON file from Minio
    # TODO: Get env vars for Minio
    # params = download_json_from_minio(bucket='aceretro', path=f"/{args.job_id}/in/input.json")

    params = {
        'smiles': 'O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO',
    }

    results = {}
    start_time = time.monotonic()
    if 'smiles' not in params:
        parser.error("❌ Error: The input JSON file must contain 'smiles' for score_from_smi")
    else: 
        results['sfscore'] = score_from_smi(params['smiles'])
        print(f"[1/3] Completed score_from_smi(). ⏰ Runtime: {(time.monotonic() - start_time):.2f} seconds")
    
    start_time = time.monotonic()
    if 'smiles' not in params:
        parser.error("❌ Error: The input JSON file must contain 'smiles' for get_chemoenzy_path_async")
    else: 
        explored_rxns, explored_nodes, start_node = get_chemoenzy_path_async(
            params['smiles'], params.get('max_depth', 10), params.get('chem_topk', 10), 
            params.get('max_num_templates', 250), params.get('max_branching', 15), 
            params.get('time_lim', 180)
        )
        results['explored_rxns'] = explored_rxns
        results['explored_nodes'] = explored_nodes
        results['start_node'] = start_node
        print(f"[2/3] Completed get_chemoenzy_path_async(). ⏰ Runtime: {(time.monotonic() - start_time):.2f} seconds")

    start_time = time.monotonic()
    if 'explored_rxns' not in results or 'explored_nodes' not in results or 'start_node' not in results:
        parser.error("❌ Error: The results dictionary must contain 'explored_rxns', 'explored_nodes', and 'start_node' for build_graph_from_async")
    else: 
        results['graph'] = build_graph_from_async_wrapper(results['explored_rxns'], results['explored_nodes'], results['start_node'])
        print(f"[3/3] build_graph_from_async(). ⏰ Runtime: {(time.monotonic() - start_time):.2f} seconds")
    
    # Upload the results to Minio
    # TODO: Get env vars for Minio
    # upload_json_to_minio(results, bucket='aceretro', path=f"/{args.job_id}/out/output.json")

if __name__ == "__main__":
    main()