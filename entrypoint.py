"""
This is called by MMLI-backend, relevant PR here: https://github.com/moleculemaker/mmli-backend/pull/56

It takes input and returns output via Minio.
"""

import argparse
import sys
import joblib
import six
import json
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

def download_json_from_minio(minio_path):
    # Initialize the Minio client
    client = Minio(
        "play.min.io",  # Replace with your Minio server address
        access_key="YOUR-ACCESSKEYID",  # Replace with your access key
        secret_key="YOUR-SECRETACCESSKEY",  # Replace with your secret key
        secure=True  # Set to False if not using HTTPS
    )

    # Parse the bucket name and object name from the minio_path
    bucket_name, object_name = os.path.split(minio_path)

    try:
        # Get the object from Minio
        response = client.get_object(bucket_name, object_name)
        # Read the data from the response
        data = response.read()
        # Parse the JSON data
        return json.loads(data)
    except S3Error as e:
        print("Error occurred.", e)
        raise

def upload_json_to_minio(data, minio_path):
    # Initialize the Minio client
    client = Minio(
        "play.min.io",  # Replace with your Minio server address
        access_key="YOUR-ACCESSKEYID",  # Replace with your access key
        secret_key="YOUR-SECRETACCESSKEY",  # Replace with your secret key
        secure=True  # Set to False if not using HTTPS
    )

    # Parse the bucket name and object name from the minio_path
    bucket_name, object_name = os.path.split(minio_path)

    try:
        # Convert data to JSON
        json_data = json.dumps(data).encode('utf-8')
        # Upload the JSON data to Minio
        client.put_object(bucket_name, object_name, data=io.BytesIO(json_data), length=len(json_data), content_type='application/json')
        print(f"Successfully uploaded results to {minio_path}")
    except S3Error as e:
        print("Error occurred while uploading.", e)
        raise

def main():
    parser = argparse.ArgumentParser(description="Run various chemical analysis functions.")
    # parser.add_argument('function', choices=['score_from_smi', 'get_chemoenzy_path_async', 'build_graph_from_async'],
    #                     help='Function to run')
    parser.add_argument('--minio_path', required=True, help='Minio path to the JSON file containing input parameters')
    args = parser.parse_args()

    # Download the JSON file from Minio
    params = download_json_from_minio(args.minio_path)

    results = {}
    if args.function == 'score_from_smi':
        if 'smiles' not in params:
            parser.error("The JSON file must contain 'smiles' for score_from_smi")
        results['sfscore'] = score_from_smi(params['smiles'])
    elif args.function == 'get_chemoenzy_path_async':
        if 'smiles' not in params:
            parser.error("The JSON file must contain 'smiles' for get_chemoenzy_path_async")
        explored_rxns, explored_nodes, start_node = get_chemoenzy_path_async(
            params['smiles'], params.get('max_depth', 10), params.get('chem_topk', 10), 
            params.get('max_num_templates', 250), params.get('max_branching', 15), 
            params.get('time_lim', 180)
        )
        results['explored_rxns'] = explored_rxns
        results['explored_nodes'] = explored_nodes
        results['start_node'] = start_node
        print("Completed get_chemoenzy_path_async()")
    elif args.function == 'build_graph_from_async':
        if 'explored_rxns' not in params or 'explored_nodes' not in params or 'start_node' not in params:
            parser.error("The JSON file must contain 'explored_rxns', 'explored_nodes', and 'start_node' for build_graph_from_async")
        results['graph'] = build_graph_from_async_wrapper(params['explored_rxns'], params['explored_nodes'], params['start_node'])
        print("Completed build_graph_from_async()")

    # Upload the results to Minio
    job_id: str = args.minio_path.split('/')[0].split('.')[0] # get ID from string like: `10b5006a-9851-43ef-b940-90305dade8c7/in/input.json`
    upload_json_to_minio(results, f"{job_id}/out/output.json")

if __name__ == "__main__":
    main()