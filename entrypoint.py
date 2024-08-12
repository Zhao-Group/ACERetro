import argparse
import sys
import joblib
import six
from sfscore import SFScore
from pathway_search_standalone.scripts.search_utils import hybridSearch, build_graph_from_async
import networkx as nx

def score_from_smi(smiles):
    PROJECT_DIRECTORY = '/app/ACERetro'
    sfscore_path = PROJECT_DIRECTORY + '/sfscore'
    sys.path.append(str(sfscore_path))

    sfscore_model = SFScore()
    sfscore_model.load()
    sfscore = sfscore_model.score_from_smi(smiles)
    print('SMILES:', smiles, ', SFScore:', sfscore)

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

    for gg in g.nodes():
        print(gg)
        print(g.nodes[gg])

def main():
    parser = argparse.ArgumentParser(description="Run various chemical analysis functions.")
    parser.add_argument('function', choices=['score_from_smi', 'get_chemoenzy_path_async', 'build_graph_from_async'],
                        help='Function to run')
    parser.add_argument('--smiles', help='SMILES string for molecule')
    parser.add_argument('--max_depth', type=int, default=10, help='Max depth for chemoenzy path search')
    parser.add_argument('--chem_topk', type=int, default=10, help='Chem topk for chemoenzy path search')
    parser.add_argument('--max_num_templates', type=int, default=250, help='Max num templates for chemoenzy path search')
    parser.add_argument('--max_branching', type=int, default=15, help='Max branching for chemoenzy path search')
    parser.add_argument('--time_lim', type=int, default=180, help='Time limit for chemoenzy path search')

    args = parser.parse_args()

    if args.function == 'score_from_smi':
        if not args.smiles:
            parser.error("--smiles is required for score_from_smi")
        score_from_smi(args.smiles)
    elif args.function == 'get_chemoenzy_path_async':
        if not args.smiles:
            parser.error("--smiles is required for get_chemoenzy_path_async")
        explored_rxns, explored_nodes, start_node = get_chemoenzy_path_async(
            args.smiles, args.max_depth, args.chem_topk, args.max_num_templates,
            args.max_branching, args.time_lim
        )
        print("Completed get_chemoenzy_path_async()")
        return explored_rxns, explored_nodes, start_node
    elif args.function == 'build_graph_from_async':
        if 'explored_rxns' not in globals() or 'explored_nodes' not in globals() or 'start_node' not in globals():
            parser.error("build_graph_from_async requires get_chemoenzy_path_async to be run first")
        build_graph_from_async_wrapper(explored_rxns, explored_nodes, start_node)
        print("Completed build_graph_from_async()")

if __name__ == "__main__":
    main()