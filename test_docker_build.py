def score_from_smi():
    PROJECT_DIRECTORY = '/app/ACERetro'
    sfscore_path = PROJECT_DIRECTORY+'/sfscore' # replace PROJECT_DIRECTORY
    import sys
    sys.path.append(str(sfscore_path))

    from sfscore import SFScore
    sfscore_model = SFScore()

    # load the model in 'process_reaction_database/saved_model/ecfp4_4096_3_layer_epoch10.pt' by default
    sfscore_model.load() 
    smiles = 'O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO'
    sfscore = sfscore_model.score_from_smi(smiles)
    print('SMILES:',smiles,', SFScore:',sfscore)
    # SMILES: O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO , SFScore: [0.4322359 0.6146086]

def get_chemoenzy_path_async():
    # ./pathway_search_standalone/
    import sys
    import joblib
    import six
    sys.modules['sklearn.externals.joblib'] = joblib
    sys.modules['sklearn.externals.six'] = six
    from pathway_search_standalone.scripts.search_utils import hybridSearch

    hybridSearch = hybridSearch()
    smiles = 'CC[C@@H](CO)NCCN[C@@H](CC)CO' # Target molecule's SMILES
    # Conduct a 3 min asynchronous chemoenzymatic synthesis planning
    explored_rxns, explored_nodes, start_node = hybridSearch.get_chemoenzy_path_async(smiles, max_depth=10, chem_topk=10, max_num_templates=250, max_branching=15, time_lim=180)
    return explored_rxns,explored_nodes,start_node

def build_graph_from_async(explored_rxns,explored_nodes,start_node):
    # ./pathway_search_standalone/
    from scripts.search_utils import build_graph_from_async
    import networkx as nx
    graph_json = build_graph_from_async(explored_rxns,explored_nodes,start_node)
    g = nx.node_link_graph(graph_json) # reaction graph

    for gg in g.nodes():
        print(gg)
        print(g.nodes[gg])

if __name__ == "__main__":
    score_from_smi()
    print("[1/3] Completed score_from_smi()")
    explored_rxns, explored_nodes, start_node = get_chemoenzy_path_async()
    print("[2/3] Completed get_chemoenzy_path_async()")
    build_graph_from_async(explored_rxns, explored_nodes, start_node)
    print("[3/3] build_graph_from_async()")