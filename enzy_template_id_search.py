import numpy as np
import pandas as pd
import re

from rxnmapper import RXNMapper # pip install rxnmapper
from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs

def add_EC_Number_to_json(json_file: dict):
    """
    Input: python dictionary object

    """
    total_count = 0
    total_successes = 0
    errors = 0
    for index, data in json_file['explored_rxns'].items():
        # data is a list of dicts. each dict is a molecule... or similar.
        for el in data:
            total_count += 1
            try:
                if el.get('is_major_precursor') and el.get('template_id') and el.get('smiles'):
                    # print("Building template thing for smile: ", el['is_major_precursor'][0])
                    res = enzy_template_id_search(el['template_id'], el['smiles'])
                    json_file['explored_rxns'][index][0]['EC_Number'] = res.EC_Number
                    total_successes += 1
            except Exception as e:
                errors += 1
                print("ERROR add_EC_Number_to_json. Error:", e, "smiles:", el['smiles'], "precursor:", el['is_major_precursor'][0])
    # print("Total errors: ", errors)
    # print("Total successes: ", total_successes)
    # print("Total count: ", total_count)
    return json_file

def v2_add_EC_Number_to_json(json_file: dict): 
    """
    Input: python dictionary object

    """
    total_count = 0
    total_successes = 0
    errors = 0
    for index, data in json_file['graph'].items():
        total_count += 1
        try: 
            # print(data)
            if data.get('major_precursor') and data.get('template_id') and data.get('smiles'):
                print("Building template thing for smile: ", data['major_precursor'][0])
                res = enzy_template_id_search(data['template_id'], data['smiles'])
                print("Adding EC: ", res.EC_Number)
                json_file['graph'][index]['EC_Number'] = res.EC_Number
                total_successes += 1
        except Exception as e: 
            errors += 1
            print("ERROR add_EC_Number_to_json. Error:", e, "smiles:", data['smiles'], "precursor:", data['is_major_precursor'][0])
    # print("Total errors: ", errors)
    # print("Total successes: ", total_successes)
    # print("Total count: ", total_count)
    return json_file


def enzy_template_id_search(template_id, target_mol):
    """
    Searches for templates matching a given template ID and finds the most similar molecule to a target molecule.

    Args:
    template_id (str/int): The template ID to search for in the template database.
    target_mol (str): The SMILES representation of the target molecule.

    Returns:
    DataFrame Row: The row from the template DataFrame corresponding to the most similar molecule.
    """

    templates = pd.read_json('bkms-retro.templates.bkms.json.gz')
    template_df = pd.read_json('bkms-templates.df.json.gz')
    index = int(templates[templates['_id']==template_id]['index'])

    ref_mol = Chem.MolFromSmiles(target_mol)
    ref_fp = Chem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=1024)

    new_template_df = template_df.loc[template_df['index']==index,['products']]
    mol_smarts_list = [new_template_df.loc[idx,'products'] for idx in new_template_df.index]
    mol_list = [Chem.MolFromSmiles(mol_smarts) for mol_smarts in mol_smarts_list]

    fp_list = [Chem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mol_list]
    sim_list = [DataStructs.TanimotoSimilarity(ref_fp,fp) for fp in fp_list]
    largest_smi_index = np.argmax(sim_list)
    return template_df.loc[new_template_df.index[largest_smi_index]]


# def add_major_precursor_to_json(json_file):
#     """
#     Here is the code to process the original JSON output. For all reaction nodes, they will have a new entry: 
#     'major_precursor': ['COc1cc2nc(N3CCNCC3)nc(N)c2cc1OC', 'O=C(Cl)c1ccco1']

#     For major precursor chemical nodes, they will have a new entry: 
#     'is_major_precursor': True. For other chemical nodes, they will NOT have this entry.
    
#     Editors note, it doesn't have `'is_major_precursor': True,` it just has the key if it exists. Doesn't have the key if not.
#     """
#     rxn_mapper = RXNMapper()
#     explored_rxns = json_file['explored_rxns']
#     for idx, rxn_list in explored_rxns.items():
#         for rxn_idx, rxn in enumerate(rxn_list):
#             rxn_smiles = rxn['id']
#             major_precursor = []
#             mapped_rxn_smiles = rxn_mapper.get_attention_guided_atom_maps([rxn_smiles])[0]['mapped_rxn']
#             reactants = mapped_rxn_smiles.split('>>')[0].split('.')
#             for reactant_idx, reactant in enumerate(reactants):
#                 if_mapped = re.search(r"\[.+:\d+\]", reactant)
#                 if if_mapped:
#                     reactant_smiles = rxn['smiles_split'][reactant_idx]
#                     major_precursor.append(reactant_smiles)
#                     for node_idx, node in enumerate(json_file['explored_nodes'][idx]):
#                         if node['id'] == reactant_smiles:
#                             node['major_precursor'] = True
#                             json_file['explored_nodes'][idx][node_idx] = node
#             json_file['explored_rxns'][idx][rxn_idx]['is_major_precursor'] = major_precursor
#     return json_file

def v2_add_major_precursor_to_json(json_file):
    """
    Process the JSON output to add major precursor information to reaction nodes in the graph.
    For reaction nodes, add a new entry: 'major_precursor': ['SMILES1', 'SMILES2', ...]
    For major precursor chemical nodes, add a new entry: 'is_major_precursor': True
    """
    print("Starting add_major_precursor_to_json function")
    rxn_mapper = RXNMapper()
    total_items = len(json_file['graph'])
    total_major_precursors = 0
    print(f"Number of items in json_file['graph']: {total_items}")
    
    for rxn_smiles, rxn_data in json_file['graph'].items():
        if rxn_data['type'] == 'reaction':
            print("rxn_smiles", rxn_smiles)
            major_precursor = []
            mapped_rxn_smiles = rxn_mapper.get_attention_guided_atom_maps([rxn_smiles])[0]['mapped_rxn']
            reactants = mapped_rxn_smiles.split('>>')[0].split('.')
            print("reactants", reactants)
            
            for reactant_idx, reactant in enumerate(reactants):
                if_mapped = re.search(r"\[.+:\d+\]", reactant)
                print("if_mapped", if_mapped)
                if if_mapped:
                    major_precursor.append(rxn_data['smiles_split'][reactant_idx])
                    total_major_precursors += 1
                    
                    # Update the corresponding chemical node
                    for smiles, node in json_file['graph'].items():
                        if node['type'] == 'chemical' and smiles == reactant:
                            node['is_major_precursor'] = True
            
            rxn_data['major_precursor'] = major_precursor
    
    print("Finished processing all nodes")
    print(f"Total items in graph: {total_items}")
    print(f"Total major precursors added: {total_major_precursors}")
    return json_file


# path = './results/fda_drug_all/CHEMBL2.json' # an ACERetro output json file
# path = 'output_with_svg_on_type_chemical.json' # an ACERetro output json file
# with open(path, 'r') as f:
#     json_file = json.load(f)
# json_file_revised = add_major_precursor_to_json(json_file)
# json_file_revised
