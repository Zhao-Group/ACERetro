import numpy as np
import pandas as pd
import pickle

from tqdm import tqdm
from time import sleep

from rdkit import Chem
import matplotlib.pyplot as plt

import sys
from pathlib import Path
root_directory = Path('__file__').parent.parent.resolve()
print('root_directory: ', root_directory)
project_directory = root_directory.parent
print('project_directory: ', project_directory)

sfscore_path = project_directory/'sfscore'
sys.path.append(str(sfscore_path))
print(sys.path)

from sfscore import SFScore
sfscore_model = SFScore()
sfscore_model.load()



def standardize_smiles_list (smi_ls):
    # Ref: https://github.com/itai-levin/hybmind/blob/main/analyze_templates/Notebooks/compare_model_outputs.ipynb
    std_ls = []
    for smi in tqdm(smi_ls):
        mol = Chem.MolFromSmiles(smi)
        for atm in mol.GetAtoms():
            atm.SetAtomMapNum(0)
        std_ls.append(Chem.MolToSmiles(mol))
    return std_ls

def calculate_sfscore_from_path_dict(path_dict):
    
    with open("../process_reaction_database/data/bio_products_uni.pkl","rb") as f:
        bio_products_uni = pickle.load(f)
    with open("../process_reaction_database/data/chem_products_uni.pkl","rb") as f:
        chem_products_uni = pickle.load(f)
    sfscore_subset = standardize_smiles_list(bio_products_uni + chem_products_uni)
    
    for name, path in path_dict.items():
        print(f"Processing {name}...")
        np.random.seed(42)
        df = pd.read_csv(path)
        if len(df) > 12000:
            rand_inds = np.random.choice(len(df), 12000)
            df = df.loc[rand_inds, :]
        try:
            df_subset = standardize_smiles_list(df['smiles'])
        except KeyError:
            df_subset = standardize_smiles_list(df['SMILES'])
        
        df_in_train = set(sfscore_subset).intersection(df_subset)
        df_not_in_train_idxs = [x not in df_in_train for x in df_subset]
        df_sfscore = sfscore_model.score_from_smi_many(df_subset)
        np.savetxt(f"../data/zinc/{name}_sfscore.csv", df_sfscore[df_not_in_train_idxs], delimiter=",")
        print(f"Saved {name} sfscore to ../data/zinc/{name}_sfscore.csv")
        
if __name__ == "__main__":
    path_dict = {
        'fda': '../data/zinc/fda.csv',
        'biogenic': '../data/zinc/biogenic_2015.csv',
        'MOSES': '../data/zinc/MOSES_chems.txt',
    }
    calculate_sfscore_from_path_dict(path_dict)
