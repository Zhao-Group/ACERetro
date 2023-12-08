import torch
import torch.nn as nn
import numpy as np
from map4 import MAP4Calculator
from rdkit import Chem
from rdkit.Chem import AllChem
import os

class NeuralNet(nn.Module):
    def __init__(self, lrate, loss_fn, in_size, out_size):
        super(NeuralNet, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(in_size, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, out_size),
            nn.Sigmoid()
        )
        self.loss_fn = loss_fn
        self.optimizer = torch.optim.Adam(self.model.parameters(), lrate)
    def forward(self, x):
        y = self.model(x)
        return y
    def step(self, x, y):
        loss = self.loss_fn(x, y)
        return loss


class NeuralNet_1(nn.Module):
    def __init__(self, lrate, loss_fn, in_size, out_size):
        super(NeuralNet_1, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(in_size, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, out_size),
            nn.Sigmoid()
        )
        self.loss_fn = loss_fn
        self.optimizer = torch.optim.Adam(self.model.parameters(), lrate)
    def forward(self, x):
        y = self.model(x)
        return y
    def step(self, x, y):
        loss = self.loss_fn(x, y)
        return loss

class NeuralNet_3(nn.Module):
    def __init__(self, lrate, loss_fn, in_size, out_size):
        super(NeuralNet_3, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(in_size, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, out_size),
            nn.Sigmoid()
        )
        self.loss_fn = loss_fn
        self.optimizer = torch.optim.Adam(self.model.parameters(), lrate)
    def forward(self, x):
        y = self.model(x)
        return y
    def step(self, x, y):
        loss = self.loss_fn(x, y)
        return loss

class NeuralNet_5(nn.Module):
    def __init__(self, lrate, loss_fn, in_size, out_size):
        super(NeuralNet_5, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(in_size, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, 300),
            nn.ReLU(),
            nn.Linear(300, out_size),
            nn.Sigmoid()
        )
        self.loss_fn = loss_fn
        self.optimizer = torch.optim.Adam(self.model.parameters(), lrate)
    def forward(self, x):
        y = self.model(x)
        return y
    def step(self, x, y):
        loss = self.loss_fn(x, y)
        return loss

project_root = os.path.dirname(os.path.dirname(__file__))

class SFScore():
    def __init__(self, path = None, fp_type = None, dim = None, layer_num = None):
        if path == None:
            self.path = os.path.join(project_root, 'process_reaction_database', 'saved_model', 'ecfp4_4096_3_layer_epoch10.pt')
        else:
            self.path = path
        
        if fp_type == None:
            self.fp_type = 'ECFP4'
        else:
            self.fp_type = fp_type

        if dim == None:
            self.dim = 4096
        else:
            self.dim = dim

        if layer_num == None:
            self.layer_num = 3
        else:
            self.layer_num = layer_num

        self.MAP4 = MAP4Calculator(dimensions=self.dim, is_folded=True)
        self.device = torch.device('cpu')

    def load(self):
        print(f'Loading model {self.path}')
        if self.layer_num == 1:
            self.sfscore_model = NeuralNet_5(lrate=1e-5, loss_fn=None, in_size=self.dim, out_size=2)
        if self.layer_num == 3:
            self.sfscore_model = NeuralNet_3(lrate=1e-5, loss_fn=None, in_size=self.dim, out_size=2)
        if self.layer_num == 5:
            self.sfscore_model = NeuralNet_5(lrate=1e-5, loss_fn=None, in_size=self.dim, out_size=2)
        self.sfscore_model.load_state_dict(torch.load(self.path, map_location=torch.device('cpu')))
        self.sfscore_model.eval()
        return self

    def get_fp(self, smi):
        mol = Chem.MolFromSmiles(smi)
        if self.fp_type == 'MAP4':
            try:
                map4_fp = self.MAP4.calculate(mol)
                map4_fp = torch.tensor(map4_fp, dtype=torch.float, device=self.device)
            except:
                map4_fp = None
                print(f"Cannot get map4 fp from {smi}, its sfscore is set as None")
            return map4_fp
        elif self.fp_type == 'ECFP4':
            ecfp4_fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=self.dim))
            ecfp4_fp = torch.tensor(ecfp4_fp, dtype=torch.float, device=self.device)
            return ecfp4_fp
    
    def score_from_smi(self, smi):
        fp = self.get_fp(smi)
        if fp == None:
            return None
        with torch.no_grad():
            out = self.sfscore_model(fp)
            out = out.to('cpu').numpy()
        return out
    
    def get_fp_many(self, smi_list):
        mol_list = [Chem.MolFromSmiles(smi) for smi in smi_list]
        if self.fp_type == 'MAP4':
            try:
                map4_fp = self.MAP4.calculate_many(mol_list)
                map4_fp = torch.tensor(map4_fp, dtype=torch.float, device=self.device)
            except:
                map4_fp = None
                print(f"Cannot get map4 fp from {smi}, its sfscore is set as None")
            return map4_fp
        elif self.fp_type == 'ECFP4':
            ecfp4_fp = np.array([AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=self.dim) for mol in mol_list])
            ecfp4_fp = torch.tensor(ecfp4_fp, dtype=torch.float, device=self.device)
            return ecfp4_fp
        
    def score_from_smi_many(self, smi_list):
        fps = self.get_fp_many(smi_list)
        with torch.no_grad():
            out = self.sfscore_model(fps)
            out = out.to('cpu').numpy()
            #out = out.tolist()
        return out
if __name__ == '__main__':
    sfscore_model = SFScore()
    sfscore_model.load()
    smi = 'O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO'
    sfscore = sfscore_model.score_from_smi(smi)
    print('SMILES:',smi,', SFScore:',sfscore.tolist())

    smi_list = ['O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO','CCCCCCO']
    sfscore_list = sfscore_model.score_from_smi_many(smi_list)
    print('SMILES list:',smi_list,', SFScore list',sfscore_list.tolist())

