
from rdkit import Chem
from rdkit.Chem import AllChem

import numpy as np
import pubchempy as pcp
import time
import heapq
import json

import sys
from pathlib import Path
from pathway_search_standalone.scripts.directory_utils import project_directory

print('project_directory: ', project_directory())
askcos_core_path = project_directory()/'pathway_search_standalone'/'askcos-core'
rxn_cluster_path = project_directory()/'pathway_search_standalone'/'rxn_cluster_token_prompt'
biocatalysis_path = project_directory()/'pathway_search_standalone'/'biocatalysis-model'
sfscore_path = project_directory()/'sfscore'
sys.path.append(str(askcos_core_path))
sys.path.append(str(rxn_cluster_path))
sys.path.append(str(sfscore_path))
rxn_cluster_parent_path = project_directory()/'pathway_search_standalone'
sys.path.append(str(rxn_cluster_parent_path))

from sfscore import SFScore
from rxn_cluster_token_prompt.model import RXNClusterTokenPrompt
import askcos.global_config as gc
import askcos.retrosynthetic.transformer as retro_trans
from askcos.utilities.buyable.pricer_sa_em import Pricer


def get_canonical_smiles(smi):
        mol = Chem.MolFromSmiles(smi)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        return smiles

def score_precursor(precursor, ppg=None, rank = 1):
    """Score a given precursor using a combination of the template relevance score and a heuristic rule
    Source: https://github.com/ASKCOS/askcos-core/blob/main/askcos/prioritization/precursors/relevanceheuristic.py
    Args:
        precursor (dict): dictionary of precursor to score

    Returns:
        float: combined relevance heuristic score of precursor
    """
    scores = []
    #necessary_reagent_atoms = precursor['necessary_reagent'].count('[')/2.
    #for smiles in precursor['smiles_split']:
    # If buyable, basically free
    if ppg:
        scores.append(- ppg / 100.0)
    # Else, use heuristic
    mol = Chem.MolFromSmiles(precursor)
    total_atoms = mol.GetNumHeavyAtoms()
    ring_bonds = sum([b.IsInRing() - b.GetIsAromatic()
                        for b in mol.GetBonds()])
    chiral_centers = len(Chem.FindMolChiralCenters(mol))

    scores.append(
        - 2.00 * np.power(total_atoms, 1.5)
        - 1.00 * np.power(ring_bonds, 1.5)
        - 2.00 * np.power(chiral_centers, 2.0)
    )
    sco = np.sum(scores) 
    #return -sco * np.power(rank, 0.5)
    return -sco * (2-np.exp((1-rank)/10))


class hybridSearch:
    def __init__(self,pricer=None,sfef=0.75,deepef=0.75,geo_iso=True):
        '''
        Args:
            sfef (int): synthesis field exploration factor
            deepef (int): depth exploration factor
        '''
        if pricer == None:
            self.pricer = Pricer(use_db=False)
            self.pricer.load()
        else:
            self.pricer = pricer
        self.askcos_enzy = retro_trans.RetroTransformer(use_db=False, template_set='bkms', template_prioritizer='bkms')
        self.askcos_enzy.load('bkms')

        self.rxn4chem_retro_model = RXNClusterTokenPrompt(n_best=1)

        self.sfscore_model = SFScore()
        self.sfscore_model.load()
        self.available_list = ['[Li]O','O[Na]','O=C(O[K])O[K]','O=[N+]([O-])O[K]',
                                'O[K]', '[Na]I', 'O=C(O)O[Na]','Cl[Ca]Cl','N#C[Na]',
                                '[O-][Cl+]O[Na]','Cl[Sn]Cl','C#C']
        self.sfef = sfef
        self.deepef = deepef
        self.geo_iso = geo_iso
    def unit_test(self):
        print('Start unit test.')
        ppg = self.pricer.lookup_smiles('[CH2-]CCC',source = ['EM','SA'],alreadyCanonical=True)
        if ppg == 99:
            print('Checking pricer...Done')
        else:
            print('Error on pricer.')
        rxn4chem_result = self.rxn4chem_retro_model.retro_predict(['Oc1c(Br)ccc2[nH]c3c(c12)CCNC3'], reorder_by_backward_likelihood=True, verbose=False, fap=0.6)
        if rxn4chem_result:
            print('Checking RXN4Chemistry retro model...Done')
        else:
            print('Error on RXN4Chemistry retro model.')
        askcos_result = self.askcos_enzy.get_outcomes('Fc1ccc(-c2ccc3[nH]c4c(c3c2)CCNC4)cn1',fast_filter_threshold=0, max_num_templates=250, 
                                                max_cum_prob=0.9999, cluster_precursors=False, use_ban_list=False,
                                                )
        if askcos_result:
            print('Checking ASKCOS(Enzy) model...Done')
        else:
            print('Error on ASKCOS(Enzy) model.')
        print('Done unit test.')
        return None

    def get_rxn4chem_result(self,smiles,fap=0.6):
        predict_dict = self.rxn4chem_retro_model.retro_predict([smiles], reorder_by_backward_likelihood=True, verbose=False, fap=fap)
        result_list = []
        for product, pred_results in predict_dict.items():
            unique_predictions = []
            for (prediction,confidence,round_trip_prediction, round_trip_confidence,rclass) in pred_results:
                if prediction not in unique_predictions:
                    single_rxn = {}
                    # single_rxn['smiles'] = product
                    single_rxn['smiles'] = smiles
                    single_rxn['bw_confidence'] = confidence
                    single_rxn['confidence'] = round_trip_confidence
                    single_rxn['rclass'] = rclass
                    single_rxn['children'] = []
                    for precursor in prediction.split('.'):
                        single_rxn['children'].append({'smiles':precursor})
                    unique_predictions.append(prediction)
                    result_list.append(single_rxn)
        return result_list
    
    def rebuild_rxn4chem_dict(self,smiles, rxn4chem_result, depth, ban_smiles,sfef=None,geo_iso=None):
        available_list = self.available_list
        rxn_list=[]
        node_list_dict = {}
        if not sfef:
            sfef = self.sfef
        if not geo_iso:
            geo_iso = self.geo_iso
        # rank = 1 
        for rank, rxn in enumerate(rxn4chem_result):
            has_ban_smiles = False
            precursor_smi_list = []

            for child in rxn['children']:
                try:
                    child_smi_break = child['smiles'].replace('~Cl','').replace('~','.')
                    child_smi_break = [get_canonical_smiles(child_smi) for child_smi in child_smi_break.split('.')]
                except:
                    child_smi_break = [get_canonical_smiles(child_smi) for child_smi in child['smiles'].split('.')]
                for child_smi in child_smi_break:
                #if '.' not in child['smiles']:
                    child_smi_std = child_smi
                    # get_canonical_smiles(child_smi)
                    # Convert '[nH+]' to 'n'
                    if '[nH+]' in child_smi_std or 'C(=O)[O-]' in child_smi_std or 'O=C([O-])' in child_smi_std:
                        child_smi_std = child_smi_std.replace('[nH+]','n').replace('C(=O)[O-]','C(=O)O').replace('O=C([O-])','O=C(O)')
                        child_smi_std = get_canonical_smiles(child_smi_std)
                    
                    precursor_smi_list.append(child_smi_std)
                    if child_smi_std == smiles or 'O+' in child_smi_std:
                        has_ban_smiles = True
                        break
                    if ban_smiles and child_smi_std in ban_smiles:
                        has_ban_smiles = True
                        break
                    #and rank < chem_topk
            if not has_ban_smiles:
                rxn_dict = dict()
                all_precursor_smi = '.'.join(precursor_smi_list)
                target_smiles = rxn['smiles']
                rxn_dict['id'] = all_precursor_smi + '>>' + target_smiles
                rxn_dict['type'] = 'reaction'
                rxn_dict['rclass'] = rxn['rclass']
                rxn_dict['bw_confidence'] = rxn['bw_confidence']
                rxn_dict['confidence'] = rxn['confidence']
                rxn_dict['model'] = ['RXN4Chem']
                rxn_dict['rank'] = rank + 1
                rxn_dict['depth'] = depth
                rxn_dict['precursor_smiles'] = all_precursor_smi
                rxn_dict['smiles_split'] = precursor_smi_list
                rxn_list.append(rxn_dict)

                score_list = []
                non_terminal_smi_list=[]
                for precursor_smi in rxn_dict['smiles_split']:
                    if precursor_smi not in node_list_dict:
                        node_dict = dict()
                        node_dict['id'] = precursor_smi
                        node_dict['depth'] = depth
                        node_dict['type'] = 'chemical'
                        mol = Chem.MolFromSmiles(precursor_smi)
                        total_atoms = mol.GetNumHeavyAtoms()
                        if precursor_smi in available_list or total_atoms<=3:
                            ppg = 1
                            node_dict['purchase_price'] = ppg
                            node_dict['terminal'] = True
                        else:
                            if 'Li' in precursor_smi or 'Na' in precursor_smi or 'Mg' in precursor_smi or 'Ca' in precursor_smi or 'K' in precursor_smi:
                                try:
                                    smi_pcp = pcp.get_properties(['CanonicalSMILES'],precursor_smi,'smiles')[0]['CanonicalSMILES']
                                    ppg_list = []
                                    
                                    for smi in smi_pcp.split('.'):
                                        ppg_sub = self.pricer.lookup_smiles(smi,source=['EM','SA'],alreadyCanonical=False)
                                        ppg_list.append(ppg_sub)
                                    ppg = min(ppg_list)
                                #    smi_modify = precursor_smi.replace('/','').replace('\\','')
                                #    ppg = self.pricer.lookup_smiles(smi_modify,source = ['EM','SA'],alreadyCanonical=False)
                                except:
                                    ppg = self.pricer.lookup_smiles(precursor_smi,source=['EM','SA'],alreadyCanonical=True)
                            else:
                                ppg = self.pricer.lookup_smiles(precursor_smi,source=['EM','SA'],alreadyCanonical=True)
                            if not geo_iso:
                                precursor_smi_nongeo = precursor_smi.replace('/','').replace('\\','')
                                ppg_nongeo = self.pricer.lookup_smiles(precursor_smi_nongeo,source=['EM','SA'],alreadyCanonical=False)
                                if (ppg_nongeo>0) and (ppg_nongeo<=100):
                                    ppg = ppg_nongeo
                            node_dict['purchase_price'] = ppg
                            node_dict['terminal'] = (ppg > 0) and (ppg <= 100)
                        #else:
                        #if precursor_smi in available_list:
                        #    node_dict['terminal'] = True
                            
                        score = score_precursor(precursor_smi, ppg, rxn_dict['rank'])
                        node_dict['score'] = score
                        if not node_dict['terminal']:
                            score_list.append(score)
                            sfscore = self.sfscore_model.score_from_smi(precursor_smi).tolist()
                            node_dict['sfscore'] = sfscore
                            #node_dict['chem_score'] = score*(1-sfscore[0])
                            #node_dict['enzy_score'] = score*(1-sfscore[1])
                            non_terminal_smi_list.append(precursor_smi)
                        #if '.' in rxn_dict['precursor_smiles']:
                        #    node_dict['score'] = 2 * score_precursor(precursor_smi, ppg, rxn_dict['rank'])
                        #else:
                        #    node_dict['score'] = score_precursor(precursor_smi, ppg, rxn_dict['rank'])
                        node_list_dict[precursor_smi]=node_dict
                if non_terminal_smi_list:
                    for non_terminal_smi in non_terminal_smi_list:
                        node_list_dict[non_terminal_smi]['score'] = sum(score_list)/len(score_list)
                        # node_list_dict[non_terminal_smi]['chem_score'] = node_list_dict[non_terminal_smi]['score']*(1-node_list_dict[non_terminal_smi]['sfscore'][0])
                        # node_list_dict[non_terminal_smi]['enzy_score'] = node_list_dict[non_terminal_smi]['score']*(1-node_list_dict[non_terminal_smi]['sfscore'][1])
                        node_list_dict[non_terminal_smi]['chem_score'] = node_list_dict[non_terminal_smi]['score']*(1-sfef*node_list_dict[non_terminal_smi]['sfscore'][0])
                        node_list_dict[non_terminal_smi]['enzy_score'] = node_list_dict[non_terminal_smi]['score']*(1-sfef*node_list_dict[non_terminal_smi]['sfscore'][1])
                        
                        # node_list_dict[non_terminal_smi]['score'] = node_list_dict[non_terminal_smi]['score']*(2-math.exp((1-(len(score_list))/5)))
        node_list = list(node_list_dict.values())
        return rxn_list, node_list
        


    def get_rxn4chem_rxn(self,smiles,
                        depth = None, fap=0.6, ai_model='12class-tokens-2021-05-14',
                        chem_topk = 10, ban_smiles=None):
        if depth == None:
            depth = 0
        try:
            rxn4chem_result = self.get_rxn4chem_result(smiles,fap=fap)
        except Exception as e:
            print("Error on RXN4Chemistry!", e)
            return [],[]
        # print('SUCCESS')
        rxn4chem_result = rxn4chem_result[:chem_topk]
        rxn_list, node_list = self.rebuild_rxn4chem_dict(smiles, rxn4chem_result, depth, 
                                                     ban_smiles=ban_smiles,
                                                    )
        return rxn_list, node_list


    def rebuild_enzy_rxn_dict(self,target_smiles, rxns, depth, ban_smiles=None, sfef=None,geo_iso=None):
        available_list = self.available_list
        rxn_list=[]
        node_list_dict = {}
        if not sfef:
            sfef = self.sfef
        if not geo_iso:
            geo_iso = self.geo_iso

        for rxn in rxns:
            has_ban_smiles = False
            for precursor_smi in rxn['smiles_split']:
                if ban_smiles and precursor_smi in ban_smiles:
                    has_ban_smiles = True
            if not has_ban_smiles:
                rxn['id'] = rxn['smiles'] + '>>' + target_smiles
                rxn['type'] = 'reaction'
                # rxn['rtype'] = 'Enzy'
                rxn['model'] = ['bkms']
                rxn['depth'] = depth
                rxn['precursor_smiles'] = rxn['smiles']
                rxn_list.append(rxn)

                score_list = []
                non_terminal_smi_list=[]

                for precursor_smi in rxn['smiles_split']:
                    if precursor_smi not in node_list_dict:
                        node_dict = dict()
                        node_dict['id'] = precursor_smi
                        node_dict['depth'] = depth
                        # ppg = self.pricer.lookup_smiles(precursor_smi,source=['EM','SA'],alreadyCanonical=True)
                        #if not geo_iso:
                        #    precursor_smi = precursor_smi.replace('/','').replace('\\','')
                        
                        ppg = self.pricer.lookup_smiles(precursor_smi,source=['EM','SA'],alreadyCanonical=True)
                        if not geo_iso:
                            precursor_smi_nongeo = precursor_smi.replace('/','').replace('\\','')
                            ppg_nongeo = self.pricer.lookup_smiles(precursor_smi_nongeo,source=['EM','SA'],alreadyCanonical=False)
                            if (ppg_nongeo>0) and (ppg_nongeo<=100):
                                ppg = ppg_nongeo
                        node_dict['purchase_price'] = ppg
                        node_dict['type'] = 'chemical'
                        if precursor_smi in available_list:
                            node_dict['terminal'] = True
                        else:
                            node_dict['terminal'] = (ppg > 0) and (ppg <= 100)

                        score = score_precursor(precursor_smi, ppg, rxn['rank'])
                        node_dict['score'] = score
                        if not node_dict['terminal']:
                            score_list.append(score)
                            sfscore = self.sfscore_model.score_from_smi(precursor_smi).tolist()
                            node_dict['sfscore'] = sfscore
                            #node_dict['chem_score'] = score*(1-sfscore[0])
                            #node_dict['enzy_score'] = score*(1-sfscore[1])
                            non_terminal_smi_list.append(precursor_smi)

                        #if '.' in rxn['smiles']:
                        #    node_dict['score'] = 2 * score_precursor(precursor_smi, ppg, rxn['rank'])
                        #else:
                        #    node_dict['score'] = score_precursor(precursor_smi, ppg, rxn['rank'])
                        
                        node_list_dict[precursor_smi]=node_dict
                if non_terminal_smi_list:
                    for non_terminal_smi in non_terminal_smi_list:
                        node_list_dict[non_terminal_smi]['score'] = sum(score_list)/len(score_list)
                        # old chem_enzy_score
                        #node_list_dict[non_terminal_smi]['chem_score'] = node_list_dict[non_terminal_smi]['score']*(1-node_list_dict[non_terminal_smi]['sfscore'][0])
                        #node_list_dict[non_terminal_smi]['enzy_score'] = node_list_dict[non_terminal_smi]['score']*(1-node_list_dict[non_terminal_smi]['sfscore'][1])
                        # new score np.exp(-0.5*)
                        node_list_dict[non_terminal_smi]['chem_score'] = node_list_dict[non_terminal_smi]['score']*(1-sfef*node_list_dict[non_terminal_smi]['sfscore'][0])
                        node_list_dict[non_terminal_smi]['enzy_score'] = node_list_dict[non_terminal_smi]['score']*(1-sfef*node_list_dict[non_terminal_smi]['sfscore'][1])
                        # node_list_dict[non_terminal_smi]['score'] = node_list_dict[non_terminal_smi]['score']*(2-math.exp((1-(len(score_list))/5)))
        node_list = list(node_list_dict.values())
        return rxn_list, node_list

    def get_enzy_rxn(self,smiles, max_num_templates = 1000, depth = None, 
                    max_branching=25, ban_smiles=None):
        if depth == None:
            depth = 0
        precursors = self.askcos_enzy.get_outcomes(smiles,fast_filter_threshold=0, max_num_templates=max_num_templates, 
                                                max_cum_prob=0.9999, cluster_precursors=False, use_ban_list=False
                                                )
        # print('t_bio',precursors)
        precursors = precursors[:max_branching]
        rxn_list, node_list = self.rebuild_enzy_rxn_dict(smiles, precursors, depth, 
                                                    ban_smiles=ban_smiles,)
        return rxn_list, node_list


    def get_onestep_paths(self,smiles,sfscore,depth,ban_smiles=None, 
                        chem_topk = 10, fully_hybrid = False, 
                        max_num_templates = 1000, max_branching = 25,
                        margin = 0.2):
        # final_result = []

        if sfscore[0]-sfscore[1]>margin and not fully_hybrid:
            search_rtype = "chem"
            print(f"Working on {smiles}, sfscore = {sfscore}, search_rtype = {search_rtype}")
            rxn_list, node_list = self.get_rxn4chem_rxn(smiles,depth=depth, chem_topk=chem_topk,ban_smiles=ban_smiles)
            print(f'Got {len(rxn_list)} Chem RXNs')

        elif sfscore[0]-sfscore[1]<-margin and not fully_hybrid:
            search_rtype = "enzy"
            print(f"Working on {smiles}, sfscore = {sfscore}, search_rtype = {search_rtype}")
            rxn_list, node_list = self.get_enzy_rxn(smiles,max_num_templates = max_num_templates, depth = depth, max_branching=max_branching,ban_smiles=ban_smiles)
            print(f'Got {len(rxn_list)} Enzy RXNs')

        else:
            search_rtype = "hybrid"
            print(f"Working on {smiles}, sfscore = {sfscore}, search_rtype = {search_rtype}")
            chem_rxn_list, chem_node_list = self.get_rxn4chem_rxn(smiles, depth=depth, chem_topk=chem_topk,ban_smiles=ban_smiles)
            print(f'Got {len(chem_rxn_list)} Chem RXNs')
            enzy_rxn_list, enzy_node_list = self.get_enzy_rxn(smiles,max_num_templates = max_num_templates, depth = depth, max_branching=max_branching,ban_smiles=ban_smiles)
            print(f'Got {len(enzy_rxn_list)} Enzy RXNs')
            #add fs
            rxn_list = chem_rxn_list + enzy_rxn_list
            all_node_dict = {}
            all_node_list = chem_node_list + enzy_node_list
            # print('total rxn and node',len(rxn_list), len(all_node_list))
            # print('all_node_list',all_node_list)
            for node_dict in all_node_list:
                if node_dict['id'] not in all_node_dict:
                    all_node_dict[node_dict['id']] = node_dict
                elif node_dict['score'] < all_node_dict[node_dict['id']]['score']:
                    all_node_dict[node_dict['id']] = node_dict
            node_list = list(all_node_dict.values())
        return rxn_list, node_list

    def get_chemoenzy_path_sync(self,smiles,
                        max_depth=10, chem_topk=10, 
                        max_num_templates=250, max_branching=15,
                        fully_hybrid=False, time_lim=180, margin=0.15, deepef=None):
        if not deepef:
            deepef = self.deepef
        timeout = time.time() + time_lim
        frontier = []
        explored_rxns = {}
        explored_nodes = {}
        node_num = 0
        smiles = get_canonical_smiles(smiles)
        start_node = {}
        sfscore = self.sfscore_model.score_from_smi(smiles).tolist()
        mol_score = score_precursor(smiles)
        # print(mol_score)
        start_node['id'] = smiles
        start_node['type'] = 'chemical'
        start_node['score'] = mol_score
        start_node['sfscore'] = sfscore
        start_node['depth'] = -1
        start_node['terminal'] = False
        heapq.heappush(frontier,(mol_score,node_num,start_node))
        '''
        sfscore = self.sfscore_model.score_from_smi(smiles).tolist()
        rxn_list, node_list = self.get_onestep_paths(smiles,sfscore,
                                                depth=0, ban_smiles=[smiles], 
                                                chem_topk=chem_topk, fully_hybrid = fully_hybrid, 
                                                max_num_templates=max_num_templates, 
                                                max_branching=max_branching)

        explored_rxns[smiles] = rxn_list
        explored_nodes[smiles] = node_list
        for node in node_list:
            if not node['terminal']:
                heapq.heappush(frontier,(node['score'],node_num,node))
                node_num += 1
        '''
        while len(frontier) != 0 and time.time() < timeout:
            score, num, frontire_node = heapq.heappop(frontier)
            #precursor = 
            # for precursor in precursors:
            depth = frontire_node['depth'] + 1
            if not frontire_node['terminal'] and frontire_node['id'] not in explored_rxns and depth < max_depth:
                try:
                    sfscore = frontire_node['sfscore']
                except:
                    sfscore = self.sfscore_model.score_from_smi(frontire_node['id']).tolist()
                rxn_list, node_list = self.get_onestep_paths(frontire_node['id'],sfscore,depth,
                                                        ban_smiles=[smiles], 
                                                        chem_topk=chem_topk,fully_hybrid = fully_hybrid, 
                                                        max_num_templates=max_num_templates, 
                                                        max_branching=max_branching, margin=margin)
                explored_rxns[frontire_node['id']] = rxn_list
                explored_nodes[frontire_node['id']] = node_list
                if rxn_list and node_list:
                    for node in node_list:
                        if not node['terminal'] and node['id'] not in explored_rxns:
                            score = node['score'] * np.power(depth, deepef)
                            heapq.heappush(frontier,(score,node_num,node))
                            node_num += 1
        return explored_rxns, explored_nodes, start_node
    
    def get_chemoenzy_path_async(self,smiles,
                        max_depth=10, chem_topk = 10, 
                        max_num_templates=250, max_branching=15,
                        time_lim=180, deepef=None, sfef=None):
        timeout = time.time() + time_lim
        if not deepef:
            deepef = self.deepef
        if not sfef:
            sfef = self.sfef
        frontier = []
        explored_status = {}
        explored_rxns = {}
        explored_nodes = {}
        node_num = 0
        smiles = get_canonical_smiles(smiles)
        start_node = {}
        sfscore = self.sfscore_model.score_from_smi(smiles).tolist()
        mol_score = score_precursor(smiles)
        # print(mol_score)
        start_node['id'] = smiles
        start_node['type'] = 'chemical'
        start_node['score'] = mol_score
        start_node['sfscore'] = sfscore
        start_node['chem_score'] = mol_score*(1-sfscore[0])
        start_node['enzy_score'] = mol_score*(1-sfscore[1])
        start_node['depth'] = -1
        start_node['terminal'] = False
        explored_status[smiles] = {'chem': False, 'enzy': False}
        heapq.heappush(frontier,(start_node['chem_score'],node_num,'chem',start_node))
        heapq.heappush(frontier,(start_node['enzy_score'],node_num,'enzy',start_node))
        while len(frontier) != 0 and time.time() < timeout:
            weighted_score, num, search_type, frontire_node = heapq.heappop(frontier)
            depth = frontire_node['depth'] + 1
            if not frontire_node['terminal'] and depth < max_depth and not explored_status[frontire_node['id']][search_type]:
                if search_type == 'chem':
                    # TODO: remove the SFScore in the print
                    print(f"Working on {frontire_node['id']}, sfscore = {frontire_node['sfscore']}, weighted_score = {weighted_score}, depth={depth}, search_rtype = chem")
                    rxn_list, node_list = self.get_rxn4chem_rxn(frontire_node['id'],depth=depth, chem_topk=chem_topk,ban_smiles=[smiles])
                    print(f'Got {len(rxn_list)} Chem RXNs')
                    
                elif search_type == 'enzy':
                    print(f"Working on {frontire_node['id']}, sfscore = {frontire_node['sfscore']}, weighted_score = {weighted_score},  depth={depth}, search_rtype = enzy")
                    rxn_list, node_list = self.get_enzy_rxn(frontire_node['id'],max_num_templates = max_num_templates, depth = depth, max_branching=max_branching,ban_smiles=[smiles])
                    print(f'Got {len(rxn_list)} Enzy RXNs')
                explored_status[frontire_node['id']][search_type] = True
                
                if rxn_list and node_list:
                    if explored_rxns.get(frontire_node['id']):
                        explored_rxns[frontire_node['id']].extend(rxn_list)
                        explored_nodes[frontire_node['id']].extend(node_list)
                    else:
                        explored_rxns[frontire_node['id']] = rxn_list
                        explored_nodes[frontire_node['id']] = node_list
                    for node in node_list:
                        if not node['terminal'] and node['id'] not in explored_status:
                            heapq.heappush(frontier,(node['chem_score']*np.power(depth+1, deepef),node_num,'chem',node))
                            node_num += 1
                            heapq.heappush(frontier,(node['enzy_score']*np.power(depth+1, deepef),node_num,'enzy',node))
                            explored_status[node['id']] = {'chem': False, 'enzy': False}
                            node_num += 1
        return explored_rxns, explored_nodes, start_node

def build_graph_from_async(explored_rxns, explored_nodes, start_node):
    '''
    return:
        dict for nx.node_link_graph
    '''
    data = {'directed': True, 'multigraph': False, 'graph': {}, 'nodes': [], 'links': []}
    data['nodes'].append(start_node)
    all_rxns = list(explored_rxns.values())
    all_nodes = list(explored_nodes.values())
    all_rxns_flatten = [rxn for sublist in all_rxns for rxn in sublist]
    all_nodes_flatten = [node for sublist in all_nodes for node in sublist if node['id']!=start_node['id']]
    data['nodes'].extend(all_rxns_flatten)
    data['nodes'].extend(all_nodes_flatten)
    for rxn in all_rxns_flatten:
        product = rxn['id'].split('>>')[-1]
        product_link = {'target':product, 'source':rxn['id']}
        data['links'].append(product_link)
        reactants = rxn['id'].split('>>')[0]
        if '.' in reactants:
            if start_node['id'] not in reactants.split('.'):
                for reactant in reactants.split('.'):
                    reactant_link = {'target':rxn['id'], 'source':reactant}
                    data['links'].append(reactant_link)
        elif reactants != start_node['id']:
            reactant_link = {'target':rxn['id'], 'source':reactants}
            data['links'].append(reactant_link)
    return data


    

