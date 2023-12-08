import numpy as np
import pandas as pd
import urllib.parse
from IPython.display import Image,display


def compare_results_with_askcos(merged_path_length_pd,prioritizer,compare_all=False):
    tbtf_path_found = 0
    askcos_path_found = 0
    tbtf_found_askcos_not_found = 0
    tbtf_not_found_askcos_found = 0
    both_found = 0 
    tbtf_shorter = 0
    equal_length = 0
    tbtf_longer = 0

    mol_count = 0
    tbtf_shorter_mols = []
    tbtf_long_pathway_mols = []

    tbtf_shorter_1 = 0
    tbtf_shorter_2 = 0
    tbtf_shorter_3 = 0
    tbtf_shorter_more = 0
    for index in merged_path_length_pd.index:
        tbtf_path_length = merged_path_length_pd.loc[index, prioritizer]
        if compare_all:
            askcos_path_length = min(merged_path_length_pd.loc[index, 'bkms,reaxys'],merged_path_length_pd.loc[index, 'reaxys'],merged_path_length_pd.loc[index, 'bkms'])
        else:
            askcos_path_length = merged_path_length_pd.loc[index, 'bkms,reaxys']
        if tbtf_path_length < np.inf:
            tbtf_path_found += 1
            if askcos_path_length == np.inf:
                tbtf_found_askcos_not_found += 1
                #if tbtf_path_length >= 2:
                #    tbtf_long_pathway_mols.append(index)
        if askcos_path_length < np.inf:
            askcos_path_found += 1
            if tbtf_path_length == np.inf:
                tbtf_not_found_askcos_found += 1

        if tbtf_path_length < np.inf and askcos_path_length < np.inf:
            both_found += 1
            if tbtf_path_length < askcos_path_length:
                tbtf_shorter += 1
                tbtf_shorter_mols.append(index)
                # Check shoter 1/2/3/mpore
                if askcos_path_length - tbtf_path_length >3:
                    tbtf_shorter_more += 1
                elif askcos_path_length - tbtf_path_length ==3:
                    tbtf_shorter_3 += 1
                elif askcos_path_length - tbtf_path_length ==2:
                    tbtf_shorter_2 += 1
                elif askcos_path_length - tbtf_path_length ==1:
                    tbtf_shorter_1 += 1
                ##
            elif tbtf_path_length == askcos_path_length:
                equal_length += 1
            elif tbtf_path_length > askcos_path_length:
                tbtf_longer += 1
            #else:
            #    print(mol_count, tbtf_path_length, askcos_path_length)
        
        mol_count += 1

    print('Totol number of molecules:', mol_count,
          f'\n{prioritizer} found:', tbtf_path_found, 
          '\nASKCOS found:',askcos_path_found, 
          f'\n{prioritizer} found but ASKCOS not found:', tbtf_found_askcos_not_found, 
          f'\n{prioritizer} not found but ASKCOS found:',tbtf_not_found_askcos_found, 
          f'\nBoth {prioritizer} and ASKCOS found:',both_found, 
          f'\n{prioritizer} found shorter:',tbtf_shorter, 
          f'\n{prioritizer} found equiv:',equal_length, 
          f'\n{prioritizer} found longer:',tbtf_longer)
    print('\n1 step shorter:',tbtf_shorter_1,
          '\n2 steps shorter:',tbtf_shorter_2,
          '\n3 steps shorter:',tbtf_shorter_3,
          '\n>=4 steps shorter:',tbtf_shorter_more)
    return None



def print_rxns(spg_dict,prioritizer,smi):
    #images = []
    for rxn in spg_dict[smi][prioritizer].nodes[smi]['shortest_pathway']:
        url_rxn = urllib.parse.quote(rxn)
        img = Image(url=f'https://askcos.mit.edu/api/v2/draw/?smiles={url_rxn}&draw_map=false&highlight=false')
        #images.append(img)
        model = spg_dict[smi][prioritizer].nodes[rxn]['model']
        product = rxn.split('>>')[-1]
        product_sfscore = spg_dict[smi][prioritizer].nodes[product]['sfscore']
        print('Product SFSCore:',np.round(product_sfscore,3),model,': ',rxn)
        display(img)

def print_bypass_rxns(spg_dict,prioritizer,smi,smi_for_search):
    for rxn in spg_dict[smi][prioritizer].nodes[smi_for_search]['shortest_pathway']:
        url_rxn = urllib.parse.quote(rxn)
        img = Image(url=f'https://askcos.mit.edu/api/v2/draw/?smiles={url_rxn}&draw_map=false&highlight=false')
        #images.append(img)
        model = spg_dict[smi][prioritizer].nodes[rxn]['model']
        product = rxn.split('>>')[-1]
        product_sfscore = spg_dict[smi][prioritizer].nodes[product]['sfscore']
        print('Product SFSCore:',np.round(product_sfscore,3),model,': ',rxn)
        display(img)


def print_rxns_list(spg_dict,prioritizer,smi,rxns_list):
    for rxn in rxns_list:
        url_rxn = urllib.parse.quote(rxn)
        img = Image(url=f'https://askcos.mit.edu/api/v2/draw/?smiles={url_rxn}&draw_map=false&highlight=false')
        #images.append(img)
        model = spg_dict[smi][prioritizer].nodes[rxn]['model']
        product = rxn.split('>>')[-1]
        product_sfscore = spg_dict[smi][prioritizer].nodes[product]['sfscore']
        print('Product SFSCore:',np.round(product_sfscore,3),model,': ',rxn)
        display(img)