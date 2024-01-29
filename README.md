# Synthetic field guided asynchronous chemoenzymatic synthesis planning 

This repository includes the scripts used for the manuscript "Synthetic field guided asynchronous chemoenzymatic synthesis planning."

**AceRetro** (**A**synchronous **C**hemo**E**nzymatic **Retro**synthesis) is a computer-aided chemoenzymatic synthesis planning tool. 

## Synthetic Field Score (SFScore)

### Usage
Dependencies
- numpy=1.21.5
- pytorch=1.12.1
- rdkit-pypi==2022.3.5
- map4==1.0 (MAP4 is not required if only ECFP4 embedding is used)

`./sfscore/sfscore.py` shows the basic usage. If you want to use SFScore elsewhere, follow the steps below.

```python
sfscore_path = [PROJECT_DIRECTORY]/'sfscore' # replace PROJECT_DIRECTORY
import sys
sys.path.append(str(sfscore_path))
from sfscore import SFScore
sfscore_model = SFScore()
# load the model in 'process_reaction_database/saved_model/ecfp4_4096_3_layer_epoch10.pt' by default
sfscore_model.load() 
smiles = 'O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO'
sfscore = sfscore_model.score_from_smi(smi)
print('SMILES:',smiles,', SFScore:',sfscore)
```

### Training the SFScore model
The scripts and checkpoints of trained models are in `./process_reaction_database`.

Create Conda environment though `conda env create -f ./process_reaction_database/environment.yml`

`process_reaction_database/fingerprint_embedding.ipynb`: Process USPTO and ECREACT database. Process ECFP4 and MAP4 fingerprint embedding.

`process_reaction_database/model_training.ipynb`: SFScore model training and results processing.

### Benchmarking SFScore

`evalueate_score/benchmark_zinc_one_step.ipynb`: SFScore results and one-step retrosynthesis(RXN4Chemistry)/retrobiosynthesis(ASKCOS-Enzy) on 11K molecules from ZINC database (Please refer to the next section for retrosynthesis/retrobiosynthesis).

`evalueate_score/process_zinc_one_step.ipynb`: Process one-step retrosynthesis results.

`evalueate_score/benchmark_askcos_results.ipynb`: SFScore results on pathways predicted by ASKCOS-hybrid (<https://github.com/itai-levin/hybmind>) in Nat Commun 13, 7747 (2022). (<https://doi.org/10.1038/s41467-022-35422-y>)

The processed molecule embedding data and ASKCOS benchmark data available in [Figshare]()

## **AceRetro** (***A***synchronous ***C***hemo***E***nzymatic ***Retro***synthesis)

### Install/Build
Create Conda environment though `conda env create -f ./pathway_search_standalone/environment.yml`. To meet the installation requirements of ASKCOS and RXN4Chemistry, the environment is different from the previous training environment.

Installation of ASKCOS ([Nat Commun, 13, 7747 (2022)](https://doi.org/10.1038/s41467-022-35422-y)):
- Download ASKCOS (<https://github.com/xuanliugit/askcos-core/tree/aceretro>) to `pathway_search_standalone/askcos-core`
- Download ASKCOS data (<https://github.com/xuanliugit/askcos-data>) to `pathway_search_standalone/askcos-core/askcos/data`

Installation of RXN4Chemistry ([Digital Discovery, 2023,2, 489-501](https://doi.org/10.1039%2Fd2dd00110a)): 
- Download RXN4Chemistry (<https://github.com/xuanliugit/rxn_cluster_token_prompt/tree/aceretro>) to `pathway_search_standalone/rxn_cluster_token_prompt`
- RXN4Chemistry's model is the intellectual property of IBM, which is not open-source but available in <https://rxn.res.ibm.com>. To use the open-source models download them in `pathway_search_standalone/rxn_cluster_token_prompt/models` following this [link](https://doi.org/10.6084/m9.figshare.20121944.v1).

### Usage
```python
# ./pathway_search_standalone/
from scripts.search_utils import hybridSearch
hybridSearch = hybridSearch()
smi = 'CC[C@@H](CO)NCCN[C@@H](CC)CO' # Target molecule's SMILES
explored_rxns, explored_nodes, start_node = hybridSearch.get_chemoenzy_path_async(smi, max_depth=10, chem_topk=10, max_num_templates=250, max_branching=15, time_lim=180)
```
To process the search result:
```python
# ./pathway_search_standalone/
from scripts.search_utils import build_graph_from_async
import networkx as nx
graph_json = build_graph_from_async(explored_rxns,explored_nodes,start_node)
g = nx.node_link_graph(graph_json) # reaction graph
# Check ./pathway_search_standalone/result_processing.ipynb to extract synthesis routes from reaction graph 
```

### Pathway search in the manuscript

`pathway_search_standalone/search_rxn4chemistry_askcos.ipynb`: Search hybrid pathways to 1k molecules from ZINC (Fig. 5) and case studies (Fig. 6,7,8).

`pathway_search_standalone/result_processing.ipynb`: Process the search results and translate graphs to readable pathways