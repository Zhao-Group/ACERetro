[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10578664.svg)](https://doi.org/10.5281/zenodo.10578664)

<div align='center'>
<picture>
  <source srcset="assets/aceretro_logo.png" width='50%'>
  <img alt="ACERetro logo" src="/assets/" width="50%" >
</picture>
</div>


# Chemoenzymatic synthesis planning guided by synthetic potential score 

This repository includes the scripts used for the manuscript **"Chemoenzymatic synthesis planning guided by synthetic potential score (formerly Synthetic Field Score, SFScore)."** 

**AceRetro** (**A**synchronous **C**hemo**E**nzymatic **Retro**synthesis) is a computer-aided chemoenzymatic synthesis planning tool. 

## Synthetic Potential Score, SPScore (formerly Synthetic Field Score, SFScore)

### Usage
Dependencies
- conda=22.9.0
- python=3.9.12
- numpy=1.21.5
- pytorch=1.12.1
- rdkit-pypi==2022.3.5
- map4==1.0 (MAP4 is not required if only ECFP4 embedding is used)

Operating systems
- Ubuntu 20.04.5 LTS

`./sfscore/sfscore.py` shows the basic usage. If you want to use SFScore elsewhere, follow the steps below.

Create Conda environment though `conda env create -f ./process_reaction_database/environment.yml`

Activate Conda environment `conda activate sfscore_train`

```python
sfscore_path = [PROJECT_DIRECTORY]/'sfscore' # replace PROJECT_DIRECTORY
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
```

### Training the SFScore model
The scripts and checkpoints of trained models are in `./process_reaction_database`.

`process_reaction_database/fingerprint_embedding.ipynb`: Process USPTO and ECREACT database. Process ECFP4 and MAP4 fingerprint embedding. (~ 2h)

`process_reaction_database/model_training.ipynb`: SFScore model training and results processing. (~ 40 min)

### Benchmarking SFScore

`evalueate_score/benchmark_zinc_one_step.ipynb`: SFScore results and one-step retrosynthesis(RXN4Chemistry)/retrobiosynthesis(ASKCOS-Enzy) on 11K molecules from ZINC database (Cost ~ 43.5 hours) (Please refer to the next section for one-step retrosynthesis/retrobiosynthesis prediction).

`evalueate_score/process_zinc_one_step.ipynb`: Process one-step retrosynthesis results. (Fig. 2)

`evalueate_score/benchmark_askcos_results.ipynb`: SFScore results on pathways predicted by ASKCOS-hybrid (<https://github.com/itai-levin/hybmind>) in Nat Commun 13, 7747 (2022). (<https://doi.org/10.1038/s41467-022-35422-y>) (Fig. 3 and ASKCOS pathways in Fig. 5)

## **AceRetro** (***A***synchronous ***C***hemo***E***nzymatic ***Retro***synthesis)

### Install/Build
Installation of ASKCOS ([Nat Commun, 13, 7747 (2022)](https://doi.org/10.1038/s41467-022-35422-y)):
- Download ASKCOS (<https://github.com/xuanliugit/askcos-core/tree/aceretro>) to `./pathway_search_standalone/askcos-core`
- Download ASKCOS data (<https://github.com/xuanliugit/askcos-data>) to `./pathway_search_standalone/askcos-core/askcos/data`

```
# Download ASKCOS
## Directory: ./pathway_search_standalone
git clone https://github.com/xuanliugit/askcos-core.git
## Change branch to aceretro
cd askcos-core
git switch aceretro

# Download ASKCOS data
## Directory: ./pathway_search_standalone/askcos-core/askcos/
git lfs clone https://github.com/xuanliugit/askcos-data.git
## Rename 'askcos-data' to 'data'
mv askcos-data data
```

Installation of RXN4Chemistry ([Digital Discovery, 2023,2, 489-501](https://doi.org/10.1039%2Fd2dd00110a)): 
- Download RXN4Chemistry (<https://github.com/xuanliugit/rxn_cluster_token_prompt/tree/aceretro>) to `pathway_search_standalone/rxn_cluster_token_prompt`
- RXN4Chemistry's model is the intellectual property of IBM, which is not open-source but available in <https://rxn.res.ibm.com>. To use the open-source models download them in `pathway_search_standalone/rxn_cluster_token_prompt/models` following this [link](https://doi.org/10.6084/m9.figshare.20121944.v1).

```
# Download RXN4Chemistry
## Directory: ./pathway_search_standalone
git clone https://github.com/xuanliugit/rxn_cluster_token_prompt.git
## Change branch to aceretro
cd rxn_cluster_token_prompt
git switch aceretro

## Model Directory: ./pathway_search_standalone/rxn_cluster_token_prompt/models
```

```
# Create environmentt
conda create -n aceretro python=3.6.13
conda activate aceretro
pip install rdkit-pypi

# Install packages for ASKCOS
pip install -r ./pathway_search_standalone/askcos-core/requirements.txt

# Install packages for RXN4Chemistry
cd pathway_search_standalone/rxn_cluster_token_prompt
pip install -e .

# Install other packages
pip install pubchempy==1.0.4
pip install git+https://github.com/reymond-group/map4@v1.0
conda install -c tmap tmap
```

### Usage
```python
# ./pathway_search_standalone/
from scripts.search_utils import hybridSearch
hybridSearch = hybridSearch()
smiles = 'CC[C@@H](CO)NCCN[C@@H](CC)CO' # Target molecule's SMILES
# Conduct a 3 min asynchronous chemoenzymatic synthesis planning
explored_rxns, explored_nodes, start_node = hybridSearch.get_chemoenzy_path_async(smiles, max_depth=10, chem_topk=10, max_num_templates=250, max_branching=15, time_lim=180)
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

`pathway_search_standalone/search_rxn4chemistry_askcos.ipynb`: Search hybrid pathways to 1k molecules from ZINC (Fig. 5) (around 6.25 days) and case studies (Fig. 6,7,8).

`pathway_search_standalone/result_processing.ipynb`: Process the search results and translate graphs to readable pathways (less than 5 min)
