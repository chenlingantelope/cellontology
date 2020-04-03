import pandas as pd
import obonet
import numpy as np
from cellontology import CellOntolgy

# load ontology
obo = 'data/cl.obo.txt'
f = open(obo, "r")
co = obonet.read_obo(f)
f.close()
celltype_dict = {}
for x in co.nodes:
    celltype_dict[co.nodes[x]['name']] = x

# define ontology
parent_choice = 'data/parent_choice.csv'
modification = 'data/modification.csv'
new_terms = 'data/new_terms.csv'

ontology = CellOntolgy(co, parent_choice, new_terms, modification)

# add new terms
ontology.add_new_terms()

# filter edges
ontology.filter_edges('is_a')

# celltype_ontology is the mapping between free annotation and ontology terms
celltype_ontology_dict = pd.read_csv('data/celltype_ontology.dict.txt', header=None)
all_ontology = np.unique(celltype_ontology_dict[1])

# modfiy ontology
ontology.modify_tree()
ontology.prune_to_tree(set(all_ontology))
ontology.graphviz_plot('ontology.pdf')

# processing for scanvi
ontology.modify_tree('data/modification_scanvi.csv')
ontology.graphviz_plot('ontology_scanvi.pdf')
ontology.flatten_tree()
ontology.graphviz_plot('ontology_scanvi_flat.pdf')

ontology.save_ontology('data/scanvi_tree.pkl')
# generate adjacency matrix from ontology
adjm = ontology.adjacency_matrix()
