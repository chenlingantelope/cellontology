cl_ontology_file = 'cl-basic.obo'

import obonet
f = open(cl_ontology_file, "r")
co = obonet.read_obo(f)
f.close()

celltype_dict = {}
for x in co.nodes:
    celltype_dict[co.nodes[x]["name"].lower()] = x
