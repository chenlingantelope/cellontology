from networkx import MultiDiGraph
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

class CellOntolgy(MultiDiGraph):

    def __init__(self, G, parent_choice_file=None, new_term_file=None, modification_file=None):
        super().__init__()
        self.__dict__ = G.__dict__.copy()
        self.all_nodes = set()
        self.labeled_nodes = set()
        self.parent_choice_file = parent_choice_file
        self.modification_file = modification_file
        self.new_term_file = new_term_file
        self.is_flat = False
        self.levels = None

    def add_new_terms(self):
        new_terms = pd.read_csv(self.new_term_file, delimiter='|')
        for i in new_terms.index:
            temp = new_terms.loc[i]
            if not temp['new_term_id'].startswith('#'):
                self.add_new_term(temp['new_term_id'], temp['new_term_name'],
                                      temp['new_term_def'], temp['parent_id'])

    def add_new_term(self, x, name, description, parent=None):
        self.add_node(x)
        self.nodes[x]['name'] = name
        self.nodes[x]['namespace'] = 'cell'
        self.nodes[x]['def'] = description
        if parent is not None:
            self.add_edge(x, parent, 'is_a')
        return None

    def filter_edges(self,edge_type):
        edges = list(self.edges)
        for x in edges:
            if x[2] != edge_type:
                self.remove_edge(x[0], x[1])
        return None

    def trace2root(self, node_id):
        self.all_nodes.add(node_id)
        parents = list(self.successors(node_id))
        parent_choice = pd.read_csv(self.parent_choice_file, delimiter='|')
        parent_choice.index = parent_choice['id']
        if len(parents) > 1:
            if node_id in parent_choice.index.values:
                parents = parent_choice.loc[node_id]['parents_id']
                parent = parents.split(';')[parent_choice.loc[node_id]['choice']]
            else:
                parents_names = [self.nodes[i]['name'] for i in parents]
                temp = [str(i) + '. ' + x for i, x in enumerate(parents_names)]
                choice = input(" for %s choose parent %s" % (self.nodes[node_id]['name'], " ".join(temp)))
                parent = parents[int(choice)]
                choice_record = [node_id, self.nodes[node_id]['name'],
                                 ';'.join(parents), ';'.join(parents_names), choice]
                f = open(self.parent_choice_file, 'a')
                f.write('|'.join(choice_record)+'\n')
                f.close()
        elif len(parents)==1:
            parent = parents[0]
        else:
            return None
        self.all_nodes.add(parent)
        self.trace2root(parent)
        return None

    def modify_tree(self, modification_file=None):
        if modification_file is None:
            modification_file = self.modification_file
        modification = pd.read_csv(modification_file, delimiter='|')
        for i in modification.index:
            if not modification.loc[i, 'node1'].startswith('#'):
                operation = modification.loc[i, 'operation']
                if operation == 'remove_node' and modification.loc[i, 'node1'] in self.nodes:
                    self.remove_node(modification.loc[i, 'node1'])
                elif operation == 'remove_edge':
                    self.remove_edge(modification.loc[i, 'node1'], modification.loc[i, 'node2'])
                elif operation == 'add_edge' and modification.loc[i, 'node1'] in self.nodes:
                    if modification.loc[i, 'node2'] not in self.nodes:
                        self.add_new_term(modification.loc[i, 'node2'], modification.loc[i, 'node2_name'],
                                          '', parent=None)
                    self.add_edge(modification.loc[i, 'node1'], modification.loc[i, 'node2'])
        return None

    def prune_to_tree(self, node_set=None):

        if len(self.labeled_nodes)==0 and node_set is None:
            raise ValueError('Have to provide a subset of nodes to prune the network to')
        elif node_set is not None and len(node_set)>0:
            self.labeled_nodes = node_set

        for x in self.labeled_nodes:
            self.trace2root(x)
        self.remove_nodes_from([node for node in self if node not in self.all_nodes])
        parent_choice = pd.read_csv(self.parent_choice_file, delimiter='|')
        parent_choice.index = parent_choice['id']
        for node_id in list(self.nodes):
            parents = list(self.successors(node_id))
            if len(parents) > 1:
                temp = parent_choice.loc[node_id]['parents_id']
                parent = temp.split(';')[parent_choice.loc[node_id]['choice']]
                for x in parents:
                    if x != parent:
                        self.remove_edge(node_id, x)
        self.remove_chains()
        return None

    def remove_chains(self):
        nodes = list(self.nodes)
        for x in nodes:
            if x not in self.labeled_nodes:
                parents = list(self.successors(x))
                children = list(self.predecessors(x))
                if len(children) == 1 and len(parents) == 1:
                    self.add_edge(children[0], parents[0])
                    self.remove_node(x)
        return None

    # find all leave cells
    def flatten_tree(self, levels=4):
        leaves = self.get_leaves()
        roots = self.get_roots()
        for root in roots:
            children = nx.ancestors(self, root)
            leave_children = [x for x in children if x in leaves]
            depth = [nx.shortest_path_length(self, x, root) for x in leave_children]
            last_root = root
            new_root = root
            assert np.max(depth) <= levels, "levels variable needs to be greater or equal " \
                                            "to the maximum number of existing levels"
            if np.max(depth) < levels:
                for i in range(levels - 1 - np.max(depth)):
                    new_root = root + "_" + str(i)
                    self.add_node(new_root)
                    self.nodes[new_root]['name'] = self.nodes[root]['name']
                    self.add_edge(last_root, new_root)
                    last_root = new_root
            for leaf in leave_children:
                dist = nx.shortest_path_length(self, leaf, new_root)
                last_node = leaf
                for i in range(levels - 1 - dist):
                    new_node = leaf + "_" + str(i)
                    parent = list(self.successors(last_node))[0]
                    self.add_node(new_node)
                    self.nodes[new_node]['name'] = self.nodes[leaf]['name']
                    # self.add_edge(new_node, last_node)
                    self.add_edge(new_node, parent)
                    self.add_edge(last_node, new_node)
                    self.remove_edge(last_node, parent)
                    last_node = new_node

        # if a leave node is different from its parent, and the parent only has one child,
        # add another child that has the same name as the parent
        for leaf in leaves:
            parent = list(self.successors(leaf))[0]
            grandparent = list(self.successors(parent))[0]
            siblings = list(self.predecessors(parent))
            if len(siblings) == 1 and parent in self.labeled_nodes:
                self.add_node(parent+"_1")
                self.nodes[parent+"_1"]['name'] = self.nodes[parent]['name']
                self.add_edge(parent+"_1", grandparent)
                self.add_edge(parent, parent + "_1")
                self.remove_edge(parent, grandparent)
                self.add_edge(leaf, parent+"_1")
                self.remove_edge(leaf, parent)

        self.is_flat = True
        self.levels = levels
        return None

    def graphviz_plot(self, filename=None):
        pos = nx.nx_agraph.graphviz_layout(self, prog='dot', root='cell')
        node_names = dict()
        node_colors = []
        for x in self.nodes:
            node_names[x] = self.nodes[x]['name']
            if x in self.labeled_nodes:
                node_colors.append('#eded11')
            else:
                node_colors.append('#11c1ed')

        plt.figure(figsize=(50, 10))
        nx.draw_networkx(self, pos, arrows=True, with_labels=False, node_color=node_colors)
        nx.draw_networkx_labels(self, pos, labels=node_names, font_size=5)
        if filename is not None:
            plt.savefig(filename)
        else:
            plt.show()

    def get_roots(self):
        roots = []
        for x in self.nodes():
            if len(list(self.successors(x))) == 0:
                roots.append(x)
        return roots

    def get_leaves(self):
        leaves = []
        for x in self.nodes():
            if len(list(self.predecessors(x))) == 0:
                leaves.append(x)
        return leaves

    def adjacency_matrix(self):
        assert self.is_flat is True,"The tree must be flattened before adjacency matrices can be generated. " \
                                    "See function flatten_tree"
        parents = self.get_roots()
        adjacency_matrix = []
        for i in range(self.levels-1):
            children = []
            children_list = []
            for i,x in enumerate(parents):
                y = list(self.predecessors(x))
                children.append(y)
                for z in y: children_list.append(z)

            matrix = pd.DataFrame(0, index=parents, columns=children_list)
            for i,p in enumerate(parents):
                matrix.loc[p,children[i]] = 1

            adjacency_matrix.append(matrix)
            parents = children_list

        return adjacency_matrix

    def save_ontology(self, filename):
        pkl.dump(self, open(filename, "wb"))
