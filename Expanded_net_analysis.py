# -*- coding: utf-8 -*-
"""

@author: jwKim

using Dynaimcs_pyBoolnet` object from Model_read_using_pyboolnet.py module,
construct expanded network and analyze it.
each expanded node can be expressed in dictionary format or string format.
string format can be changed in `Expanded_node` class.
"""

import networkx as nx
import re

class Expanded_node:
    #expanded node format setting
    format_expanded_single_node = "{original_node}_{state_suffix}"#nodename_state
    format_composite_node_connector = "__and__"
    suffix_on = "1"
    suffix_off = "0"
    suffixes = [suffix_off, suffix_on]
    
    def __init__(self, dict_form):
        self.dict_form = dict_form
        self.str_form = self.get_str_form_from_dict_form(dict_form)

    def __repr__(self):
        return self.str_form
    
    def __eq__(self, dict_form):
        """expaned_node A_1 can be described as {"A":1}"""
        return self.dict_form == dict_form
    
    def is_composite(self):
        return len(self.dict_form) >=2
    
    def get_regulators_of_composite(self):
        """if this node is composite node, this returns list of regulators of this node.
        regulators are returned in str form."""
        if not self.is_composite():
            return []
        regulators_str_forms = []
        for node, state in self.dict_form.items():
            regulators_str_forms.append(self.get_str_form_from_dict_form({node:state}))
        return regulators_str_forms
    
    @classmethod
    def get_str_form_from_dict_form(cls, dict_form):
        """this method get the expanded node expressed as dictionary format,
        and convert it to string format, then return this string format expanded node."""
        dict_form_new = {}
        for node, state in dict_form.items():
            if int(state):
                suffix = cls.suffix_on
            else:
                suffix = cls.suffix_off
            dict_form_new[node] = suffix
            
        if len(dict_form) == 1:
            return cls.format_expanded_single_node.format(original_node=node, state_suffix=suffix)
        else:
            return cls._get_str_form_of_composite_node(dict_form_new)
    
    @classmethod
    def _get_str_form_of_composite_node(cls, dict_form_composite_node):
        """this method get the composite expanded node expressed as dictionary format,
        and convert it to string format.
        while converting it, node names are ordered.
        then return this string format composite expanded node."""
        nodes_in_composite = list(dict_form_composite_node.keys())
        nodes_in_composite.sort()

        expanded_nodes_in_composite = []
        for node in nodes_in_composite:
            suffix = dict_form_composite_node[node]
            expanded_node_form = cls.format_expanded_single_node.format(original_node=node, state_suffix=suffix)
            expanded_nodes_in_composite.append(expanded_node_form)

        return cls.format_composite_node_connector.join(expanded_nodes_in_composite)

    @classmethod
    def get_dict_form_from_str_form(cls, str_form):
        format_converted = cls.format_expanded_single_node.replace('{original_node}',"(.*)" )
        format_converted = format_converted.replace('{state_suffix}',"((?:{})|(?:{}))".format(cls.suffix_off, cls.suffix_on))
        dict_form = {}
        str_forms_splited = str_form.split(cls.format_composite_node_connector)
        for str_form_splited in str_forms_splited:
            match_object = re.match(format_converted, str_form_splited)
            nodename, state_suffix = match_object.groups()
            state_int = cls.suffixes.index(state_suffix)
            dict_form[nodename] = state_int
        
        return dict_form


def make_expanded_net_using_dynamics_pyboolnet(dynamics_pyboolnet, reduction=True, perturbation={}):
    expanded_net = Expanded_Network()
    primes = dynamics_pyboolnet.get_primes(reduction, perturbation)
    single_nodes, composite_nodes, edges = _make_expanded_net_from_prime_implicants(primes)

    for expanded_node_str_form, expanded_node in single_nodes.items():
        expanded_net.add_node(expanded_node_str_form)
        expanded_net.nodes[expanded_node_str_form]["info"] = expanded_node
    
    for expanded_node_str_form, expanded_node in composite_nodes.items():
        expanded_net.add_node(expanded_node_str_form)
        expanded_net.nodes[expanded_node_str_form]["info"] = expanded_node

    for edge in edges:
        expanded_net.add_edge(edge[0], edge[1])
    
    return expanded_net

def _make_expanded_net_from_prime_implicants(primes):
    """make networkx object whose nodes are expanded node
    using pyboolnet primes implicants format."""
    single_nodes = {}#str form: dict form
    composite_nodes ={} #str form: dict form
    edges = []#[(str form source, str form target),,,]

    for node, clauses in primes.items():
        for state in [0,1]:
            # first, calculate all single expanded nodes.
            node_dict_form = {node: state}
            expanded_single_node = Expanded_node(node_dict_form)
            node_str_form = str(expanded_single_node)
            single_nodes[node_str_form] = expanded_single_node

            # find regulators for each expanded single node and record their relationships.
            for composite_node_candidate in clauses[state]:
                if composite_node_candidate == {}:
                    # This case occurs when a specific node is fixed to 1 or 0.
                    # For example, prime["A"] = [[], [{}]].
                    # In this case, A is fixed to 1, so, change A_0 to A_1
                    opposite_state = 1-state
                    opposite_expanded_node_dict_form = {node: opposite_state}
                    opposite_expanded_node_str_form = Expanded_node.get_str_form_from_dict_form(opposite_expanded_node_dict_form)
                    edges.append((opposite_expanded_node_str_form, node_str_form))
                elif len(composite_node_candidate) >= 2:#composite node -> single node case
                    expanded_composite_node = Expanded_node(composite_node_candidate)
                    composite_node_str_form = str(expanded_composite_node)
                    if composite_node_str_form not in composite_nodes:
                        composite_nodes[composite_node_str_form] = expanded_composite_node
                    edges.append((composite_node_str_form, node_str_form))
                else:# the case when a regulator of single node is single node.
                    expanded_regulating_single_node = Expanded_node.get_str_form_from_dict_form(composite_node_candidate)
                    regulator_str_form = str(expanded_regulating_single_node)
                    edges.append((regulator_str_form, node_str_form))

    # find regulators for each composite node and record their relationships.
    for composite_node_str_form, expanded_composite_node in composite_nodes.items():
        for node, state in expanded_composite_node.dict_form.items():
            single_node_str_form = Expanded_node.get_str_form_from_dict_form({node: state})
            edges.append((single_node_str_form, composite_node_str_form))
    
    return single_nodes, composite_nodes, edges


class Expanded_Network(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_expanded_node_name(self, node_dict_form:dict):
        """this method get the expanded node expressed as dictionary format,
        and convert it to string format, 
        then return this string format composite expanded node."""
        return Expanded_node.get_str_form_from_dict_form(node_dict_form)
    
    @property
    def original_nodes(self):
        original_nodes_set = set()
        for node in self.nodes:
            node_info = self.nodes[node]["info"]
            if not node_info.is_composite():
                original_nodes_set.update(node_info.dict_form.keys())
        
        return original_nodes_set
    
    def get_node_info(self, expanded_node_name):
        """this method get the string format or dictionary format expression of expanded node,
        then returns the object of the expanded node."""
        if isinstance(expanded_node_name, dict):
            expanded_node_name = self.get_expanded_node_name(expanded_node_name)
        return self.nodes[expanded_node_name]["info"]
    
    def show_names_of_composite_nodes(self):
        """this method returns tuple of composite nodes (string format) in this expanded net."""
        composite_node_names = []
        for node_name in self.nodes:
            if self.nodes[node_name]["info"].is_composite():
                composite_node_names.append(node_name)
        
        return tuple(composite_node_names)
    
    def show_names_of_single_nodes(self):
        """this method returns tuple of single expanded nodes (string format) in this expanded net."""
        single_node_names = []
        for node_name in self.nodes:
            if not self.nodes[node_name]["info"].is_composite():
                single_node_names.append(node_name)
        
        return tuple(single_node_names)
    
    def get_regulator_names_of_node(self, node_name:"str or dict format"):
        """this method get the expanded node name, and return the regulators of the node"""
        if isinstance(node_name, dict):
            node_name = self.get_expanded_node_name(node_name)

        return tuple(self.predecessors(node_name))
    
    def get_subnetwork(self, expanded_node_names:"list of str or dict forms"):
        """Given a set of expanded nodes, return an expanded network object
        consisting only of those expanded nodes and the links connecting them.
        
        If a composite node activates at leats one of the given expanded single nodes,
        and all of its regulator single expanded nodes are included in the given set,
        then that composite node is also included in the subnetwork.

        For example, if A_1, B_1, and C_1 are given expanded nodes, and a link
        A_1__and__B_1 --> C_1 exists, then this composite node and link
        are also included in the subnetwork."""
        expanded_nodes_sub = self._convert_expanded_node_names_to_str_forms(expanded_node_names)
        
        composite_nodes_sub = []
        for composite_node in self.show_names_of_composite_nodes():
            if set(self.predecessors(composite_node)).issubset(expanded_nodes_sub):
                # if condition is satisfied when
                # regulator nodes of the composite node are all contained
                # in the expanded_nodes_sub,
                if set(self.successors(composite_node)).intersection(expanded_nodes_sub):
                    # if condition holds when
                    # there exists a target of the composite node in the expanded_nodes_sub,
                    composite_nodes_sub.append(composite_node)
        
        return nx.subgraph(self, expanded_nodes_sub+composite_nodes_sub)
    
    def _convert_expanded_node_names_to_str_forms(self, expanded_node_names:"list of str or dict forms"):
        """this method gets the list of expanded node names (string format or dictionary format).
        dictionary format node names are all converted to string format node names.
        and the list of expanded node names is returned."""
        expanded_node_names_str_forms = []
        for expanded_node in expanded_node_names:
            if isinstance(expanded_node, dict):
                expanded_node = self.get_expanded_node_name(expanded_node)
            expanded_node_names_str_forms.append(expanded_node)
        
        return expanded_node_names_str_forms
    
    def get_LDOI(self, expanded_node_names:"list of str or dict forms", only_return_single_nodes=True):
        """Given a set of expanded nodes, return the LDOI induced by them.
        
        If 'only_return_single_nodes' is True, return only the single expanded nodes
        that constitute the LDOI.
        
        It is assumed that the initially provided expanded nodes contain no contradictions."""
        expanded_node_names = set(self._convert_expanded_node_names_to_str_forms(expanded_node_names))
        single_expanded_nodes_in_LDOI = set()
        composite_expanded_nodes_in_LDOI = set()
        for expanded_node in expanded_node_names:
            if self.get_node_info(expanded_node).is_composite():
                composite_expanded_nodes_in_LDOI.add(expanded_node)
            else:
                single_expanded_nodes_in_LDOI.add(expanded_node)

        node_names_in_ldoi = set()
        for expanded_node_name in single_expanded_nodes_in_LDOI:
            node_names_in_ldoi.add(list(self.get_node_info(expanded_node_name).dict_form.keys())[0])
        
        while True:
            # If at least one new node is found, restart the while loop.
            # If no new nodes can be found, exit the while loop.
            already_checked = set()
            for expanded_node in expanded_node_names:
                break_outer_for = False
                for downstream_expanded_node in self.successors(expanded_node):
                    if downstream_expanded_node in expanded_node_names:
                        continue
                    if downstream_expanded_node in already_checked:
                        continue

                    expanded_node_obj = self.get_node_info(downstream_expanded_node)
                    if expanded_node_obj.is_composite():
                        if single_expanded_nodes_in_LDOI.issuperset(self.get_regulator_names_of_node(downstream_expanded_node)):
                            composite_expanded_nodes_in_LDOI.add(downstream_expanded_node)
                            break_outer_for = True
                            break
                        else:
                            already_checked.add(downstream_expanded_node)
                    else: # expanded_node_obj is single node
                        original_node_name = list(expanded_node_obj.dict_form.keys())[0]
                        if original_node_name in node_names_in_ldoi:
                            already_checked.add(downstream_expanded_node)
                            continue
                        single_expanded_nodes_in_LDOI.add(downstream_expanded_node)
                        node_names_in_ldoi.add(original_node_name)
                        break_outer_for = True
                        break

                if break_outer_for:
                    break # for
            
            expanded_nodes_names_new =single_expanded_nodes_in_LDOI.union(composite_expanded_nodes_in_LDOI) 
            if expanded_node_names == expanded_nodes_names_new:
                # no new expanded node is added
                break # while
            else:
                expanded_node_names = expanded_nodes_names_new
        
        if only_return_single_nodes:
            return single_expanded_nodes_in_LDOI
        else:
            return expanded_node_names