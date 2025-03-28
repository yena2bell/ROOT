# -*- coding: utf-8 -*-
"""

@author: jwKim

Model_read의 Dynamics_pyBoolnet으로부터 expanded net을 만들고, 분석하는 코드
expanded node format setting 에서 만들어지는 expanded net의 format 조정 가능.
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
        """expaned_node A_1 == {"A":1}"""
        return self.dict_form == dict_form
    
    def is_composite(self):
        return len(self.dict_form) >=2
    
    def get_regulators_of_composite(self):
        """이 node가 composite node일 경우, 그 composite node가 만족하기 위한 regulators의
        str forms 를 return한다."""
        if not self.is_composite():
            return []
        regulators_str_forms = []
        for node, state in self.dict_form.items():
            regulators_str_forms.append(self.get_str_form_from_dict_form({node:state}))
        return regulators_str_forms
    
    @classmethod
    def get_str_form_from_dict_form(cls, dict_form):
        """dict form으로 주어진 expanded node를 str form의 expanded node로 바꿔서
        return한다."""
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
        """composite node 를 str form으로 만들 때, nodes 의 순서를
        고려해서, 같은 dict form 에서 다른 str form이 생기는 것을
        방지한다."""
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
    """나중에 integer form logic으로부터 만드는 함수도 추가하려 함."""
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
    """pyboolnet의 primes를 사용하여 expanded network networkx 객체를 생성"""
    single_nodes = {}#str form: dict form
    composite_nodes ={} #str form: dict form
    edges = []#[(str form source, str form target),,,]

    for node, clauses in primes.items():
        for state in [0,1]:
            #일단 single expanded nodes 먼저 다 구하고,
            node_dict_form = {node: state}
            expanded_single_node = Expanded_node(node_dict_form)
            node_str_form = str(expanded_single_node)
            single_nodes[node_str_form] = expanded_single_node

            #각 expanded single node에 대한 regulators 를 찾아주고, 그 관계를 기록함.
            for composite_node_candidate in clauses[state]:
                if composite_node_candidate == {}:
                    #특정 node를 1이나 0으로 fix한 후 구하는 경우 발생.
                    #prime["A"] = [[], [{}]] 와 같이.
                    #이 경우 A가 1로 fix 되는 것이므로
                    #A_0 --> A_1 로 해준다.
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
                else:#single node -> single node case
                    expanded_regulating_single_node = Expanded_node.get_str_form_from_dict_form(composite_node_candidate)
                    regulator_str_form = str(expanded_regulating_single_node)
                    edges.append((regulator_str_form, node_str_form))
    
    #각 composite node 의 regulator nodes를 찾아 기록한다.
    for composite_node_str_form, expanded_composite_node in composite_nodes.items():
        for node, state in expanded_composite_node.dict_form.items():
            single_node_str_form = Expanded_node.get_str_form_from_dict_form({node: state})
            edges.append((single_node_str_form, composite_node_str_form))
    
    return single_nodes, composite_nodes, edges


class Expanded_Network(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_expanded_node_name(self, node_dict_form:dict):
        """dict 형태로 넣은 expanded node (composite node 포함)
        에 대해 그 expanded node가 이 객체에서 사용하고 있는
        str 형태의 node name을 return 해준다."""
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
        """str form, 또는 dict form으로 주어진 expanded_node_name에 대해,
        그 node의 Expanded_node 객체를 return한다."""
        if isinstance(expanded_node_name, dict):
            expanded_node_name = self.get_expanded_node_name(expanded_node_name)
        return self.nodes[expanded_node_name]["info"]
    
    def show_names_of_composite_nodes(self):
        """composite node 의 이름을 tuple에 모아서 return해준다."""
        composite_node_names = []
        for node_name in self.nodes:
            if self.nodes[node_name]["info"].is_composite():
                composite_node_names.append(node_name)
        
        return tuple(composite_node_names)
    
    def show_names_of_single_nodes(self):
        """single node 의 이름을 tuple 에 모아서 return해준다."""
        single_node_names = []
        for node_name in self.nodes:
            if not self.nodes[node_name]["info"].is_composite():
                single_node_names.append(node_name)
        
        return tuple(single_node_names)
    
    def get_regulator_names_of_node(self, node_name:"str 또는 dict form"):
        """주어진 expanded node의 이름으로부터, 그것의 regulators names 을 찾아서
        tuple 로 return해준다."""
        if isinstance(node_name, dict):
            node_name = self.get_expanded_node_name(node_name)
            #dict form 이면 적절한 composite form으로 바꿔준다.
        return tuple(self.predecessors(node_name))
    
    def get_subnetwork(self, expanded_node_names:"list of str or dict forms"):
        """어떤 expanded nodes가 주어지면, 그 expanded nodes와 그 expanded nodes를 잇는 links
        만으로 이루어진 expanded network 객체를 return한다.
        
        어떠한 composite node가 넣어준 expanded single nodes를 키고, 
        그 regulators인 single expanded nodes가 전부 넣어준 expanded single nodes에 포함되면
        그 composite node로 subnetwork에 포함시켜준다.
        즉 A_1, B_1, C_1이 넣어준 값이고, A_1__and__B_1 --> C_1인 links가 존재하면
        이것도 sub network에 포함된다는 것."""
        expanded_nodes_sub = self._convert_expanded_node_names_to_str_forms(expanded_node_names)
        
        composite_nodes_sub = []
        for composite_node in self.show_names_of_composite_nodes():
            if set(self.predecessors(composite_node)).issubset(expanded_nodes_sub):
                #composite node의 regulators가 전부 expanded_nodes_sub에 들어가면서,
                if set(self.successors(composite_node)).intersection(expanded_nodes_sub):
                    #그 composite node의 target 중 expanded_nodes_sub인 것이 있을 때,
                    composite_nodes_sub.append(composite_node)
        
        return nx.subgraph(self, expanded_nodes_sub+composite_nodes_sub)
    
    def _convert_expanded_node_names_to_str_forms(self, expanded_node_names:"list of str or dict forms"):
        """list에 expanded nodes가 str form이나 dict form이 섞여서 담겨 있을 때,
        이것들을 전부 stro form으로 통일해서 return 한다."""
        expanded_node_names_str_forms = []
        for expanded_node in expanded_node_names:
            if isinstance(expanded_node, dict):
                expanded_node = self.get_expanded_node_name(expanded_node)
            expanded_node_names_str_forms.append(expanded_node)
        
        return expanded_node_names_str_forms
    
    def get_LDOI(self, expanded_node_names:"list of str or dict forms", only_return_single_nodes=True):
        """어떤 expanded nodes가 주어지면, 그로 인해 만들어지는 LDOI 를 return한다.
        'only_return_single_nodes'==True이면 LDOI 를 이루는 single expanded nodes만 return
        
        처음에 넣어주는 expanded nodes에 모순이 없다고 가정하고 진행한다."""
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
            #하나라도 새로운 node를 찾으면 다시 while 문을 시작한다.
            #새로운 node가 더는 찾아지지 않으면 while 문을 탈출한다.
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
                    else:#single node
                        original_node_name = list(expanded_node_obj.dict_form.keys())[0]
                        if original_node_name in node_names_in_ldoi:
                            already_checked.add(downstream_expanded_node)
                            continue
                        single_expanded_nodes_in_LDOI.add(downstream_expanded_node)
                        node_names_in_ldoi.add(original_node_name)
                        break_outer_for = True
                        break

                if break_outer_for:
                    break#outer for
            
            expanded_nodes_names_new =single_expanded_nodes_in_LDOI.union(composite_expanded_nodes_in_LDOI) 
            if expanded_node_names == expanded_nodes_names_new:
                #새로운 node가 추가되지 않음.
                break#while
            else:
                expanded_node_names = expanded_nodes_names_new
        
        if only_return_single_nodes:
            return single_expanded_nodes_in_LDOI
        else:
            return expanded_node_names



        

    
    # def get_downstream_single_nodes(self, nodes_perturbed:"list 안에 dict form이나 str form의 nodes"):
    #     """주어진 expanded nodes에 대해, expanded network 정보로부터 downstream single nodes를
    #     구한다. expanded_nodes는 expanded network에서 사용하는 node 이름이어야 한다.
    #     downstream single nodes는 그 nodes로부터 composite node 를 거치지 않고 
    #     연결될 수 있는 single expanded nodes를 의미한다."""
    #     nodes_perturbed_str_forms = self._convert_dict_forms_to_str_forms(nodes_perturbed)
    #     #nodes_peturbed에 dict form이 섞여있으면 다 str form으로 바꿔준다.
        
    #     if Make_Using_Prime_Implicants.expanded_nodes_group_contains_contradiction(nodes_perturbed_str_forms):
    #         print("perturbed nodes are contradict.")
        
    #     single_nodes = []
    #     composite_nodes = []
    #     for node in nodes_perturbed_str_forms:
    #         if self.nodes[node]["composite"]:
    #             composite_nodes.append(node)
    #         else:
    #             single_nodes.append(node)
        
    #     num_nodes_selected = len(single_nodes) + len(composite_nodes)
    #     num_new = 0
    #     while num_new < num_nodes_selected:#전 loop에 비해 새로 찾은 것이 있는가?
    #         num_new = num_nodes_selected
    #         for node in single_nodes+composite_nodes:
    #             for successor_node in self.successors(node):
    #                 if not Make_Using_Prime_Implicants.two_expanded_nodes_groups_are_contradict([successor_node], single_nodes+composite_nodes):
    #                     if not self.nodes[successor_node]["composite"]:
    #                         single_nodes.append(successor_node)
            
    #         single_nodes = list(set(single_nodes))
    #         num_nodes_selected = len(single_nodes) + len(composite_nodes)
        
    #     return single_nodes
        
    
    # def get_LDOI(self, nodes_perturbed:"list 안에 dict form이나 str form의 nodes"):
    #     """perturbed 된 nodes 를 넣어줄 때, 그것으로 인한 LDOI 를 계산해서 return한다."""
    #     nodes_perturbed_str_forms = self._convert_dict_forms_to_str_forms(nodes_perturbed)
    #     #nodes_peturbed에 dict form이 섞여있으면 다 str form으로 바꿔준다.
        
    #     if Make_Using_Prime_Implicants.expanded_nodes_group_contains_contradiction(nodes_perturbed_str_forms):
    #         print("perturbed nodes are contradict.")
        
    #     single_nodes = []
    #     composite_nodes = []
    #     for node in nodes_perturbed_str_forms:
    #         if self.nodes[node]["composite"]:
    #             composite_nodes.append(node)
    #         else:
    #             single_nodes.append(node)
        
    #     num_nodes_selected = len(single_nodes) + len(composite_nodes)
    #     num_new = 0
    #     while num_new < num_nodes_selected:#전 loop에 비해 새로 찾은 것이 있는가?
    #         num_new = num_nodes_selected
    #         set_composite_nodes_LDOI_candidates = set()
    #         for node in single_nodes+composite_nodes:
    #             for successor_node in self.successors(node):
    #                 if not Make_Using_Prime_Implicants.two_expanded_nodes_groups_are_contradict([successor_node], single_nodes+composite_nodes):
    #                     if self.nodes[successor_node]["composite"]:
    #                         set_composite_nodes_LDOI_candidates.add(successor_node)
    #                     else:
    #                         single_nodes.append(successor_node)
    #         #find composite node canalized by nodes
    #         for composite_candidate in set_composite_nodes_LDOI_candidates:
    #             regulators = self.get_regulator_names_of_node(composite_candidate)
    #             if all((regulator in single_nodes) for regulator in regulators):
    #                 composite_nodes.append(composite_candidate)
            
    #         single_nodes = list(set(single_nodes))
    #         composite_nodes = list(set(composite_nodes))
    #         num_nodes_selected = len(single_nodes) + len(composite_nodes)
        
    #     return single_nodes


if __name__ == "__main__":
    import Model_read
    add_test_model = os.path.join(r"D:\new canalizing kernel\우정 tumorigenesis 모델", "Fumia_logic_model.bnet")
    dynamics_test_model = Model_read.read_pyboolnet_file(add_test_model)
    expanded_toy = make_expanded_net_using_dynamics_pyboolnet(dynamics_test_model,reduction=True, perturbation={})
