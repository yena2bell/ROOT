import itertools

import Expanded_net_analysis as ENA
from Cycle_analysis import Find_cycles_containing_the_node
from Attractor_landscape_calculation import Attractor_landscape_for_specific_IC

import iATG_module

class ITP:
    def __init__(self, 
                 iCA,
                 attr_basal_irrev,
                 attr_transition,
                 attr_basal_rev):
        self.iCA = iCA
        self.attr_basal_irrev = attr_basal_irrev
        self.attr_transition = attr_transition
        self.attr_basal_rev = attr_basal_rev

        self.non_cyclic_part_basal_rev = self._get_non_cyclic_part_of_attractor(self.attr_basal_rev)
        self.non_cyclic_part_basal_irrev = self._get_non_cyclic_part_of_attractor(self.attr_basal_irrev)
        self.non_cyclic_part_transition = self._get_non_cyclic_part_of_attractor(self.attr_transition)

        self.phenotype_nodes = []

        self.irreversibility_motifs = []
        self.coherency_conditions = []

        self.pheno_basal_irrev = {}
        self.pheno_transition = {}
        self.pheno_basal_rev = {}
        

    def __repr__(self):
        return "ITP from {} to {} to {}".format(self.attr_basal_irrev, self.attr_transition, self.attr_basal_rev)
    
    def _get_non_cyclic_part_of_attractor(self, attr_tuple_form):
        return self._get_attractor_object_using_attr_tuple_form(attr_tuple_form).extract_non_cyclic_states()
    
    def set_phenotype_nodes(self, phenotype_nodes):
        self.phenotype_nodes = list(phenotype_nodes)
        self._calculate_phenotype_scores()

    def _calculate_phenotype_scores(self):
        attractor_basal_irrev = self._get_attractor_object_using_attr_tuple_form(self.attr_basal_irrev)
        attractor_transition = self._get_attractor_object_using_attr_tuple_form(self.attr_transition)
        attractor_basal_rev = self._get_attractor_object_using_attr_tuple_form(self.attr_basal_rev)
        
        self.pheno_basal_irrev = {phenonode:averstate for phenonode, averstate in attractor_basal_irrev.get_average_state().items() if phenonode in self.phenotype_nodes}
        self.pheno_transition = {phenonode:averstate for phenonode, averstate in attractor_transition.get_average_state().items() if phenonode in self.phenotype_nodes}
        self.pheno_basal_rev = {phenonode:averstate for phenonode, averstate in attractor_basal_rev.get_average_state().items() if phenonode in self.phenotype_nodes}

    def reset_the_info_of_irreversibility_kernel(self):
        self.irreversibility_motifs = []
        self.coherency_conditions = []

    def get_phenotype_score(self):
        print("attractor_basal_irrev: ", self.pheno_basal_irrev)
        print("->")
        print("attractor_transition: ", self.pheno_transition)
        print("->")
        print("attractor_basal_rev: ", self.pheno_basal_rev)

    def find_irreversibility_kernel(self, max_len=0):
        fixed_nodes = list(self.iCA.iATG.fixed_node_state_map.keys())
        expanded_subnet, expanded_nodes_basal_rev_only = self._make_expanded_subnet_of_non_cyclic_basal_rev()
        feedback_loops_in_expanded_net = self._find_feedback_loops_from_expanded_net_containing_specific_nodes(expanded_subnet, expanded_nodes_basal_rev_only, max_len, fixed_nodes)
        # these are used to find irreversibility-enhancing positive feedback loops.
        feedback_loops_copied = feedback_loops_in_expanded_net.copy()

        non_cyclic_part_basal_common = self._intersection_between_dict(self.non_cyclic_part_basal_rev, self.non_cyclic_part_basal_irrev)
        irreversibility_enhancing_subnets_and_coherency_conditions_selected = []

        i_comb = 1
        while len(feedback_loops_in_expanded_net)>=i_comb:
            feedback_loops_used_in_selected = []
            irreversibility_enhancing_subnets_and_coherency_conditions, feedback_loops_used_for_each_irreversibility_enhancing_subnet = self._find_irreversibility_enhancing_subnet_and_corresponding_coherency_conditions(feedback_loops_copied, i_comb)
            for i, irreversibility_enhancing_subnet_and_coherency_condition in enumerate(irreversibility_enhancing_subnets_and_coherency_conditions):
                coherency_condition = irreversibility_enhancing_subnet_and_coherency_condition[1]
                if self._dict1_contain_dict2(non_cyclic_part_basal_common, coherency_condition):
                    irreversibility_enhancing_subnets_and_coherency_conditions_selected.append(irreversibility_enhancing_subnet_and_coherency_condition)
                    feedback_loops_used_in_selected.extend(feedback_loops_used_for_each_irreversibility_enhancing_subnet[i])
            feedback_loops_copied = [feedback_loop for feedback_loop in feedback_loops_copied if feedback_loop not in feedback_loops_used_in_selected]
            i_comb += 1
        
        # delete duplcated irreversibility_enhancing_subnets_and_coherency_conditions
        irreversibility_enhancing_subnets_and_coherency_conditions_not_duplicated = []
        for irreversibility_enhancing_subnet_and_coherency_condition in irreversibility_enhancing_subnets_and_coherency_conditions_selected:
            if irreversibility_enhancing_subnet_and_coherency_condition not in irreversibility_enhancing_subnets_and_coherency_conditions_not_duplicated:
                irreversibility_enhancing_subnets_and_coherency_conditions_not_duplicated.append(irreversibility_enhancing_subnet_and_coherency_condition)
            
        irreversibility_enhancing_subnets_and_coherency_conditions_not_superset = self._delete_irreversibility_enhancing_subnet_which_is_superset(irreversibility_enhancing_subnets_and_coherency_conditions_not_duplicated)
        
        self._convert_selected_irreversibility_enhancing_subnets_to_irreversibility_kernel(irreversibility_enhancing_subnets_and_coherency_conditions_not_superset)
    
    def _delete_irreversibility_enhancing_subnet_which_is_superset(self, irreversibility_enhancing_subnets_and_coherency_conditions):
        """delete irreversibility_enhancing_subnet which is superset of other irreversibility_enhancing_subnet"""
        changed = True

        while changed:
            changed = False
            to_remove = []

            for i,j in itertools.permutations(range(len(irreversibility_enhancing_subnets_and_coherency_conditions)), 2):
                irreversibility_enhancing_subnet_and_coherency_condition_i = irreversibility_enhancing_subnets_and_coherency_conditions[i]
                irreversibility_enhancing_subnet_i = irreversibility_enhancing_subnet_and_coherency_condition_i[0]
                irreversibility_enhancing_subnet_and_coherency_condition_j = irreversibility_enhancing_subnets_and_coherency_conditions[j]
                irreversibility_enhancing_subnet_j = irreversibility_enhancing_subnet_and_coherency_condition_j[0]

                if irreversibility_enhancing_subnet_i == irreversibility_enhancing_subnet_j:
                    coherency_condition_i = irreversibility_enhancing_subnet_and_coherency_condition_i[1]
                    coherency_condition_j = irreversibility_enhancing_subnet_and_coherency_condition_j[1]
                    if self._dict1_contain_dict2(coherency_condition_i, coherency_condition_j):
                        to_remove.append(irreversibility_enhancing_subnet_and_coherency_condition_i)
                        break

                elif irreversibility_enhancing_subnet_i.issuperset(irreversibility_enhancing_subnet_j):
                    to_remove.append(irreversibility_enhancing_subnet_and_coherency_condition_i)
                    break
            
            if to_remove:
                irreversibility_enhancing_subnets_and_coherency_conditions = [subnet_and_condition for subnet_and_condition in irreversibility_enhancing_subnets_and_coherency_conditions if subnet_and_condition not in to_remove]
                changed = True
        
        return irreversibility_enhancing_subnets_and_coherency_conditions

    def _convert_selected_irreversibility_enhancing_subnets_to_irreversibility_kernel(self, selected_irreversibility_enhancing_subnets):
        for subnet_and_coherency_condition in selected_irreversibility_enhancing_subnets:
            expanded_node_names_in_subnet = subnet_and_coherency_condition[0]
            coherency_condition = subnet_and_coherency_condition[1]

            nodes_in_subnet = set()
            for expanded_node_name in expanded_node_names_in_subnet:
                original_node_name = self._get_original_node_name_from_expanded_node_name(expanded_node_name)
                nodes_in_subnet.add(original_node_name)
            self.irreversibility_motifs.append(nodes_in_subnet)
            
            self.coherency_conditions.append(coherency_condition)
        

    def _get_original_node_name_from_expanded_node_name(self, expanded_node_name):
        """only when the expanded node is not composite"""
        dict_form_expanded_node = self._get_info_of_expanded_node(expanded_node_name).dict_form
        return next(iter(dict_form_expanded_node.keys()))

    def _get_expanded_network(self):
        return self.iCA.iATG.get_expanded_network()
    
    def _convert_dict_to_str_form(self, dict_form):
        """dict_form is {node:state}
        convert this to str form. the str form is used to refer the node in expanded network"""
        return ENA.Expanded_node.get_str_form_from_dict_form(dict_form)
    
    @staticmethod
    def _difference_between_dict(dict1, dict2):
        """Creates and returns a new dictionary containing key-value pairs from dict1  
        that are not part of the intersection between dict1 and dict2."""
        diff_dict = {}
        for key, val in dict1.items():
            if key in dict2:
                if dict2[key] != val:
                    diff_dict[key] = val
            else:
                diff_dict[key] = val
        
        return diff_dict
    
    @staticmethod
    def _intersection_between_dict(dict1, dict2):
        """Creates and returns a new dictionary containing only the key-value pairs  
        that match in both dict1 and dict2."""  
        inter_dict = {}
        for key, val in dict1.items():
            if key in dict2:
                if dict2[key] == val:
                    inter_dict[key] = val
        
        return inter_dict
    
    @staticmethod
    def _dict1_contain_dict2(dict1, dict2):
        """Returns True if dict1 contains dict2, False otherwise."""
        for key, val in dict2.items():
            if key not in dict1:
                return False
            if dict1[key] != val:
                return False
        
        return True
    
    def _get_attractor_object_using_attr_tuple_form(self, attr_tuple_form):
        return self.iCA.iATG.get_attractor_using_attr_tuple_form(attr_tuple_form)

    def _make_expanded_subnet_of_non_cyclic_basal_rev(self):
        """Creates an expanded subnetwork corresponding to the non-cyclic  
        portion of attractor_basal_rev.  

        Feedback loops found in this subnetwork represent positive feedback loops  
        that satisfy coherency in the non-cyclic portion of attractor_basal_rev.
        (i.e. irreversibility enhancing positive feedback loops)

        When searching for feedback loops in expanded_subnet, each loop must contain  
        at least one expanded node with a node-state pair different from attractor_basal_irrev.  

        Returns both the expanded subnetwork and a list of names of such expanded nodes."""  
        expanded_network = self._get_expanded_network()

        expanded_nodes_non_cyclic_basal_rev = []
        for node_in_basal_rev, state_in_basal_rev in self.non_cyclic_part_basal_rev.items():
            dict_form_of_expanded_node = {node_in_basal_rev: state_in_basal_rev}
            str_form_expanded_node = self._convert_dict_to_str_form(dict_form_of_expanded_node)
            expanded_nodes_non_cyclic_basal_rev.append(str_form_expanded_node)
        
        expanded_subnet = expanded_network.get_subnetwork(expanded_nodes_non_cyclic_basal_rev)
        # Determine the expanded network subnetwork corresponding to the non-cyclic part of attractor_basal_rev.  
        # The feedback loops identified here are the positive feedback loops that satisfy coherency in the non-cyclic nodes of attractor_basal_rev.  

        expanded_nodes_basal_rev_only = []        
        non_cyclic_part_basal_rev_only = self._difference_between_dict(self.non_cyclic_part_basal_rev, self.non_cyclic_part_basal_irrev)

        for node_in_basal_rev_only, state_in_basal_rev_only in non_cyclic_part_basal_rev_only.items():
            dict_form_of_expanded_node = {node_in_basal_rev_only:state_in_basal_rev_only}
            str_form_expanded_node = self._convert_dict_to_str_form(dict_form_of_expanded_node)
            expanded_nodes_basal_rev_only.append(str_form_expanded_node)
        # Irreversibility-enhancing positive feedback loops must include at least one  
        # node-state pair (expanded node) from attractor_basal_rev that is different  
        # from attractor_basal_irrev.  
        # In other words, a feedback loop should be found in the expanded_subnet  
        # that contains at least one expanded node from expanded_nodes_basal_rev_only. 

        return expanded_subnet, expanded_nodes_basal_rev_only

    def _find_feedback_loops_from_expanded_net_containing_specific_nodes(self, expanded_net, specific_expanded_nodes=[], max_len=0, fixed_nodes=[]):
        """Finds and returns feedback loops in the given expanded network  
        that contain at least one expanded node from specific_expanded_nodes.
        the found feedback loops can be used to find irreversibility-enhancing positive feedback loops.
        positive feedback loops should not contain any fixed nodes in the network.
        
        max_len: maximum length of the feedback loop to be searched.
        If max_len is 0, the length of the feedback loop is not limited."""
        feedback_loops = []
        list_of_edges = []

        expanded_nodes_in_list_of_edges_wo_specific_node = set()
        for edge in expanded_net.edges():
            list_of_edges.append(edge)
            expanded_nodes_in_list_of_edges_wo_specific_node.update(edge)

        specific_expanded_nodes = set(specific_expanded_nodes)
        specific_expanded_nodes.intersection_update(expanded_nodes_in_list_of_edges_wo_specific_node)

        while specific_expanded_nodes:
            specific_expanded_node = specific_expanded_nodes.pop()
            cycle_finder = Find_cycles_containing_the_node(specific_expanded_node, list_of_edges)
            if max_len == 0:
                cycles_containing_the_node = cycle_finder.find_cycles(algorithm="Johnson", max_len=None, return_node_form=True)
            else:
                # to limit the maximum length of the feedback loop to be searched,
                # use the simple algorithm
                cycles_containing_the_node = cycle_finder.find_cycles(algorithm="simple", max_len=max_len, return_node_form=True)

            feedback_loops.extend(cycles_containing_the_node)

            #remake the list of edges without the specific_expanded_node
            list_of_edges_wo_specific_node = []
            expanded_nodes_in_list_of_edges_wo_specific_node = set()
            for edge in list_of_edges:
                if specific_expanded_node in edge:
                    continue
                list_of_edges_wo_specific_node.append(edge)
                expanded_nodes_in_list_of_edges_wo_specific_node.update(edge)
            list_of_edges = list_of_edges_wo_specific_node
            specific_expanded_nodes.intersection_update(expanded_nodes_in_list_of_edges_wo_specific_node)
            # by the process of deleting edges, specific_expanded_node not yet popped can be deleted from the list of edges

        
        # filter feedback loops that do not contain fixed nodes
        feedback_loops_not_containing_fixed_nodes = []
        for feedback_loop in feedback_loops:
            for expanded_node in feedback_loop:
                node_info = self._get_info_of_expanded_node(expanded_node)
                if not node_info.is_composite():
                    if  set(node_info.dict_form).intersection(fixed_nodes):
                        break
            else:
                feedback_loops_not_containing_fixed_nodes.append(feedback_loop)
        
        return feedback_loops_not_containing_fixed_nodes
                    
    def _find_irreversibility_enhancing_subnet_and_corresponding_coherency_conditions(self, feedback_loops_in_expanded_net, i_comb):
        """by combining i_comb number of irreversibility-enhancing positive feedback loops (feedback_loops_in_expanded_net),
        find irreversibility-enhancing subnet and corresponding coherency condition"""
        irreversibility_enhancing_subnets_and_coherency_conditions = []
        feedback_loops_used_for_each_irreversibility_enhancing_subnet = []

        for feedback_loops_to_combine in itertools.combinations(feedback_loops_in_expanded_net, i_comb):
            nodes_in_irreversibility_enhancing_subnet, coherency_condition = self._combine_irreversibility_enhancing_positive_feedback_loops(feedback_loops_to_combine)
            irreversibility_enhancing_subnets_and_coherency_conditions.append((nodes_in_irreversibility_enhancing_subnet, coherency_condition))
            feedback_loops_used_for_each_irreversibility_enhancing_subnet.append(feedback_loops_to_combine)
        
        return irreversibility_enhancing_subnets_and_coherency_conditions, feedback_loops_used_for_each_irreversibility_enhancing_subnet
        
    def _combine_irreversibility_enhancing_positive_feedback_loops(self, feedback_loops_to_combine):
        """Combines irreversibility-enhancing positive feedback loops (feedback_loops_to_combine)  
        to construct an irreversibility-enhancing subnet.  
        Returns the set of nodes that constitute the subnet and its coherency condition."""
        nodes_in_irreversibility_enhancing_subnet = set()
        coherency_condition = {}
        # get set of nodes in the irreversibility-enhancing subnet
        # each element is node name of the expanded network
        # it also contains composite nodes
        for feedback_loop in feedback_loops_to_combine:
            nodes_in_irreversibility_enhancing_subnet.update(feedback_loop)
        
        composite_nodes = set()
        for node_name in nodes_in_irreversibility_enhancing_subnet:
            expanded_node_info = self._get_info_of_expanded_node(node_name)
            if expanded_node_info.is_composite():
                composite_nodes.add(node_name)                
                regulators_of_composite = set(expanded_node_info.get_regulators_of_composite())
                for regulator in regulators_of_composite.difference(nodes_in_irreversibility_enhancing_subnet):
                    coherency_condition = {**coherency_condition, **self._get_info_of_expanded_node(regulator).dict_form}
        
        nodes_in_irreversibility_enhancing_subnet.difference_update(composite_nodes)
        # through this process, the composite node is removed from the set of nodes 
        # in the set of nodes in irreversibility-enhancing subnet
        
        return nodes_in_irreversibility_enhancing_subnet, coherency_condition
    
    def _get_info_of_expanded_node(self, expanded_node_name):
        return self._get_expanded_network().nodes()[expanded_node_name]['info']
    
    def get_resetting_efficacy_score_of(self, control_configuration, verbose=False):
        """control configuration is given as a dictionary in the form of {node: state}.
        When the model state is in attractor_basal_rev, 
        the system applies the control configuration to let the state converge to a new attractor. 
        Then, after removing control configuration, 
        it calculates the reversion probability, 
        which is the probability that the system returns to attractor_basal_irrev.
        
        If attractor_basal_rev is a cyclic attractor or 
        the new attractor obtained after applying control configuration is a cyclic attractor, 
        the reversion probability is calculated under the assumption 
        that the system is equally likely to be in any state of the cyclic attractor."""
        dynamics_pyboolnet = self.iCA.iATG.dynamics_Boolean_net
        IC_basal = self.iCA.iATG.IC_basal
        fixed_node_state_map = self.iCA.iATG.fixed_node_state_map

        fixed_node_state_map_with_control_configuration = {**fixed_node_state_map, **control_configuration}
        attractor_landscape_wo_control = self.iCA.iATG.attractor_landscape_basal
        attractor_landscape_controlled = Attractor_landscape_for_specific_IC(dynamics_pyboolnet, IC_basal, fixed_node_state_map_with_control_configuration)

        attractor_basal_rev_object = self._get_attractor_object_using_attr_tuple_form(self.attr_basal_rev)
        attractor_basal_irrev_object = self._get_attractor_object_using_attr_tuple_form(self.attr_basal_irrev)

        attractor_states_of_basal_rev = attractor_basal_rev_object.get_attractor_states()
        num_of_attractor_basal_rev_states = len(attractor_states_of_basal_rev)
        if verbose:
            print("this attractor has {} states".format(num_of_attractor_basal_rev_states))
        reversion_probability = 0
        for state_object in attractor_states_of_basal_rev:
            if verbose:
                print("state_object: ")
                print(state_object)
            attractor_object_under_control = attractor_landscape_controlled.converge_network_state_to_attractor(state_object)
            attractor_states_converged_after_control = attractor_object_under_control.get_attractor_states()
            num_of_attractor_converged_after_control_states = len(attractor_states_converged_after_control)
            if verbose:
                print("this state converges to attractor with {} states after applying the temporary control".format(num_of_attractor_converged_after_control_states))
            num_of_states_going_to_basal_irrev = 0
            for state_converged_after_control in attractor_states_converged_after_control:
                attractor_object_after_removing_control = attractor_landscape_wo_control.converge_network_state_to_attractor(state_converged_after_control)
                if attractor_object_after_removing_control == attractor_basal_irrev_object:
                    num_of_states_going_to_basal_irrev += 1
                    if verbose:
                        print("the state in converged attractor")
                        print(state_converged_after_control)
                        print("arrived to attractor_basal_irrev")
            reversion_probability_for_this_state = num_of_states_going_to_basal_irrev / num_of_attractor_converged_after_control_states
            reversion_probability += (reversion_probability_for_this_state / num_of_attractor_basal_rev_states)
        
        return reversion_probability
    
    def get_the_average_state_for_each_IC_by_applying(self, control_configuration, verbose=False):
        """control_configuration is given as a dictionary in the form of {node: state}.
        When the model state is in attractor_basal_rev,
        applying the control_configuration causes the system to converge to a new attractor.
        From there, we identify the iCA(s) in the iATG (constructed under the control_configuration)
        that the newly converged attractor would reach via attractor transition.
        Then, we calculate the average state per input configuration (IC) for those iCAs.
        
        If the newly converged attractor leads to multiple iCAs via transition,
        we use the transition probabilities (TPs) in the iATG to determine the contribution of each iCA,
        and compute a weighted average accordingly."""
        dynamics_pyboolnet = self.iCA.iATG.dynamics_Boolean_net
        node_names = dynamics_pyboolnet.get_node_names()
        IC_basal = self.iCA.iATG.IC_basal
        IC_transition = self.iCA.iATG.IC_transition
        fixed_node_state_map = self.iCA.iATG.fixed_node_state_map
        fixed_node_and_control_configuration = {**fixed_node_state_map, **control_configuration}

        attractor_basal_rev_object = self._get_attractor_object_using_attr_tuple_form(self.attr_basal_rev)
        attractor_states_of_basal_rev = attractor_basal_rev_object.get_attractor_states()
        num_of_attractor_basal_rev_states = len(attractor_states_of_basal_rev)
        if verbose:
            print("attractor_basal_rev has {} states".format(num_of_attractor_basal_rev_states))            

        basal_average_state = {node_name: 0 for node_name in node_names}
        transition_average_state = {node_name: 0 for node_name in node_names}
        for state_object in attractor_states_of_basal_rev:
            if verbose:
                print("state_object: ")
                print(state_object)
            iATG_controlled = iATG_module.iATG(dynamics_pyboolnet, IC_basal, IC_transition, fixed_node_and_control_configuration)
            iATG_controlled.set_empty_attractor_landscape_for_each_IC_wo_calculation()
            attractor_landscape_basal_in_controlled = iATG_controlled.attractor_landscape_basal
            attractor_landscape_transition_in_controlled = iATG_controlled.attractor_landscape_transition

            attractor_converged = attractor_landscape_basal_in_controlled.converge_network_state_to_attractor(state_object)
            attractor_landscape_basal_in_controlled.attractor_index_map[0] = attractor_converged
            if verbose:
                print("in the control configuration, this state converged to attractor")
                print(attractor_converged.get_average_state())

            iATG_controlled.get_attractor_transitions_induced_by_IC_change_and_calculate_TPs()
            num_of_att_in_att_land_basal_in_controlled = len(attractor_landscape_basal_in_controlled.attractor_index_map)
            num_of_att_in_att_land_transition_in_controlled = len(attractor_landscape_transition_in_controlled.attractor_index_map)
            
            # To calculate the proportion of flow from attractor_converged into each iCA,
            # we use the code for computing iCA_sizes.
            # At this point, by overwriting the attindex_basinratio_map and 
            # setting the basin ratio of only attractor_converged to 1,
            # the size of each iCA corresponds to the proportion of flow 
            # that originates from attractor_converged and 
            # enters the iCA according to the transition process (TP).
            attractor_landscape_basal_in_controlled.attindex_basinratio_map = {i:0 for i in range(num_of_att_in_att_land_basal_in_controlled)}
            attractor_landscape_transition_in_controlled.attindex_basinratio_map = {i:0 for i in range(num_of_att_in_att_land_transition_in_controlled)}
            attractor_landscape_basal_in_controlled.attindex_basinratio_map[0] = 2
            # in the method 'find_iCAs_and_calculate_iCA_sizes' basin ratio for each attractor landscape is divided by 2.
            # to avoid by this division, we set the basin ratio of attractor_converged to 2.
            iATG_controlled.find_iCAs_and_calculate_iCA_sizes()
            if verbose:
                print("starting from this newly converged attractor,")
                print("attractor transitions induce the following iCAs:")
                iATG_controlled.show_iCAs_and_their_sizes()

            basal_average_state_for_att = {node_name: 0 for node_name in node_names}
            transition_average_state_for_att = {node_name: 0 for node_name in node_names}
            for iCA_object in iATG_controlled.iCAs:
                iCA_object.set_phenotype_nodes(node_names)
                iCA_size = iCA_object.get_iCA_size()
                basal_average_state_of_iCA = iCA_object.get_phenotype_for_IC("basal")
                transition_average_state_of_iCA = iCA_object.get_phenotype_for_IC("transition")

                if verbose:
                    print("iCA: ", iCA_object)
                    print("iCA size: ", iCA_size)
                    print("basal average state of iCA: ", basal_average_state_of_iCA)
                    print("transition average state of iCA: ", transition_average_state_of_iCA)

                for node_name in node_names:
                    basal_average_state_for_att[node_name] += iCA_size * basal_average_state_of_iCA[node_name]
                    transition_average_state_for_att[node_name] += iCA_size * transition_average_state_of_iCA[node_name]

            for node_name in node_names:               
                basal_average_state[node_name] += basal_average_state_for_att[node_name] / num_of_attractor_basal_rev_states
                transition_average_state[node_name] += transition_average_state_for_att[node_name] / num_of_attractor_basal_rev_states
        
        IC_averagestate_map = {"basal": basal_average_state, "transition": transition_average_state}
        return IC_averagestate_map
    
    def get_reversing_efficacy_score_of(self, control_configuration,
                                        phenotype_nodes=None):
        """Calculate the Hamming distance between the phenotype node states of the attractor_basal_irrev in this ITP
        and the average phenotype node states of 'basal' after applying control_configuration.
        
        Also calculate the Hamming distance between the phenotype node states of the attractor_transition in this ITP
        and the average phenotype node states of 'transition' after applying control_configuration.
        
        Return the sum of these two Hamming distances.
        If phenotype_nodes is None, use the phenotype nodes of this ITP."""
        
        IC_averagestate_map = self.get_the_average_state_for_each_IC_by_applying(control_configuration)
        
        desired_phenotype_score = 1/(1+self.get_phenotype_score_compared_to_average_state_for_each_IC(IC_averagestate_map, phenotype_nodes))
        return desired_phenotype_score
    
    def get_phenotype_score_compared_to_average_state_for_each_IC(self, 
                                                                      IC_averagestate_map, 
                                                                      phenotype_nodes=None):
        """Compute the Hamming distance between the phenotype node states 
        of the attractor_basal_irrev in this ITP 
        and the average phenotype node states of 'basal' in IC_averagestate_map.
        
        Also compute the Hamming distance between the phenotype node states 
        of the attractor_transition in this ITP 
        and the average phenotype node states of 'transition' in IC_averagestate_map.

        Return the sum of these two Hamming distances.
        
        if phenotype_nodes is None,
        use the phenotype nodes of this ITP."""
        if phenotype_nodes is None:
            phenotype_nodes = self.phenotype_nodes
        
        attractor_basal_irrev = self._get_attractor_object_using_attr_tuple_form(self.attr_basal_irrev)
        attractor_transition = self._get_attractor_object_using_attr_tuple_form(self.attr_transition)
        
        Hamming_distance_basal = self._get_Hamming_distance(attractor_basal_irrev.get_average_state(),
                                                            IC_averagestate_map["basal"],
                                                            phenotype_nodes)
        Hamming_distance_transition = self._get_Hamming_distance(attractor_transition.get_average_state(),
                                                                  IC_averagestate_map["transition"],
                                                                    phenotype_nodes)
        
        return Hamming_distance_basal + Hamming_distance_transition
    
    @staticmethod
    def _get_Hamming_distance(state1:dict, state2:dict, nodes_considered):
        """Calculate the Hamming distance between two states.
        Only consider the nodes in phenotype_nodes."""
        distance = 0
        for node in nodes_considered:
            distance += abs(state1[node] - state2[node])
        
        return distance
    
    def get_control_strategiess_having_n_control_targets(self, n_nodes=1):
        """Control candidates consisting of n control target nodes are identified.
        From the nodes included in the irreversibility motifs, 
        n nodes are selected as control targets,
        and each target node is assigned a state opposite to its state in attractor_basal_rev.

        Each control candidate must share at least one node with each irreversibility motif 
        and the set of control targets.
        
        This process should be performed after the irreversibility kernel has been analyzed."""
        nodes_in_irreversibility_motifs = set()
        for irreversibility_motif in self.irreversibility_motifs:
            nodes_in_irreversibility_motifs.update(irreversibility_motif)
        
        nodes_in_coherency_conditions = set()
        for coherency_condition in self.coherency_conditions:
            nodes_in_coherency_conditions.update(coherency_condition)
        
        control_candidates = []
        for control_targets in itertools.combinations(nodes_in_irreversibility_motifs.union(nodes_in_coherency_conditions), n_nodes):
            
            check_intersection_to_all_motifs_or_conditions = True
            for i, irreversibility_motif in enumerate(self.irreversibility_motifs):
                coherency_condition = self.coherency_conditions[i]
                condition_and_motif = irreversibility_motif.union(coherency_condition)
                if not condition_and_motif.intersection(control_targets):
                    check_intersection_to_all_motifs_or_conditions = False
                    break
            if not check_intersection_to_all_motifs_or_conditions:
                continue

            control_candidate = {}
            for node in control_targets:
                control_candidate[node] = 1 - self.non_cyclic_part_basal_rev[node]
                
            control_candidates.append(control_candidate)
        return control_candidates

    def find_reverse_controls(self):
        reverse_controls = []
        for i, irreversibility_motif in enumerate(self.irreversibility_motifs):
            coherency_condition = self.coherency_conditions[i]
            reverse_controls.extend(self._find_reverse_controls_for_given_irreversibility_kernel(irreversibility_motif, coherency_condition))
        
        return reverse_controls

    def _find_reverse_controls_for_given_irreversibility_kernel(self, irreversibility_motif, coherency_condition):
        controls_first = []
        for nodes in power_set(irreversibility_motif):
            control = {}
            for node in nodes:
                control[node] = 1 - self.non_cyclic_part_basal_rev[node]
            controls_first.append(control)
        
        controls_second = []
        for control in controls_first:
            for node, state_value in control.items():
                if self.non_cyclic_part_basal_irrev.get(node, -1) != state_value:
                    break
            else:
                controls_second.append(control)
                
        controls_third = []
        phenotype_nodes_in_transition_non_cyclic = {pheno_node: self.non_cyclic_part_transition[pheno_node] for pheno_node in self.phenotype_nodes if pheno_node in self.non_cyclic_part_transition}

        expanded_net = self._get_expanded_network()
        IC_transition = self.iCA.iATG.IC_transition
        fixed_node_state_map = self.iCA.iATG.fixed_node_state_map
        expanded_nodes_for_IC_transition_and_fixed_nodes = []
        for node, state in IC_transition.items():
            expanded_nodes_for_IC_transition_and_fixed_nodes.append(self._convert_dict_to_str_form({node:state}))
        for node, state in fixed_node_state_map.items():
            expanded_nodes_for_IC_transition_and_fixed_nodes.append(self._convert_dict_to_str_form({node:state}))

        for control in controls_second:
            IC_transition_fixed_nodes_and_control = expanded_nodes_for_IC_transition_and_fixed_nodes.copy()
            for node, state in control.items():
                IC_transition_fixed_nodes_and_control.append(self._convert_dict_to_str_form({node:state}))
            LDOI_by_control = expanded_net.get_LDOI(IC_transition_fixed_nodes_and_control, only_return_single_nodes=True)
            
            phenotype_node_state_map_in_LDOI = {}
            for expanded_node_in_LDOI in LDOI_by_control:
                original_node_name = self._get_original_node_name_from_expanded_node_name(expanded_node_in_LDOI)
                if original_node_name in self.phenotype_nodes:
                    state_in_LDOI = self._get_info_of_expanded_node(expanded_node_in_LDOI).dict_form[original_node_name]
                    phenotype_node_state_map_in_LDOI[original_node_name] = state_in_LDOI
            
            if phenotype_node_state_map_in_LDOI == phenotype_nodes_in_transition_non_cyclic:
                controls_third.append(control)
    
        return controls_third
    
    def check_the_effect_of_reverse_control(self, reverse_control, 
                                            use_all_initials=False,
                                                   waiting_num=1000, 
                                                   difference_threshold=0.00001,
                                                   verbose=False):

        fixed_node_state_map = self.iCA.iATG.fixed_node_state_map
        fixed_node_state_map_with_control = {**fixed_node_state_map, **reverse_control}
        IC_basal = self.iCA.iATG.IC_basal
        IC_transition = self.iCA.iATG.IC_transition
        dynamics_Boolean_network = self.iCA.iATG.dynamics_Boolean_net
        phenotype_nodes = self.phenotype_nodes

        iatg_controlled = iATG_module.iATG(dynamics_Boolean_network, IC_basal, IC_transition, fixed_node_state_map_with_control)
        iatg_controlled.set_phenotype_nodes(phenotype_nodes)
        iatg_controlled.calculate_attractor_landscapes_for_each_IC(use_all_initials,
                                                   waiting_num=waiting_num, 
                                                   difference_threshold=difference_threshold,
                                                   verbose=False)
        iatg_controlled.get_attractor_transitions_induced_by_IC_change_and_calculate_TPs()
        iatg_controlled.find_iCAs_and_calculate_iCA_sizes()
        
        attractor_landscape_basal_in_controlled = iatg_controlled.attractor_landscape_basal
        attractor_basal_irrev = self._get_attractor_object_using_attr_tuple_form(self.attr_basal_irrev)
        for i, attractor_in_controlled in attractor_landscape_basal_in_controlled.attractor_index_map.items():
            if attractor_basal_irrev == attractor_in_controlled:
                break
        else:
            raise ValueError("The attractor_basal_irrev in the controlled iATG is not found.")
        # this ensures that the condition of preserving attractor`_basal_irrev
        attr_tuple_form_basal_irrev_in_controlled = ('basal',i)

        for ica in iatg_controlled.iCAs:
            if attr_tuple_form_basal_irrev_in_controlled in ica.attractors_in_iCA:
                break
        else:
            raise ValueError("The attractor_basal_irrev in the controlled iATG is not contained in iCA.")
        
        ica.set_phenotype_nodes(phenotype_nodes)
        phenotype_of_ica_in_IC_transition_in_controlled = ica.get_phenotype_for_IC("transition")
        phenotype_of_ica_in_IC_transition_in_controlled = {phenonode:averstate for phenonode, averstate in phenotype_of_ica_in_IC_transition_in_controlled.items() if averstate in (0,1)}
        #extract non-cyclic nodes among phenotype nodes

        phenotype_nodes_in_transition_non_cyclic = {pheno_node: self.non_cyclic_part_transition[pheno_node] for pheno_node in self.phenotype_nodes if pheno_node in self.non_cyclic_part_transition}

        if phenotype_of_ica_in_IC_transition_in_controlled == phenotype_nodes_in_transition_non_cyclic:
            return True
        else:
            raise ValueError("The phenotype under IC_transition of iCA \ncontaining the attr_basal_irrev\n is different form phenotype of attr_transition before control")
        





                






def power_set(input_set, except_empty_set=True):
    """Generate the power set of the given set."""
    input_list = list(input_set)
    power_set_result = []
    
    if except_empty_set:
        start_comb = 1
    else:
        start_comb = 0
        
    for r in range(start_comb, len(input_list) + 1):
        power_set_result.extend(itertools.combinations(input_list, r))
    return power_set_result