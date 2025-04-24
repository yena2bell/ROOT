from itertools import zip_longest

from Attractor_landscape_calculation import Attractor_landscape_for_specific_IC
from SCC_decomposition import SCC_decomposition
from iCA_module import iCA
from Expanded_net_analysis import make_expanded_net_using_dynamics_pyboolnet

class iATG:
    """To construct the input configuration change-derived attractor transition graph (iATG), 
    a specific Boolean network must be given, 
    along with a basal input configuration (IC_basal) 
    and a transition input configuration (IC_transition).
    
    input configuration change-derived attractor transition graph (iATG)
    이것을 구축하기 위해서는 특정한 Boolean network가 주어지고,
    basal input configuration(IC_basal)과 transition input configuration(IC_transition)이 주어져야 한다."""
    def __init__(self, dynamics_Boolean_net:"Dynamics_pyBoolnet", 
                 IC_basal:dict, IC_transition:dict,
                 fixed_node_state_map:dict={}):
        self.dynamics_Boolean_net = dynamics_Boolean_net
        self.IC_basal = IC_basal
        self.IC_transition = IC_transition
        self.fixed_node_state_map = fixed_node_state_map
        self.phenotype_nodes = []

        self.attractor_landscape_basal = None
        self.attractor_landscape_transition = None

        self.attractor_transitions_induced_by_IC_change = []
        self.TPs = {}
        #key is (IC_of_source, source_attractor_index),
        #value is a dictionary whose key is (IC_of_target, target_attractor_index)
        #and value is the transition probability (TP) from the source attractor to the target attractor.
        #here, IC_of_source, IC_of_target is either string "basal" or "transition"
        # such (IC_of_source, source_attractor_index) is named as attr_tuple_form

        self.iCAs = []
        self.expanded_net = None
        self.indexes_of_iCAs_descending_order = None
        self.major_iCAs = []

    def get_attr_tuple_forms(self, get_only_specific_IC=None):
        attr_tuple_forms_basal = [("basal", index) for index in self.attractor_landscape_basal.attractor_index_map.keys()]
        attr_tuple_forms_transition = [("transition", index) for index in self.attractor_landscape_transition.attractor_index_map.keys()]
        if get_only_specific_IC == "basal":
            return attr_tuple_forms_basal
        elif get_only_specific_IC == "transition":
            return attr_tuple_forms_transition
        else:
            return attr_tuple_forms_basal + attr_tuple_forms_transition


    def get_attractor_using_attr_tuple_form(self, attr_tuple_form):
        IC, index = attr_tuple_form
        if IC == "basal":
            return self.attractor_landscape_basal.attractor_index_map[index]
        elif IC == "transition":
            return self.attractor_landscape_transition.attractor_index_map[index]
    
    def get_expanded_network(self):
        if self.expanded_net is None:
            self.expanded_net = make_expanded_net_using_dynamics_pyboolnet(self.dynamics_Boolean_net, reduction=False)
        return self.expanded_net
    
    def set_phenotype_nodes(self, phenotype_nodes):
        self.phenotype_nodes = phenotype_nodes

    def get_phenotype_score_for_IC(self, IC:"basal or transition", iCA_to_consider=[]):
        """Calculates the phenotype score of this iATG for a given input configuration.
        iCA_to_consider에 속한 iCA만 고려하여 phenotype score를 계산한다.
        단, iCA_to_consider가 비어있으면 모든 iCA를 고려한다.
        """
        if not iCA_to_consider:
            iCA_to_consider = self.iCAs
        
        phenotype_score_for_IC = {phenonode:0 for phenonode in self.phenotype_nodes}
        size_sum = 0
        for ica in iCA_to_consider:
            phenotype_score_for_IC_of_ica = ica.get_phenotype_score_for_IC(IC)
            for phenonode in self.phenotype_nodes:
                phenotype_score_for_IC[phenonode] += phenotype_score_for_IC_of_ica[phenonode]*ica.size
            size_sum += ica.size
        
        for phenonode in self.phenotype_nodes:
            phenotype_score_for_IC[phenonode] /= size_sum
            
        
        return phenotype_score_for_IC


    def calculate_attractor_landscapes_for_each_IC(self, use_all_initials=True,
                                                   waiting_num=10000, 
                                                   difference_threshold=0.0001,
                                                   verbose=False):
        """"Calculate the attractor landscape for IC_basal and IC_transition.
        If use_all_initials=True, the attractor landscape is computed using all initial network state values.
        Otherwise, random initial network state values are used.
        In this case, waiting_num and difference_threshold are utilized.
        
        IC_basal과 IC_transition에 대한 attractor landscape를 계산한다.
        이 때, use_all_initials=True이면 모든 initial network state values를 사용하여 attractor landscape를 계산한다.
        그렇지 않으면 random initial network state values를 사용하여 attractor landscape를 계산한다.
        이때, waiting_num과 difference_threshold를 활용한다."""
        if use_all_initials:
            print("Calculating attractor landscapes using all initial states.")
            self.calculate_attractor_landscapes_for_each_IC_using_all_initials()
        else:
            # 이 경우에는 모든 attractors 를 다른 package에서 계산하고,
            # 그 결과를 이용하여 attractor landscape를 approximation 기반으로 계산하는 걸로 하자.
            print("Calculating attractor landscapes using random initial states.")
            self.calcuate_attractor_landscape_for_each_IC_using_random_initials(waiting_num, difference_threshold, verbose)
            pass
    
    def calculate_attractor_landscapes_for_each_IC_using_all_initials(self):
        """Calculate the attractor landscape for IC_basal and IC_transition 
        using all network state values.

        IC_basal과 IC_transition에 대한 attractor landscape를 모든 network state values를 활용하여 계산한다."""
        print("Calculating attractor landscape on input configuration {} using all initial states.".format(self.IC_basal))
        self.attractor_landscape_basal = self.calculate_attractor_landscape_for_specific_IC_using_all_initials(self.IC_basal)
        print("Calculating attractor landscape on input configuration {} using all initial states.".format(self.IC_transition))
        self.attractor_landscape_transition = self.calculate_attractor_landscape_for_specific_IC_using_all_initials(self.IC_transition)

    def calcuate_attractor_landscape_for_each_IC_using_random_initials(self, waiting_num, difference_threshold, verbose):
        """Calculate the attractor landscape for IC_basal and IC_transition 
        using random initial network state values.

        IC_basal과 IC_transition에 대한 attractor landscape를 random initial network state values를 활용하여 계산한다."""
        print("Calculating attractor landscape on input configuration {} using random initial states.".format(self.IC_basal))
        self.attractor_landscape_basal = self.calculate_attractor_landscape_for_specific_IC_using_random_initials(self.IC_basal, waiting_num, difference_threshold, verbose)
        print("Calculating attractor landscape on input configuration {} using random initial states.".format(self.IC_transition))
        self.attractor_landscape_transition = self.calculate_attractor_landscape_for_specific_IC_using_random_initials(self.IC_transition, waiting_num, difference_threshold, verbose)

    def calculate_attractor_landscape_for_specific_IC_using_all_initials(self, IC:dict):
        """Calculate the attractor landscape for a specific input configuration (IC)."""
        attractor_landscape_for_IC = Attractor_landscape_for_specific_IC(self.dynamics_Boolean_net, IC, self.fixed_node_state_map)
        attractor_landscape_for_IC.calculate_attractor_basinratios_from_all_initial_states()

        return attractor_landscape_for_IC
    
    def calculate_attractor_landscape_for_specific_IC_using_random_initials(self, IC:dict, waiting_num, difference_threshold, verbose):
        """Calculate the attractor landscape for a specific input configuration (IC)."""
        attractor_landscape_for_IC = Attractor_landscape_for_specific_IC(self.dynamics_Boolean_net, IC, self.fixed_node_state_map)
        attractor_landscape_for_IC.calculate_attractor_basinratios_from_random_initial_states(waiting_num, difference_threshold, [], verbose)

        return attractor_landscape_for_IC
    
    def get_TP_for_attractor_transition_induced_by_IC_change(self, attr_source_tuple_form, attr_target_tuple_form):
        TPs_from_source_attractor = self.TPs[attr_source_tuple_form]
        TP = TPs_from_source_attractor.get(attr_target_tuple_form, 0)
        return TP

    
    def get_attractor_transitions_induced_by_IC_change_and_calculate_TPs(self):
        """Get the attractor transitions induced by the IC changes.
        Then, calculate the transition probabilities."""
        #initialize
        self.attractor_transitions_induced_by_IC_change = []
        self.TPs = {}
        num_of_attractors_in_IC_basal = len(self.attractor_landscape_basal.attractor_index_map)
        num_of_attractors_in_IC_transition = len(self.attractor_landscape_transition.attractor_index_map)

        while True:
            # While calculating attractor transitions, if new attractors are found, 
            # they are added to the attractor landscape.
            # In such cases, a while loop is used to recompute attractor transitions.
            self._get_attractor_transitions_induced_by_IC_change_and_calculate_TPs(self.attractor_landscape_basal)
            self._get_attractor_transitions_induced_by_IC_change_and_calculate_TPs(self.attractor_landscape_transition)

            num_of_attractors_in_IC_basal_after_calculation = len(self.attractor_landscape_basal.attractor_index_map)
            num_of_attractors_in_IC_transition_after_calculation = len(self.attractor_landscape_transition.attractor_index_map)

            if num_of_attractors_in_IC_basal_after_calculation == num_of_attractors_in_IC_basal:
                if num_of_attractors_in_IC_transition_after_calculation == num_of_attractors_in_IC_transition:
                    break
            
            self.attractor_transitions_induced_by_IC_change = []
            self.TPs = {}
            num_of_attractors_in_IC_basal = num_of_attractors_in_IC_basal_after_calculation
            num_of_attractors_in_IC_transition = num_of_attractors_in_IC_transition_after_calculation


    def _get_attractor_transitions_induced_by_IC_change_and_calculate_TPs(self, attractor_landscape_to_be_source):
        """Get the attractor transitions induced by the IC changes.
        Then, calculate the transition probabilities.
        only consider attractor transitions whose source attractor is in attractor_landscape_to_be_source."""
    
        IC_of_the_attractor_landscape = attractor_landscape_to_be_source.IC
        if IC_of_the_attractor_landscape == self.IC_basal:
            IC_source = "basal"
            IC_target = "transition"
        elif IC_of_the_attractor_landscape == self.IC_transition:
            IC_source = "transition"
            IC_target = "basal"

        for source_attractor_index in attractor_landscape_to_be_source.attractor_index_map.keys():
            attractor_transitions_induced_by_IC_change_from_the_index = self.get_target_attractor_for_given_source_attractor_induced_by_IC_change(source_attractor_index, IC_of_the_attractor_landscape)
            attractor_transition_count_from_source_attractor = {}
            for attractor_transition in attractor_transitions_induced_by_IC_change_from_the_index:
                IC_source_att_index = attractor_transition[0]
                attr_source_tuple_form = (IC_source, IC_source_att_index[1])
                
                IC_target_att_index = attractor_transition[1]
                attr_target_tuple_form = (IC_target, IC_target_att_index[1])

                attractor_transition_hashable = (attr_source_tuple_form, attr_target_tuple_form)
                if attractor_transition_hashable not in self.attractor_transitions_induced_by_IC_change:
                    self.attractor_transitions_induced_by_IC_change.append(attractor_transition_hashable)
                
                attractor_transition_count_from_source_attractor[attr_target_tuple_form] = attractor_transition_count_from_source_attractor.get(attr_target_tuple_form,0) +1
            
            TPs_from_source_attractor = {attr_target_tuple_form: count/len(attractor_transitions_induced_by_IC_change_from_the_index)
                                                              for attr_target_tuple_form, count in attractor_transition_count_from_source_attractor.items()}
            self.TPs[attr_source_tuple_form] = TPs_from_source_attractor
    
    def get_target_attractor_for_given_source_attractor_induced_by_IC_change(self, source_attractor_index, IC_of_source):
        """Given a source attractor, find the target attractor 
        that is induced by the IC change from IC_of_source to IC_of_target.

        source_attractor_index (target_attractor_index) is the key value 
        of 'attractor_index_map' in the self.attractor_landscape_basal or self.attractor_landscape_transition
        when the value is the source_attractor (target_attractor).
        
        the result is returned as a form of 
        [((IC_of_source, source_attractor_index),(IC_of_target, target_attractor_index)),...]
        each tuple in the list ((IC_of_source, source_attractor_index),(IC_of_target, target_attractor_index)) is attractor_transition_induced_by_IC_change
        and the list is attractor_transitions_induced_by_IC_change
        
        The obtained attractor_transitions_induced_by_IC_change may 
        contain duplicate attractor_transition_induced_by_IC_change.
        These duplicate entries will be used later when calculating the transition probability."""
        attractor_transitions_induced_by_IC_change = []
        
        if IC_of_source == self.IC_basal:
            IC_of_target = self.IC_transition
            attractor_landscape_of_source = self.attractor_landscape_basal
            attractor_landscape_of_target = self.attractor_landscape_transition
        elif IC_of_source == self.IC_transition:
            IC_of_target = self.IC_basal
            attractor_landscape_of_source = self.attractor_landscape_transition
            attractor_landscape_of_target = self.attractor_landscape_basal

        attractor = attractor_landscape_of_source.attractor_index_map[source_attractor_index]

        for attractor_state in attractor.attractor_states:
            attractor_state_IC_changed = attractor_state.apply_perturbation(IC_of_target, True)
            target_attractor = attractor_landscape_of_target.converge_network_state_value_to_attractor(attractor_state_IC_changed)
            target_attractor_index = attractor_landscape_of_target.get_attractor_index_from_attractor(target_attractor)

            attractor_transition_induced_by_IC_change = ((IC_of_source, source_attractor_index),(IC_of_target, target_attractor_index))
            attractor_transitions_induced_by_IC_change.append(attractor_transition_induced_by_IC_change)
        
        return attractor_transitions_induced_by_IC_change

    def find_iCAs_and_calculate_iCA_sizes(self):
        """find input configuration change-derived converged attractors (iCAs)
        iCA is a SCC in iATG.
        before find this, self.attractor_transitions_induced_by_IC_change should be calculated.
        
        after finding iCAs, calculate the size of each iCA."""
        SCCs_in_iATG = SCC_decomposition(self.attractor_transitions_induced_by_IC_change)
        #filter SCC composed of single attractor
        for SCC in SCCs_in_iATG:
            if len(SCC) > 1:
                ica = iCA(self, SCC)
                ica.set_phenotype_nodes(self.phenotype_nodes)
                self.iCAs.append(ica)
        
        self._calculate_iCA_size_and_steady_state_probabilities()
    
    def show_iCAs_and_their_sizes(self):
        """Displays iCAs in descending order of iCA size.
        Additionally, shows the phenotype for each IC in each iCA.
        Also presents the cumulative sum of iCA sizes."""
        if self.indexes_of_iCAs_descending_order is None:
            self.indexes_of_iCAs_descending_order = sorted(range(len(self.iCAs)), key=lambda i: self.iCAs[i].get_iCA_size(), reverse=True)

        col_widths = [8,15, 13, 12, 13]
        column_names = ["iCA size", "cumulative sum of iCA sizes", "attractors in this iCA",
                        "phenotype in IC_basal", "phenotype in IC_transition"]
        data = []
        data.append(column_names)
        cumulative_sum_of_iCA_sizes = 0
        for ica_index in self.indexes_of_iCAs_descending_order:
            ica = self.iCAs[ica_index]
            row_of_ica = []

            ica_size = ica.get_iCA_size()
            row_of_ica.append("{:.3f}".format(ica_size))
            cumulative_sum_of_iCA_sizes += ica.get_iCA_size()
            row_of_ica.append("{:.3f}".format(cumulative_sum_of_iCA_sizes))
            row_of_ica.append(str(ica.attractors_in_iCA).strip('[]'))
            row_of_ica.append(str(ica.get_phenotype_score_for_IC('basal')))
            row_of_ica.append(str(ica.get_phenotype_score_for_IC('transition')))
            
            data.append(row_of_ica)

        wrapped_data = [[[value[i:i+col_widths[w_index]] for i in range(0, len(value), col_widths[w_index])] for w_index, value in enumerate(row)] for row in data]
        
        row_heights = [max(len(cell) for cell in row) for row in wrapped_data]
        for row, height in zip(wrapped_data, row_heights):
            for i in range(height):  # Print each wrapped line
                print(" | ".join((cell[i] if i < len(cell) else "").ljust(width)
                                for cell, width in zip(row, col_widths)))
            print("-" * (sum(col_widths) + 3 * (len(col_widths) - 1)))  # Print separator
    
    def select_major_iCAs_according_to_threshold(self, threshold=1):
        """Selects iCAs in descending order of iCA size 
        until the cumulative sum of iCA sizes reaches or exceeds the threshold.
        The selected iCAs are defined as major iCAs."""
        if self.indexes_of_iCAs_descending_order is None:
            self.indexes_of_iCAs_descending_order = sorted(range(len(self.iCAs)), key=lambda i: self.iCAs[i].get_iCA_size(), reverse=True)
            
        cumulative_sum_of_iCA_sizes = 0
        for ica_index in self.indexes_of_iCAs_descending_order:
            ica = self.iCAs[ica_index]
            cumulative_sum_of_iCA_sizes += ica.get_iCA_size()
            self.major_iCAs.append(ica)
            if cumulative_sum_of_iCA_sizes >= threshold:
                break

        

    def _calculate_iCA_size_and_steady_state_probabilities(self):
        """Calculate the size of each iCA.
        the pseudocode is in the supplementary materials of corresponding paper."""
        #initialize the basin_ratio_now and basin_ratio_tmp for each attractor in iATG
        basin_ratios_now = {}
        basin_ratios_tmp = {}
        for index, basin_ratio in self.attractor_landscape_basal.attindex_basinratio_map.items():
            attr_basal_tuple_form = ("basal", index)
            basin_ratios_now[attr_basal_tuple_form] = basin_ratio * 0.5
            basin_ratios_tmp[attr_basal_tuple_form] = 0
        for index, basin_ratio in self.attractor_landscape_transition.attindex_basinratio_map.items():
            attr_transition_tuple_form = ("transition", index)
            basin_ratios_now[attr_transition_tuple_form] = basin_ratio * 0.5
            basin_ratios_tmp[attr_transition_tuple_form] = 0
        
        attr_tuple_forms_not_in_iCA = set(basin_ratios_now.keys())
        for ica in self.iCAs:
            attr_tuple_forms_not_in_iCA.difference_update(ica.attractors_in_iCA)
        

        basin_ratios_not_zero = {key:val for key, val in basin_ratios_now.items() if val != 0}
        while attr_tuple_forms_not_in_iCA.intersection(basin_ratios_not_zero):
            for attr_tuple_form, TPs_from_source_attractor in self.TPs.items():
                for attr_target_tuple_form, TP in TPs_from_source_attractor.items():
                    basin_ratios_tmp[attr_target_tuple_form] += basin_ratios_now[attr_tuple_form]*TP
            
            for attr_tuple_form in basin_ratios_now.keys():
                basin_ratios_now[attr_tuple_form] = 0 #to avoid very small value but not zero
                basin_ratios_now[attr_tuple_form] += basin_ratios_tmp[attr_tuple_form]
                basin_ratios_tmp[attr_tuple_form] = 0
            
            basin_ratios_not_zero = {key:val for key, val in basin_ratios_now.items() if val != 0}
        
        for ica in self.iCAs:
            ica._calculate_size(basin_ratios_now)
            
            #calculate_steady_state_probabilites_in_iCAs
            ica._calculate_steady_state_probabilities()
        
        

