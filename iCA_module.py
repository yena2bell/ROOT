import numpy as np

from ITP_module import ITP

class iCA:
    def __init__(self, iATG, attr_tuple_forms_in_iCA):
        self.iATG = iATG
        self.phenotype_nodes = []
        self.attractors_in_iCA = attr_tuple_forms_in_iCA
        self.transitions_in_iCA = []
        self.get_transitions_in_iCA(self.attractors_in_iCA)
        self.size = 0

        self.attractor_steadystateprob_map = {}
        #key is attr_tuple_form, value is steady-state probability
        self.ITPs = []

    
    def __repr__(self):
        return "iCA composed of {}".format(self.attractors_in_iCA)
    
    def set_phenotype_nodes(self, phenotype_nodes):
        self.phenotype_nodes = phenotype_nodes
    
    def get_transitions_in_iCA(self, attr_tuple_forms_in_iCA):
        for attractor_transition in self.iATG.attractor_transitions_induced_by_IC_change:
            if attractor_transition[0] in attr_tuple_forms_in_iCA and attractor_transition[1] in attr_tuple_forms_in_iCA:
                self.transitions_in_iCA.append(attractor_transition)
    
    def get_iCA_size(self):
        return self.size
    
    def _calculate_size(self, attrtupleform_basinratio_map):
        for attr_tuple_form in self.attractors_in_iCA:
            self.size += attrtupleform_basinratio_map[attr_tuple_form]
    
    def _calculate_steady_state_probabilities(self):
        """For each attractor in iCA, 
        assuming a Markov process where the transition probabilities (TPs) in iATG 
        represent transition probabilities, the steady-state probability is calculated."""
        num_of_attrs = len(self.attractors_in_iCA)
        TP_matrix = np.matrix(np.zeros((num_of_attrs,num_of_attrs)))
        for transition in self.transitions_in_iCA:
            start_index = self.attractors_in_iCA.index(transition[0])
            end_index = self.attractors_in_iCA.index(transition[1])
            TP_matrix[end_index, start_index] = self.iATG.get_TP_for_attractor_transition_induced_by_IC_change(transition[0], transition[1])
        
        eigenvalues, eigenvectors = np.linalg.eig(TP_matrix)
        index_of_eignevalue1 = np.argmin(np.abs(eigenvalues -1))
        # Due to floating-point precision, values like 0.99999999 appear instead of 1,  
        # causing np.where(eigenvalues == 1) to fail.  
        # To resolve this, assuming that at least one eigenvalue is exactly 1,  
        # the code selects the eigenvalue closest to 1.  
        eigenvector_with_eignvalue1 = eigenvectors[:,index_of_eignevalue1]
        steady_state_probabilites = np.transpose(eigenvector_with_eignvalue1).tolist()[0]
        sum_of_probs = sum(steady_state_probabilites)
        steady_state_probabilites = [prob/sum_of_probs for prob in steady_state_probabilites]

        for i, attr_tuple_form in enumerate(self.attractors_in_iCA):
            self.attractor_steadystateprob_map[attr_tuple_form] = steady_state_probabilites[i]
    
    def get_phenotype_for_IC(self, IC:"basal or transition"):
        """Calculates the phenotype of this iCA for a given input configuration.
        To obtain this, the '_calculate_steady_state_probabilities' method must be executed first.

        Among the attractors belonging to this iCA, the attractor associated with the given IC is selected.
        For each such attractor, the phenotype is determined by computing the average node state value of the phenotype node.
        Finally, the phenotypes of all selected attractors are weighted by their steady-state probabilities and averaged.
        """
        attractor_phenotype_map = {}
        for attr_tuple_form in self.attractors_in_iCA:
            if attr_tuple_form[0] == IC:
                attractor = self.iATG.get_attractor_using_attr_tuple_form(attr_tuple_form)
                phenotypenode_average_map_of_attr = {phenonode:averstate for phenonode, averstate in attractor.get_average_state().items() if phenonode in self.phenotype_nodes}
                attractor_phenotype_map[attr_tuple_form] = phenotypenode_average_map_of_attr
        
        phenotypenode_average_map = {phenonode:0 for phenonode in self.phenotype_nodes}
        sum_of_steadystate_probs = 0
        for attr_tuple_form, phenotypenode_average_map_of_attr in attractor_phenotype_map.items():
            steady_state_prob_of_attr = self.attractor_steadystateprob_map[attr_tuple_form]
            sum_of_steadystate_probs += steady_state_prob_of_attr
            for phenonode, averstate in phenotypenode_average_map_of_attr.items():
                phenotypenode_average_map[phenonode] += averstate * steady_state_prob_of_attr
        
        phenotypenode_average_map = {phenonode:averstate/sum_of_steadystate_probs for phenonode, averstate in phenotypenode_average_map.items()}

        return phenotypenode_average_map
            


    def search_ITPs(self):
        """This searches for irreversible transition paths (ITPs) that converge to this iCA.
        An ITP is a path consisting of three vertices (attractors) in the iATG, where:

            The starting attractor is called attractor_basal_irrev.
            The next attractor is called attractor_transition.
            The final attractor is called attractor_basal_rev.

        Attractor_basal_irrev is an attractor in IC_basal that does not belong to this iCA,
        while attractor_transition and attractor_basal_rev are attractors within this iCA."""
        for edge_in_iATG in self.iATG.attractor_transitions_induced_by_IC_change:
            if edge_in_iATG[1] in self.attractors_in_iCA:
                if edge_in_iATG[0] not in self.attractors_in_iCA:
                    if edge_in_iATG[0][0] == "basal":
                        attractor_basal_irrev = edge_in_iATG[0]
                        attractor_transition = edge_in_iATG[1]

                        for transition_in_iCA in self.transitions_in_iCA:
                            if transition_in_iCA[0] == attractor_transition:
                                attractor_basal_rev = transition_in_iCA[1]
                                itp = ITP(self, attractor_basal_irrev, attractor_transition, attractor_basal_rev)
                                itp.set_phenotype_nodes(self.phenotype_nodes)
                                self.ITPs.append(itp)
                        
                        



