"""Usage.
Use a .bnet file to obtain the primes representation of a model.

Use the `find_attractors_asynchronous_update` function 
to obtain an `Attractor_Repertoire` object.

Use the `get_attractors` method of the `Attractor_Repertoire` object 
to obtain `Attractor` objects.

Information can then be extracted from these objects 
using their methods."""

import pystablemotifs as sm
from pyboolnet import prime_implicants
#pystablemotifs.__version__ == '3.0.6'
#pyboolnet.version.read_version() == '3.0.11'

import networkx as nx
#nx.__version__ == '3.1'

class Attractor_Repertore:
    # Stores a pyStableMotifs `AttractorRepertoire` object 
    # and provides methods for extracting the information 
    # that we are interested in from that object.

    #pystablemotifs의 attractorrepertore 객체를 담고,
    #그 객체에서 내가 원하는 정보를 추출해주는 methods로 구성됨.
    def __init__(self, att_rep_obj_in_pystablemotifs):
        self.att_rep_obj_in_pystablemotifs = att_rep_obj_in_pystablemotifs

    def get_attractors(self):
        return [Attractor(att_obj) for att_obj in self.att_rep_obj_in_pystablemotifs.attractors]

class Succession_Diagram:
    # Stores a pyStableMotifs `SuccessionDiagram`` object 
    # and provides methods for extracting the information 
    # that we are interested in from that object.
    def __init__(self, succession_diagram_obj):
        self.succession_diagram_obj = succession_diagram_obj
    
    def get_succession_diagram_graph(self):
        graph_edges = list(self.succession_diagram_obj.digraph.edges)
        return graph_edges
    
    def get_information_of_ith_node_in_succession_diagram(self, index):
        """Logically fixed nodes: the nodes and their states that are fixed at the current stage.
        
        Selectable stable motifs: the stable motifs that are active and selectable at the current stage. 
        Selecting one of these motifs creates a branch.

        Motif history: the stable motifs that have been selected to reach the current stage."""
        motif_reduction_obj = self.succession_diagram_obj.motif_reduction_dict[index]
        information_summary = {}
        information_summary["logically fixed nodes"] = motif_reduction_obj.logically_fixed_nodes
        information_summary["stable motifs selectible"] = motif_reduction_obj.stable_motifs
        information_summary["motif history"] = motif_reduction_obj.motif_history
        
        
        return information_summary
    
    def get_selected_stable_motif_for_edge(self, edge):
        """When moving from one node in the succession diagram to the next by selecting
        a stable motif, display the stable motif that was selected."""
        infor_of_ith = self.get_information_of_ith_node_in_succession_diagram(edge[0])
        infor_of_jth = self.get_information_of_ith_node_in_succession_diagram(edge[1])
        for stable_motif in infor_of_jth["motif history"]:
            if stable_motif not in infor_of_ith["motif history"]:
                return stable_motif




class Attractor:
    # Stores a pyStableMotifs `Attractor` object and provides methods
    # for extracting the information 
    # that we are interested in from that object.
    def __init__(self, att_obj_in_pystablemotifs):
        self.att_obj_in_pystablemotifs = att_obj_in_pystablemotifs
            
    def is_point_attractor(self):
        return self.att_obj_in_pystablemotifs.n_unfixed == 0
    
    def get_states_in_att(self):
        """Return the states contained in the attractor, 
        ignoring their ordering for now."""
        if self.is_point_attractor():
            state_dict_form = self.att_obj_in_pystablemotifs.attractor_dict.copy()
            return [state_dict_form]
        
        else:
            stg_of_states = self.att_obj_in_pystablemotifs.stg
            # Here, a state is represented simply as a string 
            # consisting of 0s and 1s.
            node_order = sorted(self.att_obj_in_pystablemotifs.reduced_primes)
            att_state_dict_forms = []
            for state in stg_of_states.nodes:
                state_dict_form = dict(zip(node_order, (int(i) for i in state)))
                # The states recorded in the STG nodes appear 
                # to include only the fluctuating nodes.
                # Therefore, the state information of fixed nodes is added separately.
                
                state_dict_form = {**state_dict_form, **self.att_obj_in_pystablemotifs.fixed_nodes}
                att_state_dict_forms.append(state_dict_form)
            
            return att_state_dict_forms
    
    def get_stg_of_att_states(self):
        """Return the result as an nx.DiGraph.

        Each node is assigned an index, 
        and that index corresponds to the position
        of the state in the list returned by get_states_in_att()."""
        if self.is_point_attractor():
            # In the case of a point attractor, 
            # return a single-node graph with a self-loop.
            nx_stg = nx.DiGraph()
            nx_stg.add_edge(0,0)
        
        else:
            state_index_map = {}
            for i, state in enumerate(self.att_obj_in_pystablemotifs.stg.nodes):
                state_index_map[state] = i
            
            nx_stg = nx.DiGraph()
            for source, target in self.att_obj_in_pystablemotifs.stg.edges:
                source_i = state_index_map[source]
                target_i = state_index_map[target]
                nx_stg.add_edge(source_i,target_i)
            
        return nx_stg


    





def find_attractors_asynchronous_update(primes, node_fixation:dict["node","state"]):
    """Given the PyBoolNet primes representation of a Boolean model,
    apply node fixation and then compute all attractors under asynchronous update.

    Among the many pieces of information available in the pyStableMotifs attractor objects,
    extract and use only those that are relevant for our purposes."""
    primes_with_node_fixation = prime_implicants.create_constants(primes, node_fixation, True)
    attractor_repertore = sm.AttractorRepertoire.from_primes(primes_with_node_fixation)
     
    return Attractor_Repertore(attractor_repertore)

def find_succession_diagram(primes):
    return Succession_Diagram(sm.succession.build_succession_diagram(primes))

if __name__ == "__main__":
    pass