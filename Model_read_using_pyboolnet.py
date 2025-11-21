import os

#pyboolnet.VERSION ==  3.0.11
from pyboolnet.file_exchange import bnet2primes, primes2bnet
import pyboolnet

def read_pyboolnet_file(address_file, encoding="utf8"):
    with open(address_file, 'r', encoding=encoding) as f:
        logic_text = f.read()
    logic_lines = [line for line in logic_text.split('\n') if line!='']
    logic_text = '\n'.join(logic_lines)
    primes = bnet2primes(logic_text)
    dynamics_pyboolnet = Dynamics_pyBoolnet(primes)

    return dynamics_pyboolnet

class Dynamics_pyBoolnet:
    def __init__(self, primes):
        self.primes = primes
        self.nodes = sorted(list(self.primes.keys()))

    def get_primes(self, reduction=True, perturbations={}):
        """Returns the primes.

        If reduction == True, after computing the LDOI, 
        the logic associated with nodes whose values are already fixed 
        is excluded from the returned result."""
        primes_perturbed = pyboolnet.prime_implicants.percolate(self.primes,
                                                                add_constants=perturbations, 
                                                                remove_constants=reduction,
                                                                copy=True)
            #side effect에 의해 primes_perturbed 가 바뀌게 된다.
        
        return primes_perturbed

    def get_node_names(self):
        """return all nodes in the network model as a list form"""
        return self.nodes
    
    def get_source_node_names(self):
        """Returns the source nodes as a list.
        
        A node is considered a source node if its primes value is of the form:
        [[{'n1': 0}], [{'n1': 1}]] for node 'n1'.
        
        Needs further review to ensure correct behavior under perturbed models."""
        source_node_names = []
        for nodename, prime in self.primes.items():
            if prime[0] == [{nodename:0}]:
                if prime[1] == [{nodename:1}]:
                    source_node_names.append(nodename)
        
        return source_node_names
    
    def print_cytoscape_file(self):
        """Returns a text-formatted file used for visualizing the network in Cytoscape.
        
        Functionality to append state values for each node according to its attributes (keys) will be added later.
        """
        nx_net = pyboolnet.interaction_graphs.primes2igraph(self.primes)
        edges_txt = "source\ttarget\tmodality\n"
        for edge in nx_net.edges:
            if nx_net.edges[edge]['sign'] == {1}:
                modality = '+'
            elif nx_net.edges[edge]['sign'] == {-1}:
                modality = '-'
            else:
                raise ValueError("unknown case of interaction sign information")
            edges_txt += "{}\t{}\t{}\n".format(edge[0], edge[1], modality)
        
        return edges_txt
    
    def get_edges_with_modalities(self, node_cut=[], ignore_self_loop_on_source_nodes=True):
        """Extracts edge information using the primes data and returns it as a list.
        Each edge is represented as a tuple (source, modality, target).
        
        Edges involving any node included in `node_cut` are excluded from the result.
        """
        links = []
        for target_node, prime_implicant in self.primes.items():
            if target_node in node_cut:
                continue
            if ignore_self_loop_on_source_nodes and (target_node in self.get_source_node_names()):
                continue
            prime_implicant_off = prime_implicant[0]
            for canalizing_condition in prime_implicant_off:
                for source_node, state in canalizing_condition.items():
                    if source_node in node_cut:
                        continue
                    
                    if state == 0:
                        modality = "+"
                    else:
                        modality = '-'
                    edge = (source_node, modality, target_node)
                    if edge not in links:
                        links.append(edge)
        
        return links
        

if __name__ == "__main__":
    add_test_model = os.path.join(r"D:\new canalizing kernel\우정 tumorigenesis 모델", "Fumia_logic_model.bnet")
    dynamics_test_model = read_pyboolnet_file(add_test_model)