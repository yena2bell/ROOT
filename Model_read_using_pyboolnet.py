import os

#pyboolnet.VERSION ==  3.0.11 사용.
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
        """primes 를 return한다.
        reduction==True이면 이미 값이 고정되어 있는 것들에 대해, LDOI 계산 후,
        그것들을 제외한 logic 을 return"""
        primes_perturbed = pyboolnet.prime_implicants.percolate(self.primes,
                                                                add_constants=perturbations, 
                                                                remove_constants=reduction,
                                                                copy=True)
            #side effect에 의해 primes_perturbed 가 바뀌게 된다.
        
        return primes_perturbed

    def get_node_names(self):
        """network 가 가지는 모든 nodes의 name 을 list에 담아서
        return한다."""
        return self.nodes
    
    def get_source_node_names(self):
        """primes 가 'n1': [[{'n1': 0}], [{'n1': 1}]] 형태인 것을 골라서 
        list 형태로 return한다. source nodes를 골라내는 용도이지만 perturbed model에서도
        잘 작동할 지 확신할 수 없음."""
        source_node_names = []
        for nodename, prime in self.primes.items():
            if prime[0] == [{nodename:0}]:
                if prime[1] == [{nodename:1}]:
                    source_node_names.append(nodename)
        
        return source_node_names
    
    def print_cytoscape_file(self):
        """cytoscape에서 네트워크를 보여주는 용도의 파일을 text로 return한다.
        
        각 속성(key)에 따른 node의 상태값을 추가로 적어주는 기능을 추가하기."""
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
        """primes 정보를 활용하여 links를 추출.
        각 link는 (source, modality, target)으로 구성된 tuple이고,
        이러한 tuple로 구성된 list를 return
        
        node_cut에 들어있는 node의 경우 그 node가 포함된 links는 제외하고 return"""
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