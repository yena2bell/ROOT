# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 14:41:04 2021

@author: jwKim
"""
import numpy as np

def get_nodes_in_links(links:"[(node1,... ,node2), (node1,..., node3)...]"):
    """links에 포함된 nodes를 찾아서 return 한다."""
    nodes_in_links = set()
    for link in links:
        nodes_in_links.add(link[0])
        nodes_in_links.add(link[-1])
    return list(nodes_in_links)
    

def SCC_decomposition(ltLinks_data):
    """
    ltLinks_data = [(node1, node2), (node1, node3)...]
    (node1, node2) means node1 interacte to node2 i.e. node1 -> node2
    """
    lRemained_links = list(ltLinks_data)
    #print("the number of links to analyze is", len(lRemained_links))
    #copy the list data to conserve original data

    lRemained_nodes = get_nodes_in_links(lRemained_links)
    #print("the number of nodes to analyze is", len(lRemained_nodes))
    llSCC = []

    while(lRemained_nodes):
        llSCC += _find_SCC_under_startnode(lRemained_nodes.pop(0), lRemained_nodes, lRemained_links)
        """
        choose one node and make it start node.
        find SCC containing that start node and SCC whose hierarchy is lower than SCC containing start node
        repeat until find all SCCs
        """

    return llSCC

def get_SCC_containing_the_node(s_node, l_t_links_data):
    """
    ltLinks_data = [(node1, node2), (node1, node3)...]
    (node1, node2) means node1 interacte to node2 i.e. node1 -> node2
    s_node 가 포함된 최소의 SCC를 list 형태로 추출한다.
    """
    l_remained_links = l_t_links_data.copy()
    l_remained_nodes = get_nodes_in_links(l_remained_links)
    
    l_remained_nodes.remove(s_node)
    #print("the number of nodes to analyze is", len(lRemained_nodes))
    llSCC = _find_SCC_under_startnode(s_node, l_remained_nodes, l_remained_links)
    return llSCC[node_position_finding(llSCC, s_node)]
    


def _find_SCC_under_startnode(flownode, lRemained_nodes, ltRemained_links):
    """
    sub function of 'SCC_decomposition'
    nodes data = [node name1, node name2, ..... node name k]
    interaction data = [(node1, node2), (node1, node3)...]
    (node1, node2) means node1 interacte to node2 i.e. node1 -> node2
    remained nodes set don't contain start node
    """

    lFlow_of_node = [flownode]
    llSet_of_SCC = []
    ltCycle_positions = []
    
    while lFlow_of_node:
        ltRemained_links = list(filter(lambda x: x != (None,None),ltRemained_links))
        
        for tNum_tNode1tonode2 in enumerate(ltRemained_links):
            if tNum_tNode1tonode2[1][0] == flownode:
                nextnode = tNum_tNode1tonode2[1][1]
                ltRemained_links.pop(tNum_tNode1tonode2[0])
                ltRemained_links.insert(tNum_tNode1tonode2[0],(None,None))
                if nextnode in lFlow_of_node:
                    tCycle = (lFlow_of_node.index(nextnode),len(lFlow_of_node)-1)
                    ltCycle_positions = _evaluate_SCC_inclusion(ltCycle_positions, tCycle) 
                    break 
                elif nextnode in lRemained_nodes:
                    lFlow_of_node.append(nextnode)
                    lRemained_nodes.pop(lRemained_nodes.index(nextnode))
                    break
                else:
                    continue 
                
        else:
            iPosi_of_flownode = lFlow_of_node.index(flownode)
            if not(ltCycle_positions):
                llSet_of_SCC.append([lFlow_of_node.pop(-1)]) 
            elif iPosi_of_flownode >ltCycle_positions[-1][1]:
                llSet_of_SCC.append([lFlow_of_node.pop(-1)])
            elif iPosi_of_flownode == ltCycle_positions[-1][0]:      
                tSCC = ltCycle_positions.pop(-1)
                llSet_of_SCC.append(lFlow_of_node[tSCC[0]:tSCC[1]+1])
                lFlow_of_node = lFlow_of_node[:tSCC[0]]
            if lFlow_of_node:
                flownode = lFlow_of_node[iPosi_of_flownode-1]
                continue
            else:
                continue
        
        flownode = lFlow_of_node[-1] 
       
    return llSet_of_SCC


def _evaluate_SCC_inclusion(ltCycles,tNewCycle):
    """
    sub function of '_find_SCC_under_startnode'
    ltCycles = [(3,8},(11, 20) ....] if ltCycle is [(a,b), (c,d)...] then a,b,c,d satisfy a<b<c<d 
    tNewCycle = (50, 78) this means that node 50, node51.... node 78 is SCC
    tCycle is the record of feedback in the flow. if tCycle==(a,b) then a<b
    """

    for tCycle in ltCycles:
        if tCycle[1] < tNewCycle[0]:
            continue #tCycle has no common nodes with tNewCycle. so go to next cycle
        elif tCycle[0] <= tNewCycle[0]:
            tNewCycle = (tCycle[0],tNewCycle[1])
            iPosi_of_tCycle = ltCycles.index(tCycle) 
            break
        else: #tCycle[0]>tNewCycle
            iPosi_of_tCycle = ltCycles.index(tCycle)
            break
    else:
        ltCycles.append(tNewCycle)
        return ltCycles

    ltCycles = ltCycles[0:iPosi_of_tCycle]
    ltCycles.append(tNewCycle)
    return ltCycles


def net_of_SCCs(llSCC, ltLinks_data):
    """
    llSCC = [[node1,node2,node3... nodes in the SCC1],[node6,node7.... nodes in the SCC2],....]
    ltLinks_data = [(node1, node2), (node1, node3)...]
    give each SCC the number as llSCC.index
    calculate interactions between SCCs.
    return the list [(SCC1's number,SCC2's number),...] this means that link connecting SCC1 -> SCC2 exists
    """
    l_remained_links = list(ltLinks_data)
    lt_SCClinks = []
    
    for i_SCC, lSCC in enumerate(llSCC):
        for i_link in range(len(l_remained_links)-1,-1,-1):
            if l_remained_links[i_link][0] in lSCC:
                if l_remained_links[i_link][1] in lSCC:
                    l_remained_links.pop(i_link)
                else:
                    t_tmplink = l_remained_links.pop(i_link)
                    t_SCClink = (i_SCC,node_position_finding(llSCC, t_tmplink[1]))
                    lt_SCClinks.append(t_SCClink)
            elif l_remained_links[i_link][1] in lSCC:
                t_tmplink = l_remained_links.pop(i_link)
                t_SCClink = (node_position_finding(llSCC, t_tmplink[0]), i_SCC)
                lt_SCClinks.append(t_SCClink)

    lt_SCClinks =list(set(lt_SCClinks))

    return lt_SCClinks


def node_position_finding(llSCC, s_node):
    """
    sub function of 'net_of_SCCs'
    llSCC = [[node1,node2,node3... nodes in the SCC1],[node6,node7.... nodes in the SCC2],....]
    s_node is node name
    """
    for i, lSCC in enumerate(llSCC):
        if s_node in lSCC:
            return i
    else:
        raise(ValueError("error in node_position_finding function."))


def highest_SCCs_finding(llSCC, lt_SCClinks, b_onenodeSCC = False):
    """
    llSCC = [[node1,node2,node3... nodes in the SCC1],[node6,node7.... nodes in the SCC2],....]
    lt_SCClinks = [(0,2), (3,4), ...] (a,b) means that SCC which is llSCC[a] has direct links to SCC which is llSCC[b]
    lt_SCClinks can be result of net_of_SCCs function
    this function returns SCC number which has no upper SCCs
    if b_onenodeSCC == False, then highest SCCs should have more than 2 nodes
    """
    l_ihighest_hierarchy = []
    for i in range(len(llSCC)):
        for t_link in lt_SCClinks:
            if t_link[1] == i:
                break
        else:
            l_ihighest_hierarchy.append(i)

    if b_onenodeSCC:
        return(l_ihighest_hierarchy)
    else:#highest nodes should have more than or equal to 2 nodes
        l_candidate = []
        for i in range(len(l_ihighest_hierarchy)-1,-1,-1):
            if len(llSCC[l_ihighest_hierarchy[i]]) == 1:
                itmp = l_ihighest_hierarchy.pop(i)
                for j in range(len(lt_SCClinks)-1,-1,-1):
                    if lt_SCClinks[j][0] == itmp:
                        l_candidate.append(lt_SCClinks[j][1])
                        lt_SCClinks.pop(j)
        l_candidate = list(set(l_candidate))

        while l_candidate:
            for i in l_candidate:
                for t_link in lt_SCClinks:
                    if t_link[1] == i:
                        break
                else:
                    l_ihighest_hierarchy.append(i)

            l_candidate = []
            for i in range(len(l_ihighest_hierarchy)-1,-1,-1):
                if len(llSCC[l_ihighest_hierarchy[i]]) == 1:
                    itmp = l_ihighest_hierarchy.pop(i)
                    for j in range(len(lt_SCClinks)-1,-1,-1):
                        if lt_SCClinks[j][0] == itmp:
                            l_candidate.append(lt_SCClinks[j][1])
                            lt_SCClinks.pop(j)
            l_candidate = list(set(l_candidate))

    return [llSCC[i] for i in l_ihighest_hierarchy]

def decompose_to_SCC_from_matrix(matrix):
    l_l_SCC = []
    l_node_flow = SCC_algorithm_Kosaraju_stack_calculation(matrix)
    
    array_visited = np.zeros(len(matrix),dtype=bool)
    
    matrix_transpose = np.transpose(matrix.copy())
    
    while l_node_flow:
        i_index = l_node_flow.pop()
        l_SCC = []
        if not array_visited[i_index]:
            l_l_SCC.append(_DFSUtil(i_index, array_visited, 
                                    matrix_transpose, l_SCC))
    
    return l_l_SCC

def is_SCC_unsigned_graph(graph):
    return is_SCC(graph.show_unsigned_graph_matrix_form())
            
def is_SCC(matrix):
    l_node_flow = SCC_algorithm_Kosaraju_stack_calculation(matrix)
    
    array_visited = np.zeros(len(matrix),dtype=bool)
    
    matrix_transpose = np.transpose(matrix.copy())
    while l_node_flow:
        i_index = l_node_flow.pop()
        l_SCC=[]
        if not array_visited[i_index]:
            l_SCC = _DFSUtil(i_index, array_visited, 
                             matrix_transpose, l_SCC)
            if len(l_SCC) < len(matrix):
                return False
            else:
                return True

def SCC_algorithm_Kosaraju_stack_calculation(matrix):
    l_node_flow = []
    matrix_directed = matrix.copy()
    array_visited = np.zeros(len(matrix),dtype=bool)
    
    for i_index in range(len(matrix)):
        if not array_visited[i_index]:
            _fill_order(i_index, array_visited, 
                        l_node_flow, matrix_directed)
    
    return l_node_flow#at least one node of the upper hierarchy SCC is located latter in the l_node_flow than the lower hierarchy SCC
            
def _DFSUtil(i_index, array_visited, matrix_transpose, l_SCC):
    array_visited[i_index] = True
    l_SCC.append(i_index)
    for i_downstream in np.nonzero(matrix_transpose[:, i_index])[0]:
        if not array_visited[i_downstream]:
            l_SCC = _DFSUtil(i_downstream, array_visited, 
                                  matrix_transpose, l_SCC)
    return l_SCC
            
            
def _fill_order(i_index, array_visited, l_node_flow, matrix_directed):
    array_visited[i_index] = True
    for i_downstream in np.nonzero(matrix_directed[:,i_index])[0]:
        if not array_visited[i_downstream]:
            _fill_order(i_downstream, array_visited, 
                        l_node_flow, matrix_directed)
    l_node_flow.append(i_index)