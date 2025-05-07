import bisect
import time
import random
import math

from Network_state_and_attractor import Network_state, Attractor 

class Attractor_landscape_for_specific_IC:
    """calculate attractor landscape for a specific input configuraiton,
    and save the results in the form of attractor index and basin states."""
    def __init__(self, dynamics_Boolean_net:"Dynamics_pyBoolnet", IC:dict, fixed_node_state_map:dict={}):
        self.dynamics_Boolean_net = dynamics_Boolean_net
        self.IC = IC
        #input configuration
        self.fixed_node_state_map = fixed_node_state_map
        #node fixation by mutation

        self.fixed_node_state_map_and_IC = self.fixed_node_state_map.copy()
        self.fixed_node_state_map_and_IC.update(self.IC)

        self.attractor_index_map = {} #{0:attractor_0, 1:attractor_1,...}
        self.attindex_basinstates_map = {}
        #특정 attractor에 대해 basin으로 밝혀진 network state values의 list를 그 attractor의 index와 mapping
    
    def __repr__(self):
        return "Attractor landscape for IC {}".format(self.IC)
    
    @property
    def attindex_basinratio_map(self):
        """self.attindex_basinstates_map에서의 basin의 크기의 비율을 계산하여 return
        각각의 attractor index가 의미하는 attractor는 self.attractor_index_map에서 찾을 수 있다."""
        attindex_basinsize_map = {index:len(basin) for index, basin in self.attindex_basinstates_map.items()}
        total_basin_size = sum(attindex_basinsize_map.values())
        attindex_basinratio_map = {index:basin_size/total_basin_size for index, basin_size in attindex_basinsize_map.items()}
        return attindex_basinratio_map
    
    @property
    def num_of_states_searched(self):
        attindex_basinsize_map = {index:len(basin) for index, basin in self.attindex_basinstates_map.items()}
        return sum(attindex_basinsize_map.values())
    
    def get_num_of_not_fixed_nodes(self):
        """return 
        num of all nodes - num of (union(fixed ndoes, input nodes))"""
        return len(self.dynamics_Boolean_net.get_node_names()) - len(set(self.IC).union(set(self.fixed_node_state_map)))

    def get_attractor_index_from_attractor(self, attractor):
        for index, attractor_in_map in self.attractor_index_map.items():
            if attractor_in_map == attractor:
                return index
        else:
            # when the new attractor is found during iATG construction,
            # the new attractor is not contained in self.attractor_index_map.
            # then, add the new attractor in the self.attractor_index_map,
            # and basin states composed of network states in the new attractor 
            # is added to self.attindex_basinstates_map
            index_of_new_attractor = len(self.attractor_index_map)
            self.attractor_index_map[index_of_new_attractor] = attractor
            basin_states_of_new_attractor = []
            _put_int_values_in_list(attractor.attractor_states, basin_states_of_new_attractor)
            self.attindex_basinstates_map[index_of_new_attractor] = basin_states_of_new_attractor

            return index_of_new_attractor


            
    def calculate_attractor_basinratios_from_random_initial_states(self, waiting_num=10000, 
                                                   difference_threshold=0.0001,
                                                   attractors_exist=[],
                                                   verbose=False):
        num_of_not_fixed_nodes = len(self.dynamics_Boolean_net.get_node_names()) - len(self.fixed_node_state_map_and_IC)
        num_of_all_states = pow(2, num_of_not_fixed_nodes)
        time_start = time.time()

        attractors_exist_not_yet_found = attractors_exist.copy()

        list_of_network_state_values_detected = []
        # Stores states that have been explored at least once.  
        # Managed using bisect.


        i_no_more_att_discovery_count = 0
        attractor_basinratio_map_prev = {}
        num_of_attractors_found = len(self.attractor_index_map)
        # if i_no_more_att_discovery_count >= waiting_num, 
        # and attractors in attractors_exist are all found (i.e. attractors_exist_not_yet_found is empty),
        # then check the convergence of attractor basin ratios.
        # if converged, break the loop.
        while True:
            if i_no_more_att_discovery_count >= waiting_num:
                # first codition to stop the while loop
                if verbose:
                    num_of_searched_initial_states = len(list_of_network_state_values_detected)
                    percentage_of_searched_initial_states = 100*num_of_searched_initial_states/num_of_all_states

                    print("initial states searching proceed {:.2f}% ({}/{})".format(percentage_of_searched_initial_states, 
                                                                                num_of_searched_initial_states, 
                                                                                num_of_all_states))
                    time_elapsed_until_now = time.time() - time_start
                    hour_elapsed = time_elapsed_until_now//3600
                    min_elapsed = time_elapsed_until_now//60
                    sec_elapsed = time_elapsed_until_now%60

                    print("{}h {}min {:.2f}sec have elapsed".format(hour_elapsed, min_elapsed, sec_elapsed))
        

                attractors_exist_not_yet_found = [attractor for attractor in attractors_exist_not_yet_found if attractor not in self.attractor_index_map.values()]
                
                if not attractors_exist_not_yet_found:
                    # second condition to stop the while loop

                    if self._attractor_landscape_is_converged(attractor_basinratio_map_prev, difference_threshold):
                        # third condition to stop the while loop
                        break
                    else:
                        if verbose:
                            print("Attractors are found but not converged yet.")
                        i_no_more_att_discovery_count = 0
                        
                else:
                    if verbose:
                        print("Attractors not found yet: ", attractors_exist_not_yet_found)
                    i_no_more_att_discovery_count = 0
                
                attractor_basinratio_map_prev = self.attindex_basinratio_map
                # compare this attractor_basinratio_map_prev with the current one
                    
            else:
                network_state = get_random_state(self.dynamics_Boolean_net, self.fixed_node_state_map_and_IC, list_of_network_state_values_detected)
                #network state not in list_of_network_state_values_detected

                self._analyze_trajectory_from_network_state(network_state, list_of_network_state_values_detected)

                if len(self.attractor_index_map) > num_of_attractors_found:
                    #new attractor is found
                    i_no_more_att_discovery_count = 0
                    num_of_attractors_found = len(self.attractor_index_map)
                else:
                    i_no_more_att_discovery_count += 1


        print("time for {}: ".format(str(self)), time.time()-time_start)
    
    def _analyze_trajectory_from_network_state(self, network_state,list_of_network_state_values_detected):
        """
        Finds the trajectory starting from network_state.
        Stops when either a new attractor is formed or 
        it reaches one of the states in list_of_network_state_values_detected.

        - If a new attractor is formed, the attractor and its trajectory are stored in self.attractor_index_map and self.attindex_basinstates_map.
        - If it reaches one of the detected states, the trajectory is added to the basin of the corresponding existing attractor.
        """
        trajectory, network_state_value_next = self._calculate_trajectory_and_check(network_state, list_of_network_state_values_detected)

        if network_state_value_next in trajectory:
            #this trajectory makes a new attractor
            attractor_new = self._make_new_attractor_from_trajectory(trajectory, network_state_value_next)
            index_of_attractor = len(self.attractor_index_map)
            self.attractor_index_map[index_of_attractor] = attractor_new

            network_state_values_in_basin = []
            _put_int_values_in_list(trajectory, network_state_values_in_basin)
            self.attindex_basinstates_map[index_of_attractor] = network_state_values_in_basin
        
        else:
            # this trajectory reaches one of the existing attractors
            for basin in self.attindex_basinstates_map.values():
                if network_state_value_next in basin:
                    _put_int_values_in_list(trajectory, basin)
                    break
        
        _put_int_values_in_list(trajectory, list_of_network_state_values_detected)
        #add trajectory to detected list



    def _attractor_landscape_is_converged(self, attractor_basinratio_map_prev, difference_threshold):
        """compare the previous attractor basinratio map to the attractor basinratio map now.
        if the difference is less than the difference_threshold, return True."""
        attractor_basinratio_map_now = self.attindex_basinratio_map
        if len(attractor_basinratio_map_now) != len(attractor_basinratio_map_prev):
            return False
        
        difference_measure = calculate_JSD_on_discrete_distribution(attractor_basinratio_map_prev, attractor_basinratio_map_now)
        if difference_threshold > difference_measure:
            return True
        else:
            return False


    def calculate_attractor_basinratios_from_all_initial_states(self):
        """주어진 IC에 대해, 모든 initial state로부터 attractor basin ratio를 계산한다.
        
        검사하기. attindex_basinstates_map의 basins 사이에 겹치는 것이 없는가?
        개수는 일치하는가?
        fixed node state가 정해진 경우에도?"""
        time_start = time.time()
        num_of_all_initials = pow(2, self.get_num_of_not_fixed_nodes())
        percentage_to_report = 0
        reporting_interval = 0.2

        list_of_network_state_values_detected = []
        # Stores states that have been explored at least once.  
        # Managed using bisect.  
       

        for network_state in get_all_state(self.dynamics_Boolean_net, self.fixed_node_state_map_and_IC, list_of_network_state_values_detected):

            if len(list_of_network_state_values_detected) >= num_of_all_initials*percentage_to_report:
                percentage_checked = len(list_of_network_state_values_detected)/num_of_all_initials*100
                print("{}% of all initial states are detected.".format(percentage_checked))
                print("used time: ", time.time()-time_start)
                percentage_to_report += reporting_interval

            
            #잘 돌아가면 이 부분 지우기
            # trajectory, network_state_value_next = self._calculate_trajectory_and_check(network_state, list_of_network_state_values_detected)

            # if network_state_value_next in trajectory:
            #     #새로운 attractor가 생성됨
            #     attractor_new = self._make_new_attractor_from_trajectory(trajectory, network_state_value_next)
            #     index_of_attractor = len(self.attractor_index_map)
            #     self.attractor_index_map[index_of_attractor] = attractor_new
                
            #     network_state_values_in_basin = []
            #     _put_int_values_in_list(trajectory, network_state_values_in_basin)
            #     self.attindex_basinstates_map[index_of_attractor] = network_state_values_in_basin

            # else:
            #     #기존 attractor 중 하나의 basin에 network_state_value_next이 포함됨.
            #     for basin in self.attindex_basinstates_map.values():
            #         if network_state_value_next in basin:
            #             _put_int_values_in_list(trajectory, basin)
            #             break
            
            # _put_int_values_in_list(trajectory, list_of_network_state_values_detected)
            # #add trajectory to detected list
            self._analyze_trajectory_from_network_state(network_state, list_of_network_state_values_detected)
        
        print("time for {}: ".format(str(self)), time.time()-time_start)
        
    
    def _calculate_trajectory_and_check(self, 
                                        network_state_value_start:Network_state,
                                        list_of_network_state_values_detected:[]):
        """network_state_value_start로 시작하는 trajectory를 구한 뒤, 
        그것이 새로운 attractor를 만들거나, 
        list_of_network_state_values_detected 중 하나에 닿게 되면 멈추고
        그 결과를 return"""
        trajectory = [network_state_value_start]
        #trajectory[i]는 시작 state로부터 i step 후의 state
        while True:
            network_state_value_next = model_state_synchronous_update_using_pyboolnet(self.dynamics_Boolean_net, 
                                                                        trajectory[-1], 
                                                                        self.fixed_node_state_map_and_IC)
            
            if network_state_value_next in trajectory:
                #새로운 attractor가 생성됨
                return trajectory, network_state_value_next
            elif network_state_value_next in list_of_network_state_values_detected:
                #이미 찾은 state에 도달함
                return trajectory, network_state_value_next
            else:
                trajectory.append(network_state_value_next)
    
    def converge_network_state_value_to_attractor(self, network_state_value:Network_state):
        """Return the attractor reached as a result of network state transitions 
        from the given network state value, 
        under the input configuration and fixed node state map of this attractor landscape.
        
        주어진 network state value로부터, 
        이 attractor landscape의 input configuration과 fixed node state map 상태에서 
        network state transition의 결과 도달하는 attractor 를 return한다. """
        trajectory, network_state_value_next = self._calculate_trajectory_and_check(network_state_value, [])
        return self._make_new_attractor_from_trajectory(trajectory, network_state_value_next)

    
    def _make_new_attractor_from_trajectory(self, trajectory, network_state_value_next):
        """network_state_value_next는 trajectory 안에 이미 존재하는 network state value이다.
        따라서 trajectory의 network_state_value_next가 있는 곳부터 마지막까지의 network state values가
        atttractor states가 된다."""
        index = trajectory.index(network_state_value_next)
        attractor_states = trajectory[index:]
        attractor = Attractor(self.dynamics_Boolean_net)
        attractor.put_attractor_states_in_synchro_using_network_state_forms(attractor_states, self.fixed_node_state_map)

        return attractor

            


#
# attractor landscape를 구하기 위해 다수의 intial states를 찾는 것과 관련된 함수들 모음
#
def get_random_state(dynamics_pyboolnet, fixed_node_state_map={}, states_except=[]):
    """주어진 dynamics_pyboolnet에서 perturbation을 고려한 상태로 임의의 state를 뽑는다.
    
    단 states_except에 있는 것들은 제외하고.
    states_except에는 State_in_pyboolnet 객체나 int 값들이 들어간다.
    이미 찾은 state는 빼고 뽑는 용도로 사용함.
    !!!!!!!!!!!!!!이 때, states_except는 bisect를 통해서 관리되도록 할 것.!!!!!!!!!!!!!!!!!!!!1"""
    num_of_nodes = len(dynamics_pyboolnet.get_node_names())
    num_of_possible_states = pow(2, num_of_nodes)
    while True:
        random_int = random.randint(0, num_of_possible_states-1)
        if fixed_node_state_map:
            dict_form = Network_state._convert_int_form_to_dict_form(dynamics_pyboolnet, random_int)
            dict_form = Network_state.apply_perturbation_to_dict_form_state(dict_form, fixed_node_state_map)
            random_int = Network_state._convert_dict_form_to_int_form(dynamics_pyboolnet, dict_form)

        position_to_put = _check_existence_of_int_in_list(random_int, states_except)
        #random_int가 states_except에 이미 있으면 None을 return
        #없으면 넣어야 하는 위치를 return
        if position_to_put is not None:
            state_obj = Network_state(dynamics_pyboolnet)
            state_obj.put_state_int_form(random_int)

            return state_obj

def get_all_state(dynamics_pyboolnet, fixed_node_state_map={}, states_except=[]):
    """주어진 dynamics_pyboolnet에서 perturbation을 고려한 상태로 모든 state를 순차적으로 뽑는다.
    iter 객체로 사용한다.
    
    단 states_except에 있는 것들은 제외하고.
    states_except에는 State_in_pyboolnet 객체나 int 값들이 들어간다.
    이미 찾은 state는 빼고 뽑는 용도로 사용함.

    !!!!!!!!!!!!!!이 때, states_except는 bisect를 통해서 관리되도록 할 것.!!!!!!!!!!!!!!!!!!!!1"""
    num_of_nodes = len(dynamics_pyboolnet.get_node_names())
    num_of_possible_states = pow(2, num_of_nodes)
    if fixed_node_state_map:
        for int_form in range(num_of_possible_states):
            dict_form = Network_state._convert_int_form_to_dict_form(dynamics_pyboolnet, int_form)
            dict_form = Network_state.apply_perturbation_to_dict_form_state(dict_form, fixed_node_state_map)
            int_form_perturbed = Network_state._convert_dict_form_to_int_form(dynamics_pyboolnet, dict_form)
        
            position_to_put = _check_existence_of_int_in_list(int_form_perturbed, states_except)
            #int_form_perturbed가 states_except에 이미 있으면 None을 return
            #없으면 넣어야 하는 위치를 return
            if position_to_put is not None:
                network_state_value = Network_state(dynamics_pyboolnet)
                network_state_value.put_state_int_form(int_form_perturbed)
                yield network_state_value
    
    else:
        #굳이 fixed node state 판정을 for 문 밖으로 뺀 것은
        # fixed ndoes가 없을 경우, 이렇게 하는 것이 조금이라도 더 빨라져서
        for int_form in range(num_of_possible_states):
            position_to_put = _check_existence_of_int_in_list(int_form, states_except)
            #int_form가 states_except에 이미 있으면 None을 return
            #없으면 넣어야 하는 위치를 return
            if position_to_put is not None:
                network_state_value = Network_state(dynamics_pyboolnet)
                network_state_value.put_state_int_form(int_form)
                yield network_state_value


#
# dynamics에 따른 node, state update 관련 함수들 모음.
# attractor까지 converging하는 함수도 포함.
#
def model_state_synchronous_update_using_pyboolnet(dynamics_pyboolnet, 
                                                   network_state_value_now:Network_state, 
                                                   node_fixed_state_map={}):
    """현재의 model state에서 
    먼저 node_fixed_state_map을 반영한 뒤,
    synchronous update를 수행했을 때 나오는 next network state value를 return
    
    따라서 not fixed nodes가 update 될 때, regulator 중 fixed node가 있으면 그 fixed value로 영향을 받음."""
    network_state_value_dict_form = network_state_value_now.get_state_dict_form()
    network_state_value_dict_form = network_state_value_now.apply_perturbation_to_dict_form_state(network_state_value_dict_form, node_fixed_state_map)
    #node_fixed_state_map 덮어 씌우기
    
    network_state_value_dict_form_next = {}
    for node in network_state_value_dict_form:
        if node not in node_fixed_state_map:
            network_state_value_dict_form_next[node] = node_state_update_using_pyboolnet(dynamics_pyboolnet, network_state_value_dict_form, node)
        else:
            network_state_value_dict_form_next[node] = node_fixed_state_map[node]
    
    network_state_value_next = Network_state(dynamics_pyboolnet)
    network_state_value_next.put_state_dict_form(network_state_value_dict_form_next)
    
    return network_state_value_next

def node_state_update_using_pyboolnet(dynamics_pyboolnet, network_state_value_dict_form, nodename_to_update):
    """network_state_value_dict_form에서 특정 node의 next state를 구한다.
    network_state_value_dict_form는 그 node의 regulator nodes만 다 포함하면 충분하다.
    """
    prime_of_node = dynamics_pyboolnet.primes[nodename_to_update]
    prime_for_0 = prime_of_node[0]
    prime_for_1 = prime_of_node[1]
    if len(prime_for_0) >= len(prime_for_1):
        state_to_check = 1
        primes_to_check = prime_for_1
    else:
        state_to_check = 0
        primes_to_check = prime_for_0
    #0에 대한 primes 나 1에 대한 primes 중 하나만 보면 되기 때문에
    #더 볼 것이 적은 primes 를 골라서 진행한다.

    for prime in primes_to_check:
        for node, state in prime.items():
            if int(network_state_value_dict_form[node]) != state:
                break
        else:
            return state_to_check
        #state_of_model 이 prime 을 온전히 포함할 경우,
        #이 node 는 state_to_check 로 update 된다.
    else:
        #주어진 primes_to_check 에서 전부 조건을 만족하지 않으면
        #반대 값으로 update 된다.
        return 1-state_to_check


###############################################################################################
### bisect related functions
###############################################################################################

def _check_existence_of_int_in_list(int_value, list_of_int):
    """list_of_int 는 오름차순으로 정렬된 int 값들의 list
    
    int_value가 list_of_int 안에 없을 경우, 어디에 삽입하면 될지 좌표를 return
    이미 있을 경우, None을 return
    
    주의할 점은 i==len(list_of_int)가 나올 경우, append를 통해 추가해 줘야 한다는 것."""
    i = bisect.bisect_left(list_of_int, int_value)
    if i != len(list_of_int) and list_of_int[i] == int_value:
        return None
    else:
        return i

def _put_int_value_in_list(value, list_to_put, index):
    """index == line(list_to_put)이면 append로 넣어준다."""
    if index == len(list_to_put):
        list_to_put.append(value)
    else:
        list_to_put.insert(index, value)

def _put_int_values_in_list(values:iter, list_to_put):
    """values의 값을 bisect를 사용해서 이미 오름차순으로 정리된 list_to_put에 넣어준다."""
    for value in values:
        bisect.insort_left(list_to_put, value)

###############################################################################################
### Distribution comparison measures related functions
###############################################################################################


def calculate_KLD_on_discrete_distribution(P,Q, base_of_log=math.e, epsilon=0):
    """KLD는 P대신 Q를 사용할 경우 entropy 변화를 의미
    
    어떤 P와 Q의 key 값 중 하나의 분포에만 존재하는 경우가 있으면,
    다른 분포에서 그 key 값의 확률값은 0으로 취급
    
    Q[key] = 0인데 P[key]!=0인 경우가 있으면 정의가 안된다.
    
    이를 방지하기 위해 epsilon에 매우 작은 값을 넣어서 Q[key]==0, P[key]!=0인 경우도 고려하는
    Smoothed KLD라는 기법이 있다고 함."""
    KLD = 0
    keys = set(list(P) + list(Q))
    if epsilon:
        #Smoothed KLD
        for key in keys:
            if P.get(key,0) == 0:
                pass
            else:
                KLD += P.get(key,0)*math.log(P.get(key,0)/(Q[key]+epsilon), base_of_log)
    
    else:
        #not smooted KLD
        for key in keys:
            if P.get(key,0) == 0:
                #P[key]*log(P[key]/Q[key]) = 0*log(0/Q[key]) = 0
                pass
            elif Q.get(key,0) == 0:
                raise ValueError("Q[{}] is zero although P[{}] is not zero".format(key))
            else:
                KLD += P.get(key,0)*math.log(P.get(key,0)/Q[key], base_of_log)
        
    return KLD

def calculate_JSD_on_discrete_distribution(P,Q, base_of_log=math.e):
    """Jensen-Shannon Divergence을 계산.
    JSD(P∥Q)=0.5(KLD(P∥M)+KLD(Q∥M))
    M=0.5(P+Q)로 계산된다.
    
    symmetric하다.
    Q[key] = 0인데 P[key]!=0인 경우가 있어도 상관 없음.
    """
    keys = set(list(P) + list(Q))
    M = {}
    for key in keys:
        M[key] = 0.5*P.get(key,0) + 0.5*Q.get(key,0)
    
    KLD_PM = calculate_KLD_on_discrete_distribution(P,M,base_of_log)
    KLD_QM = calculate_KLD_on_discrete_distribution(Q,M,base_of_log)
    
    return 0.5*(KLD_PM + KLD_QM)
