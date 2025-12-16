import bisect
import time
import random
import math

from Network_state_and_attractor import Network_state, Attractor 

class Attractor_landscape_for_specific_IC:
    """calculate attractor landscape for a specific input configuraiton,
    and save the results in the form of attractor index and basin states."""
    def __init__(self, dynamics_Boolean_net:"Dynamics_pyBoolnet", 
                IC:dict, 
                fixed_node_state_map:dict={}):
        self.dynamics_Boolean_net = dynamics_Boolean_net
        self.IC = IC
        #input configuration
        self.fixed_node_state_map = fixed_node_state_map
        #node fixation by mutation

        self.fixed_node_state_map_and_IC = self.fixed_node_state_map.copy()
        self.fixed_node_state_map_and_IC.update(self.IC)

        self.attractor_index_map = {} #{0:attractor_0, 1:attractor_1,...}
        self.attindex_basinstates_map = {}
        # Map the list of network states identified as belonging to a basin
        # to the index of the corresponding attractor.
        self._attindex_basinstates_map_manual_override = None

    def __repr__(self):
        return "Attractor landscape for IC {}".format(self.IC)
    
    @property
    def attindex_basinratio_map(self):
        """Compute and return the ratio of basin sizes from self.attindex_basinstates_map.
        The attractor corresponding to each attractor index can be found in
        self.attractor_index_map."""
        if self._attindex_basinstates_map_manual_override is None:
            attindex_basinsize_map = {index:len(basin) for index, basin in self.attindex_basinstates_map.items()}
            total_basin_size = sum(attindex_basinsize_map.values())
            attindex_basinratio_map = {index:basin_size/total_basin_size for index, basin_size in attindex_basinsize_map.items()}
            return attindex_basinratio_map
        else:
            # if there exists attractor basinstates map which is manually overridden, return that value
            return self._attindex_basinstates_map_manual_override
    
    @attindex_basinratio_map.setter
    def attindex_basinratio_map(self, value):
        """this setter is needed for method 'get_the_average_state_for_each_IC_by_permanent_control'
        of 'ITP' object."""
        self._attindex_basinstates_map_manual_override = value
    
    @property
    def num_of_states_searched(self):
        attindex_basinsize_map = {index:len(basin) for index, basin in self.attindex_basinstates_map.items()}
        return sum(attindex_basinsize_map.values())
    
    def get_num_of_not_fixed_nodes(self):
        """return (num of all nodes - num of (union(fixed ndoes, input nodes)))"""
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
        """For a given IC, compute the attractor basin ratios
        starting from all initial states."""
        time_start = time.time()
        num_of_all_initials = pow(2, self.get_num_of_not_fixed_nodes())
        percentage_to_report = 0
        reporting_interval = 0.2

        list_of_network_state_values_detected = []
        # Stores states that have been explored at least once.  
        # Managed using bisect.  
       

        for network_state in get_all_state(self.dynamics_Boolean_net, 
                                           self.fixed_node_state_map_and_IC, 
                                           list_of_network_state_values_detected):

            if len(list_of_network_state_values_detected) >= num_of_all_initials*percentage_to_report:
                percentage_checked = len(list_of_network_state_values_detected)/num_of_all_initials*100
                print("{}% of all initial states are detected.".format(percentage_checked))
                print("used time: ", time.time()-time_start)
                percentage_to_report += reporting_interval

            self._analyze_trajectory_from_network_state(network_state, list_of_network_state_values_detected)
        
        print("time for {}: ".format(str(self)), time.time()-time_start)
        
    
    def _calculate_trajectory_and_check(self, 
                                        network_state_value_start:Network_state,
                                        list_of_network_state_values_detected:[]):
        """Generate a trajectory starting from `network_state_value_start`.
        Stop when the trajectory either forms a new attractor
        or reaches one of the network states in list_of_network_state_values_detected,
        and return the result.
        """
        trajectory = [network_state_value_start]
        # trajectory[i] means the network state after i steps from the starting network state
        # trajectory[0] is the starting network state
        while True:
            network_state_value_next = model_state_synchronous_update_using_pyboolnet(self.dynamics_Boolean_net, 
                                                                        trajectory[-1], 
                                                                        self.fixed_node_state_map_and_IC)
            
            if network_state_value_next in trajectory:
                # new attractor is found
                return trajectory, network_state_value_next
            elif network_state_value_next in list_of_network_state_values_detected:
                # this trajectory reaches one of the existing attractor basins
                return trajectory, network_state_value_next
            else:
                trajectory.append(network_state_value_next)
    
    def converge_network_state_to_attractor(self, network_state:Network_state):
        """Return the attractor reached as a result of network state transitions 
        from the given network state value, 
        under the input configuration and fixed node state map of this attractor landscape."""
        trajectory, network_state_next = self._calculate_trajectory_and_check(network_state, [])
        return self._make_new_attractor_from_trajectory(trajectory, network_state_next)

    
    def _make_new_attractor_from_trajectory(self, trajectory, network_state_value_next):
        """network_state_value_next already exists within the trajectory.
        Therefore, the network state values from the position of
        network_state_value_next to the end of the trajectory
        constitute the attractor states.
        """
        index = trajectory.index(network_state_value_next)
        attractor_states = trajectory[index:]
        attractor = Attractor(self.dynamics_Boolean_net)
        attractor.put_attractor_states_in_synchro_using_network_state_forms(attractor_states, self.fixed_node_state_map)

        return attractor

            


##########################################################################
# functions for finding multiple initial states to compute the attractor landscape
##########################################################################
def get_random_state(dynamics_pyboolnet, fixed_node_state_map={}, states_except=[]):
    """Randomly samples a network state from the given `dynamics_pyboolnet`,
    taking node state fixations into account.
    
    States included in `states_except` are excluded.
    `states_except` may contain State_in_pyboolnet objects or integer values.
    It is used to avoid sampling states that have already been discovered.
    IMPORTANT: `states_except` must be managed using bisect."""
    num_of_nodes = len(dynamics_pyboolnet.get_node_names())
    num_of_possible_states = pow(2, num_of_nodes)
    while True:
        random_int = random.randint(0, num_of_possible_states-1)
        if fixed_node_state_map:
            dict_form = Network_state._convert_int_form_to_dict_form(dynamics_pyboolnet, random_int)
            dict_form = Network_state.apply_perturbation_to_dict_form_state(dict_form, fixed_node_state_map)
            random_int = Network_state._convert_dict_form_to_int_form(dynamics_pyboolnet, dict_form)

        position_to_put = _check_existence_of_int_in_list(random_int, states_except)
        # if this `random_int` is already in `states_except`, return None
        # else, return the position to put it
        if position_to_put is not None:
            state_obj = Network_state(dynamics_pyboolnet)
            state_obj.put_state_int_form(random_int)

            return state_obj

def get_all_state(dynamics_pyboolnet, fixed_node_state_map={}, states_except=[]):
    """Iterates over all network states in the given `dynamics_pyboolnet`,
    taking node state fixations into account.
    This method is intended to be used as an iterator.
    
    States included in `states_except` are excluded.
    `states_except` may contain State_in_pyboolnet objects or integer values.
    It is used to avoid sampling states that have already been discovered.
    IMPORTANT: `states_except` must be managed using bisect."""
    num_of_nodes = len(dynamics_pyboolnet.get_node_names())
    num_of_possible_states = pow(2, num_of_nodes)
    if fixed_node_state_map:
        for int_form in range(num_of_possible_states):
            dict_form = Network_state._convert_int_form_to_dict_form(dynamics_pyboolnet, int_form)
            dict_form = Network_state.apply_perturbation_to_dict_form_state(dict_form, fixed_node_state_map)
            int_form_perturbed = Network_state._convert_dict_form_to_int_form(dynamics_pyboolnet, dict_form)
        
            position_to_put = _check_existence_of_int_in_list(int_form_perturbed, states_except)
            # if this `int_form_perturbed` is already in `states_except`, return None
            # else, return the position to put it
            if position_to_put is not None:
                network_state_value = Network_state(dynamics_pyboolnet)
                network_state_value.put_state_int_form(int_form_perturbed)
                yield network_state_value
    
    else:
        # the reason for separating this for loop without fixed nodes from the for loop is
        # if there are no fixed nodes, this makes it a little faster
        for int_form in range(num_of_possible_states):
            position_to_put = _check_existence_of_int_in_list(int_form, states_except)
            # if this `int_form` is already in `states_except`, return None
            # else, return the position to put it
            if position_to_put is not None:
                network_state_value = Network_state(dynamics_pyboolnet)
                network_state_value.put_state_int_form(int_form)
                yield network_state_value


# ##########################################################################
# functions related to node and state update by dynamics.
# ##########################################################################
def model_state_synchronous_update_using_pyboolnet(dynamics_pyboolnet, 
                                                   network_state_value_now:Network_state, 
                                                   node_fixed_state_map={}):
    """From the current model state (`network_state_value_now`),
    first applies `node_fixed_state_map`,
    then performs a synchronous update and returns
    the resulting next network state value.
    
    Thus, when non-fixed nodes are updated,
    any fixed regulator nodes influence the update
    through their fixed values."""
    network_state_value_dict_form = network_state_value_now.get_state_dict_form()
    network_state_value_dict_form = network_state_value_now.apply_perturbation_to_dict_form_state(network_state_value_dict_form, node_fixed_state_map)
    # `network_state_value_now` is overridden by `node_fixed_state_map

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
    """Computes the next state of a specific node (whose name is given by `nodename_to_update`)
    from `network_state_value_dict_form`.
    `network_state_value_dict_form` only needs to include
    the regulator nodes of that node."""
    prime_of_node = dynamics_pyboolnet.primes[nodename_to_update]
    prime_for_0 = prime_of_node[0]
    prime_for_1 = prime_of_node[1]
    if len(prime_for_0) >= len(prime_for_1):
        state_to_check = 1
        primes_to_check = prime_for_1
    else:
        state_to_check = 0
        primes_to_check = prime_for_0
    # it is enough to check only one of the primes for 0 or 1
    # so we proceed with the primes that have less to check.

    for prime in primes_to_check:
        for node, state in prime.items():
            if int(network_state_value_dict_form[node]) != state:
                break
        else:
            return state_to_check
        # if the `state_of_model` contains the `prime`,
        # this `node` is updated to `state_to_check`.
    else:
        # if not, this `node` is updated to the opposite value.
        return 1-state_to_check


###############################################################################################
### bisect related functions
###############################################################################################

def _check_existence_of_int_in_list(int_value, list_of_int):
    """`list_of_int` is a list of integers sorted in ascending order.
    
    If `int_value` is not present in `list_of_int`, return the index
    at which it should be inserted.
    If it already exists in the list, return None.
    
    Note that if i == len(list_of_int), the value should be added
    using append() outside of this function."""
    i = bisect.bisect_left(list_of_int, int_value)
    if i != len(list_of_int) and list_of_int[i] == int_value:
        return None
    else:
        return i

def _put_int_value_in_list(value, list_to_put, index):
    """if index == line(list_to_put), then `value` is appended to `list_to_put`."""
    if index == len(list_to_put):
        list_to_put.append(value)
    else:
        list_to_put.insert(index, value)

def _put_int_values_in_list(values:iter, list_to_put):
    """Insert the `values` into `list_to_put`, which is already sorted in ascending order,
    using the bisect method."""
    for value in values:
        bisect.insort_left(list_to_put, value)

###############################################################################################
### Distribution comparison measures related functions
###############################################################################################


def calculate_KLD_on_discrete_distribution(P,Q, base_of_log=math.e, epsilon=0):
    """KLD (Kullback–Leibler Divergence) represents the change in entropy
    when Q is used instead of P.
    
    If a key exists in only one of the distributions P or Q,
    the probability of that key in the other distribution
    is treated as 0.
    
    However, if Q[key] == 0 while P[key] != 0, the KLD is undefined.
    
    To address this issue, a technique called Smoothed KLD is used,
    where a very small value (epsilon) is added so that cases with
    Q[key] == 0 and P[key] != 0 can be handled."""
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
    """Compute the Jensen–Shannon Divergence (JSD).
    JSD(P || Q) = 0.5 * (KLD(P || M) + KLD(Q || M)),
    where M = 0.5 * (P + Q).
    
    JSD is symmetric.
    It is well-defined even when Q[key] == 0 and P[key] != 0."""
    keys = set(list(P) + list(Q))
    M = {}
    for key in keys:
        M[key] = 0.5*P.get(key,0) + 0.5*Q.get(key,0)
    
    KLD_PM = calculate_KLD_on_discrete_distribution(P,M,base_of_log)
    KLD_QM = calculate_KLD_on_discrete_distribution(Q,M,base_of_log)
    
    return 0.5*(KLD_PM + KLD_QM)
