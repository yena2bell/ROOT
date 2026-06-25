import bisect
import time
import random
import math
import itertools

from Network_state_and_attractor import Network_state, Attractor
import SCC_decomposition
import Wrap_pystablemotifs as wp

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


class Asynchro_att_landscape_for_specific_IC(Attractor_landscape_for_specific_IC):
    """This class inherits from Attractor_landscape_for_specific_IC, 
    but the method of calculating the attractor landscape is based on asynchronous updates."""
    def __init__(self, dynamics_Boolean_net:"Dynamics_pyBoolnet", 
                IC:dict, 
                fixed_node_state_map:dict={},
                repeat_for_each_state=50):
        super().__init__(dynamics_Boolean_net, IC, fixed_node_state_map)
        self.repeat_for_each_state = repeat_for_each_state

        self.attindex_initial_state_arrival_count_map = {}
    
    @property
    def attindex_basinratio_map(self):
        """Compute and return the ratio of basin sizes from self.attindex_initial_state_arrival_count_map.
        The attractor corresponding to each attractor index can be found in
        self.attractor_index_map."""
        if self._attindex_basinstates_map_manual_override is None:
            all_arrival_count = sum(self.attindex_initial_state_arrival_count_map.values())
            attindex_basinratio_map = {index:arrival_count/all_arrival_count for index, arrival_count in self.attindex_initial_state_arrival_count_map.items()}
            return attindex_basinratio_map
        else:
            # if there exists attractor basinstates map which is manually overridden, return that value
            return self._attindex_basinstates_map_manual_override

    @attindex_basinratio_map.setter
    def attindex_basinratio_map(self, value):
        """this setter is needed for method 'get_the_average_state_for_each_IC_by_permanent_control'
        of 'ITP' object."""
        self._attindex_basinstates_map_manual_override = value
    
    def _network_state_is_in_attractor_found(self, network_state_obj):
        """Check if the given network state value is in any of the attractors found so far.
        If it is, return the index of the attractor. 
        If not, return the new index value (length of self.attractor_index_map)."""
        for att_index, attractor in self.attractor_index_map.items():
            if network_state_obj in attractor.get_attractor_states():
                return att_index
        return False

    def calculate_asynchro_att_basinratios_from_all_initial_states(self):
        """For a given IC, compute the asynchronous attractor basin ratios
        starting from all initial states."""
        time_start = time.time()
        num_of_all_initials = pow(2, self.get_num_of_not_fixed_nodes())
        percentage_to_report = 0
        reporting_interval = 0.2

        list_of_network_state_values_detected = []
        # Stores states that have been explored at least once.  
        # Managed using bisect.
        # different from synchronous version, it only takes initial states checked
       
        for network_state in get_all_state(self.dynamics_Boolean_net, 
                                           self.fixed_node_state_map_and_IC, 
                                           list_of_network_state_values_detected):

            if len(list_of_network_state_values_detected) >= num_of_all_initials*percentage_to_report:
                percentage_checked = len(list_of_network_state_values_detected)/num_of_all_initials*100
                print("{}% of all initial states are detected.".format(percentage_checked))
                print("used time: ", time.time()-time_start)
                percentage_to_report += reporting_interval
            
            for _ in range(self.repeat_for_each_state):
                self._analyze_asynchro_trajectory_from_network_state(network_state)
        
        print("time for {}: ".format(str(self)), time.time()-time_start)
    
    def calculate_asynchro_atts_through_time_reversal(self):
        primes = self.dynamics_Boolean_net.get_primes()
        attractor_repertore = wp.find_attractors_asynchronous_update(primes, self.fixed_node_state_map_and_IC)
        attractors = attractor_repertore.get_attractors()

        for i, attractor in enumerate(attractors):
            att_states = attractor.get_states_in_att()
            att_state_objs = []
            for att_state in att_states:
                network_state_obj = Network_state(self.dynamics_Boolean_net)
                network_state_obj.put_state_dict_form(att_state)
                att_state_objs.append(network_state_obj)

            attractor_obj = Attractor(self.dynamics_Boolean_net)
            attractor_obj.put_attractor_states_in_asynchro_using_network_state_forms(att_state_objs, self.fixed_node_state_map)
            self.attractor_index_map[i] = attractor_obj

            # this arrival count is meaninless in this case
            self.attindex_initial_state_arrival_count_map[i] = 1

    def _analyze_asynchro_trajectory_from_network_state(self, network_state):
        """Extend the trajectory from the given `network_state` 
        and determine the attractor that it eventually reaches.

        Use this information to update both `attractor_index_map` 
        and `attindex_initial_state_arrival_count_map`."""
        
        trajectory, att_obj, att_index = self._calculate_asynchro_trajectory_and_check(network_state)

        if att_index not in self.attractor_index_map:
            # this trajectory makes a new attractor
            self.attractor_index_map[att_index] = att_obj

            self.attindex_initial_state_arrival_count_map[att_index] = 1
        else:
            # this trajectory reaches one of the existing attractor basins
            self.attindex_initial_state_arrival_count_map[att_index] += 1

    def _calculate_asynchro_trajectory_and_check(self, 
                                        network_state_value_start:Network_state):
        """For a given state, perform asynchronous updates 
        to generate a trajectory.

        Unless the trajectory reaches a point attractor, 
        the same state cannot appear in consecutive steps.

        The computation terminates when either:
        the trajectory reaches an already known attractor, or
        a new attractor is discovered.

        Returns: 
        (trajectory, associated attractor object, index of that attractor in self.attractor_index_map)

        Algorithm description
        
        1. If the trajectory reaches a previously identified attractor, 
        terminate the computation.
        Find the attractor's index in `self.attractor_index_map` 
        and return it together with the trajectory and attractor object.
        
        2. If the trajectory reaches a newly discovered attractor, 
        return the next available index that would be assigned 
        in `self.attractor_index_map`.
        Since this index does not yet exist 
        in `self.attractor_index_map`, 
        it can be used to determine that a new attractor has been found.

        2.1. If the newly added network state in the trajectory is 
        identical to the immediately preceding network state, 
        that state is a point attractor.
        This is because the asynchronous update function is designed 
        to always produce a state different from the previous one 
        unless no update is possible.

        2.2 If the newly added network state matches a state 
        that has already appeared earlier in the trajectory, 
        there may exist a complex attractor containing that state.
        Identify the feedback loops in the trajectory 
        that lead back to the revisited state and treat them 
        as a single SCC (Strongly Connected Component).
        Then check whether additional states should be included 
        in that SCC.
        If the SCC has no outgoing edges, it is considered 
        a complex attractor.
        If the SCC has outgoing edges, append the SCC states to the trajectory, 
        select one of the outgoing edges, and continue extending the trajectory.
                
        2.3. If the trajectory revisits a previously identified SCC 
        that has outgoing edges, expand the SCC to include 
        the newly discovered feedback paths.
        Then reanalyze downstream states reachable from the SCC 
        and re-evaluate whether the SCC satisfies the conditions 
        of a complex attractor."""
        trajectory = [network_state_value_start]
        # trajectory[i] means the network state after i steps from the starting network state
        # trajectory[0] is the starting network state
        
        SCC_nodes_in_trajectory = []
        # (a,b) in this list means that trajectory[a:b+1] makes a SCC
        # this not assure taht trajectory[a:b+1] is a maximal SCC
        # SCC is considered as a single network state in asynchronous state transition graph
        # `SCC_nodes_in_trajectory` shoule be ordered such that
        # (a,b) = SCC_nodes_in_trajectory[i] and (c,d) = SCC_nodes_in_trajectory[j], 
        # if i<j, then a<b<c<d
        indices_all_downstreams_are_checked = []
        # for SCC `SCC_nodes_in_trajectory[i] = (a,b)`
        # index `indices_all_downstreams_are_checked[i] = k`
        # means that from network state index a to network state k,
        # all downstream nodes of network state is in this SCC

        # test_on = True
        # att_pickle_save = r"D:\new canalizing kernel\ROOT framework\ROOT_application_examples\Asynchronous update approach analysis\tmp_for_att_save"
        # # print("trajectory analysis starts")
        # time_start = time.time()

        network_state_value_next = model_state_asynchronous_update_using_pyboolnet(self.dynamics_Boolean_net, 
                                                                        trajectory[-1], 
                                                                        self.fixed_node_state_map_and_IC)

        while True:            
            index_att_having_state = self._network_state_is_in_attractor_found(network_state_value_next)
            if index_att_having_state is not False:
                # this trajectory reaches one of the existing attractor basins
                # print("")
                # print("converged to already found attractor")
                return trajectory, self.attractor_index_map[index_att_having_state], index_att_having_state
            elif network_state_value_next == trajectory[-1]:
                # this trajectory reaches a new point attractor
                attractor = Attractor(self.dynamics_Boolean_net)
                attractor.put_attractor_states_in_asynchro_using_network_state_forms(trajectory[-1:], self.fixed_node_state_map)
                # if test_on: #for test
                #     print("")
                #     print("converged to new point attractor")
                #     # change self.IC to string that can be used file name
                #     IC_str_for_file_name = "_".join([f"{node}_{state}" for node, state in self.IC.items()])
                #     name_of_att_file = "{}_point_att_{}.pkl".format(IC_str_for_file_name, len(self.attractor_index_map))
                #     pickle.dump(attractor, open(os.path.join(att_pickle_save, name_of_att_file), "wb"))
                return trajectory, attractor, len(self.attractor_index_map)
            elif network_state_value_next in trajectory:
                # this trajectory reaches a state already visited in the same trajectory, 
                # which suggests the possibility of a complex attractor.
                # Therefore, check if there is a complex attractor around the re-visited state.
                # If there is a complex attractor, return the trajectory and the complex attractor.
                index_of_revisited_state = trajectory.index(network_state_value_next)
                SCC_formed = (index_of_revisited_state, len(trajectory)-1)
                index_all_downstream_checked = index_of_revisited_state-1

                # if test_on: #for test
                #     print("\n")
                #     trajectory_int_form = [state.get_state_int_form() for state in trajectory]
                #     print('in trajectory:', trajectory_int_form)
                #     print("feedback loop(SCC) occurs")
                #     print("SCC_formed: ", SCC_formed)
                #     print("net_ints in this SCC is:", trajectory_int_form[SCC_formed[0]:SCC_formed[1]+1])
                #     print("SCC_nodes_in_trajectory before update is")
                #     print(SCC_nodes_in_trajectory)
                #     print("indices_all_downstreams_are_checked before update")
                #     print(indices_all_downstreams_are_checked)

                # newly found feedback (SCC_formed) is combinded 
                # with the previous SCC (SCC_previous) if they are overlapped, 
                # and the combined one is considered as the new SCC_formed.
                SCC_nodes_in_trajectory_new =[]
                indices_all_downstreams_are_checked_new = []

                for SCC_previous, index_checked in zip(SCC_nodes_in_trajectory, indices_all_downstreams_are_checked):
                    if SCC_previous[1] < SCC_formed[0]:
                        SCC_nodes_in_trajectory_new.append(SCC_previous)
                        indices_all_downstreams_are_checked_new.append(index_checked)
                    elif SCC_previous[0] <= SCC_formed[0] and SCC_previous[1] >= SCC_formed[0]:
                        SCC_formed = (SCC_previous[0], SCC_formed[1])
                        index_all_downstream_checked = index_checked
                        # the downstream nodes of SCC_previous are also the downstream nodes of new SCC_formed
                        # so, `index_all_downstream_checked` of this new SCC_formed is same as `index_checked` of SCC_previous
                        SCC_nodes_in_trajectory_new.append(SCC_formed)
                        # `index_all_downstream_checked` of this SCC_formed will be
                        # added to `indices_all_downstreams_are_checked` later
                        break
                    elif SCC_previous[0] > SCC_formed[0]:
                        SCC_nodes_in_trajectory_new.append(SCC_formed)
                        # `index_all_downstream_checked` of this SCC_formed will be
                        # added to `indices_all_downstreams_are_checked` later
                        break
                else:
                    # at first, if SCC_nodes_in_trajectory is empty,
                    # for loop is not done
                    # directly go to this else
                    SCC_nodes_in_trajectory_new.append(SCC_formed)
                SCC_nodes_in_trajectory = SCC_nodes_in_trajectory_new
                indices_all_downstreams_are_checked = indices_all_downstreams_are_checked_new
                
                # SCC_nodes_in_trajectory = [SCC_previous for SCC_previous in SCC_nodes_in_trajectory if SCC_previous[0] < SCC_formed[0]]
                # SCC_nodes_in_trajectory.append(SCC_formed)
                # SCC_previous which are contained in new SCC_formed are all deleted

                # network_states_in_extended_SCC, network_states_directly_downstream_of_SCC = self._search_downstream_of_prev_SCC(trajectory[SCC_formed[0]:])
                network_state_directly_downstraem_of_SCC, index_all_downstream_checked = self._select_downstream_of_prev_SCC(trajectory, SCC_formed, index_all_downstream_checked)
                indices_all_downstreams_are_checked.append(index_all_downstream_checked)
                # if such downstream state does not exist (this SCC is complex attractr)
                # this `network_state_directly_downstraem_of_SCC` becomes None
                # time_relapsed = time.time() - time_start
                # if test_on: #for test
                #     print(f"time relapsed: {time_relapsed:10.3f}, SCCs found: {str(SCC_nodes_in_trajectory):30}", end="\r")

                # if test_on: #for test
                #     print("SCC_nodes_in_trajectory after update")
                #     print(SCC_nodes_in_trajectory)
                #     print("indices_all_downstreams_are_checked after update")
                #     print(indices_all_downstreams_are_checked)
                #     print("SCC_formed_new", SCC_formed)
                #     print("the number of states in this SCC is", SCC_formed[1]-SCC_formed[0]+1)
                
                # if test_on:#for test
                #     network_states_in_extended_SCC_int_form = [state.get_state_int_form() for state in network_states_in_extended_SCC]
                #     #print(network_states_in_extended_SCC)
                #     print("network_states_in_extended_SCC", network_states_in_extended_SCC_int_form)
                #     print("the number of states in extended SCC is ", len(network_states_in_extended_SCC))
                #     print("this SCC has {} number of outgoing edge".format(len(network_states_directly_downstream_of_SCC)))

                # if not network_states_directly_downstream_of_SCC:
                #     # this trajectory reaches a new complex attractor
                #     comp_attractor = Attractor(self.dynamics_Boolean_net)
                #     comp_attractor.put_attractor_states_in_asynchro_using_network_state_forms(network_states_in_extended_SCC, self.fixed_node_state_map)

                #     if test_on:#for test
                #         print("complex attractor is found: ", comp_attractor)
                #         comp_attractor.show_states()
                #         raise ValueError("complex attractor is found: ", comp_attractor)

                #     return trajectory, comp_attractor, len(self.attractor_index_map)
                # else:
                #     if test_on:#for test
                #         print("this SCC is not complex attractor. trajectory is more extended")
                #         print("previous length of trajectory was",len(trajectory))
                #     # trajectory is more extended
                #     trajectory = trajectory[:SCC_formed[0]] + network_states_in_extended_SCC
                #     if test_on:#for test
                #         print("the SCC newly found is added to this trajectory,")
                #         print("so, the length of trajectory becomes",len(trajectory))
                #     network_state_value_next = random.choice(network_states_directly_downstream_of_SCC)

                #     SCC_nodes_in_trajectory = [SCC_previous for SCC_previous in SCC_nodes_in_trajectory if SCC_previous[0] < SCC_formed[0]]
                #     SCC_formed = (SCC_formed[0], len(trajectory)-1)
                #     SCC_nodes_in_trajectory.append(SCC_formed)
                #     # SCC_previous which are contained in new SCC_formed are all deleted

                #     trajectory.append(network_state_value_next)

                if network_state_directly_downstraem_of_SCC is None:
                    # this trajectory reaches a new complex attractor
                    network_states_in_comp_att = trajectory[SCC_formed[0]:SCC_formed[1]+1]
                    comp_attractor = Attractor(self.dynamics_Boolean_net)
                    comp_attractor.put_attractor_states_in_asynchro_using_network_state_forms(network_states_in_comp_att, self.fixed_node_state_map)
                    # if test_on:#for test
                    #     print("")
                    #     print("converged to new complex attractor with {} states".format(len(network_states_in_comp_att)))
                    #     IC_str_for_file_name = "_".join([f"{node}_{state}" for node, state in self.IC.items()])
                    #     name_of_att_file = "{}_comp_att_{}.pkl".format(IC_str_for_file_name, len(self.attractor_index_map))
                    #     pickle.dump(comp_attractor, open(os.path.join(att_pickle_save, name_of_att_file), "wb"))

                    # if test_on:#for test
                    #     print("complex attractor is found: ", comp_attractor)
                    #     comp_attractor.show_states()
                    #     raise ValueError("complex attractor is found: ", comp_attractor)

                    return trajectory, comp_attractor, len(self.attractor_index_map)
                else:
                    # if test_on:#for test
                    #     print("this SCC is not complex attractor. trajectory is more extended")
                    # trajectory is more extended
                    
                    network_state_value_next = network_state_directly_downstraem_of_SCC

            else:
                trajectory.append(network_state_value_next)
                network_state_value_next = model_state_asynchronous_update_using_pyboolnet(self.dynamics_Boolean_net, 
                                                                        trajectory[-1], 
                                                                        self.fixed_node_state_map_and_IC)
    
    def _select_downstream_of_prev_SCC(self, trajectory, SCC_nodes_postions, index_all_downstream_checked):
        networks_states_in_SCC = trajectory[SCC_nodes_postions[0]:SCC_nodes_postions[1]+1]
        for index in range(index_all_downstream_checked+1, SCC_nodes_postions[1]+1):
            downstream_states = all_possible_asynchro_next_states(self.dynamics_Boolean_net, trajectory[index], self.fixed_node_state_map_and_IC)
        
            for downstream_state in downstream_states:
                if downstream_state not in networks_states_in_SCC:
                    return downstream_state, index-1
        else:
            return None, None
    
    def _search_downstream_of_prev_SCC(self, networks_states_in_SCC):
        """this method is not used anymore

        calculate downstream network states for each network state forming SCC
        and check that these downstream network states are also contained
        in the SCC containing the previous SCC
        only check the `self.complex_att_search` number of downstream network states
        
        return list of network states composing extended SCC
        and list of network states not contained the SCC but directly downstream of the SCC"""
        edges_for_SCC_calculation = []
        stack_not_checked_queue = []
        already_checked = set()
        test_on = True
        for state in networks_states_in_SCC:
            if state not in stack_not_checked_queue:
                stack_not_checked_queue.append(state)
                already_checked.add(state.get_state_int_form())
        
        if test_on:
            states_int_copy = already_checked.copy()
            if len(networks_states_in_SCC) > len(already_checked):
                print("###############################")
                print("this SCC has duplicated states!")
                print("###############################")
        
        # set edges within this SCC as very simple feedback loop.
        # to reduce the complexity
        for i_state in range(len(stack_not_checked_queue)-1):
            state_from = stack_not_checked_queue[i_state]
            state_to = stack_not_checked_queue[i_state + 1]
            edge = (state_from.get_state_int_form(), state_to.get_state_int_form())
            edges_for_SCC_calculation.append(edge)
        edges_for_SCC_calculation.append((stack_not_checked_queue[-1].get_state_int_form(), stack_not_checked_queue[0].get_state_int_form()))
        

        while stack_not_checked_queue:
            network_state = stack_not_checked_queue.pop(0)
            int_form = network_state.get_state_int_form()
            downstream_states = all_possible_asynchro_next_states(self.dynamics_Boolean_net, network_state, self.fixed_node_state_map_and_IC)

            for downstream_state in downstream_states:
                int_form_downstream = downstream_state.get_state_int_form()
                edge = (int_form, int_form_downstream)
                edges_for_SCC_calculation.append(edge)
                if int_form_downstream not in already_checked:
                    stack_not_checked_queue.append(downstream_state)
                    already_checked.add(int_form_downstream)
            
            if len(already_checked) >= self.complex_att_search+len(networks_states_in_SCC):
                break
        
        SCC_extended = SCC_decomposition.get_SCC_containing_the_node(networks_states_in_SCC[0].get_state_int_form(), edges_for_SCC_calculation)
        # its elements are all int form

        if test_on: #for test
            set_SCC_extended = set(SCC_extended)
            if not states_int_copy.issubset(set_SCC_extended):
                print("###############################")
                print("some states in this SCC are not included in the extended SCC!")
                print("###############################")
                print("previous states")
                print(states_int_copy)
                print("extended SCC states")
                print(set_SCC_extended)
                print("from\tto")
                for edge in edges_for_SCC_calculation:
                    print("{}\t{}".format(edge[0],edge[1]))
                raise ValueError
        SCC_extended_network_state_form = []
        for int_form in SCC_extended:
            network_state_in_SCC = Network_state(self.dynamics_Boolean_net)
            network_state_in_SCC.put_state_int_form(int_form)
            SCC_extended_network_state_form.append(network_state_in_SCC)
        
        int_form_outgoing_from_SCC = []
        for edge_int_form in edges_for_SCC_calculation:
            if edge_int_form[0] in SCC_extended and edge_int_form[1] not in SCC_extended:
                int_form_outgoing_from_SCC.append(edge_int_form[1])
        network_states_outgoing_from_SCC = []
        for int_form in int_form_outgoing_from_SCC:
            network_state_outgoing_from_SCC = Network_state(self.dynamics_Boolean_net)
            network_state_outgoing_from_SCC.put_state_int_form(int_form)
            network_states_outgoing_from_SCC.append(network_state_outgoing_from_SCC)
        
        return SCC_extended_network_state_form, network_states_outgoing_from_SCC

    def asynchro_converge_network_state_to_attractor(self, network_state:Network_state):
        """Return the attractor reached as a result of network state transitions 
        from the given network state value, 
        under the input configuration and fixed node state map of this attractor landscape."""
        trajectory, att_arrrived, index_of_att = self._calculate_asynchro_trajectory_and_check(network_state)

        if index_of_att not in self.attractor_index_map:
            # this trajectory makes a new attractor
            self.attractor_index_map[index_of_att] = att_arrrived

            self.attindex_initial_state_arrival_count_map[index_of_att] = 1

        return att_arrrived, index_of_att
    
    def asynchro_repeated_converge_network_state_to_attractors(self, network_state:Network_state):
        """Return the attractors reached as a result of repeated network state transitions 
        from the given network state value, 
        under the input configuration and fixed node state map of this attractor landscape.
        This method performs the convergence process multiple times 
        and returns the attractors reachable and probabilities to reach them"""
        index_count_map = {}
        index_att_arrived_map = {}
        for _ in range(self.repeat_for_each_state):
            att_arrived, index_of_att = self.asynchro_converge_network_state_to_attractor(network_state)
            if index_of_att not in index_count_map:
                index_count_map[index_of_att] = 1
                index_att_arrived_map[index_of_att] = att_arrived
            else:
                index_count_map[index_of_att] += 1
        
        index_probability_map = {index:count/self.repeat_for_each_state for index, count in index_count_map.items()}
        return index_att_arrived_map, index_probability_map

    def asynchro_repeated_converge_network_states_to_attractors(self, network_states):
        """return attractors reached as a result of repeated network state transitions 
        from the given network states, 
        under the input configuration and fixed node state map of this attractor landscape.
        the probability of starting each network state is assumed to be uniform."""
        index_att_arrived_map, index_probability_map = {}, {}
        for network_state in network_states:
            index_att_arrived_map_of_this_state, index_probability_map_of_this_state = self.asynchro_repeated_converge_network_state_to_attractors(network_state)
            for index in index_att_arrived_map_of_this_state:
                if index not in index_att_arrived_map:
                    index_att_arrived_map[index] = index_att_arrived_map_of_this_state[index]
                    index_probability_map[index] = index_probability_map_of_this_state[index]/len(network_states)
                else:
                    index_probability_map[index] += index_probability_map_of_this_state[index]/len(network_states)
        return index_att_arrived_map, index_probability_map

    def calculate_asynchro_att_basinratios_from_random_initial_states(self, waiting_num=10000, 
                                                   difference_threshold=0.0001,
                                                   attractors_exist=[],
                                                   verbose=False):
        """For a given IC, compute the asynchronous attractor basin ratios
        starting from random initial states."""
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

                ###########################################################
                # different region from synchronous version
                bisect.insort_left(list_of_network_state_values_detected, network_state)

                for _ in range(self.repeat_for_each_state):
                    self._analyze_asynchro_trajectory_from_network_state(network_state)
                if verbose:
                    print("{} initial states are analyzed.".format(len(list_of_network_state_values_detected)), end="\r")
                ###########################################################

                if len(self.attractor_index_map) > num_of_attractors_found:
                    #new attractor is found
                    i_no_more_att_discovery_count = 0
                    num_of_attractors_found = len(self.attractor_index_map)
                else:
                    i_no_more_att_discovery_count += 1


        print("time for {}: ".format(str(self)), time.time()-time_start)


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
    all_nodes_names =dynamics_pyboolnet.get_node_names()
    
    if fixed_node_state_map:
        all_nodes_not_fixed = tuple([node for node in all_nodes_names if node not in fixed_node_state_map])
        num_of_nodes = len(all_nodes_not_fixed)
        for all_state_combination in itertools.product((0,1), repeat=num_of_nodes):
            dict_form = {node:state for node, state in zip(all_nodes_not_fixed, all_state_combination)}
            dict_form = dict(**dict_form, **fixed_node_state_map)
            int_form_perturbed = Network_state._convert_dict_form_to_int_form(dynamics_pyboolnet, dict_form)

            position_to_put = _check_existence_of_int_in_list(int_form_perturbed, states_except)
            # if this `int_form_perturbed` is already in `states_except`, return None
            # else, return the position to put it
            if position_to_put is not None:
                network_state_value = Network_state(dynamics_pyboolnet)
                network_state_value.put_state_int_form(int_form_perturbed)
                yield network_state_value
    
    else:
        num_of_nodes = len(all_nodes_names)
        num_of_possible_states = pow(2, num_of_nodes)
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

def model_state_asynchronous_update_using_pyboolnet(dynamics_pyboolnet, 
                                                   network_state_value_now:Network_state, 
                                                   node_fixed_state_map={}):
    """From the current model state (`network_state_value_now`),
    first applies `node_fixed_state_map`,
    then performs a asynchronous update and returns
    the resulting next network state value.
    
    Thus, when non-fixed nodes are updated,
    any fixed regulator nodes influence the update
    through their fixed values."""
    network_state_value_dict_form = network_state_value_now.get_state_dict_form()
    network_state_value_dict_form = network_state_value_now.apply_perturbation_to_dict_form_state(network_state_value_dict_form, node_fixed_state_map)
    # `network_state_value_now` is overridden by `node_fixed_state_map
    
    # 1. Filter only the non-fixed nodes.
    non_fixed_nodes = [node for node in network_state_value_dict_form.keys() if node not in node_fixed_state_map]
    
    # 2. Identify nodes that would actually cause a state change if updated (candidates).
    changing_candidates = []
    updated_values = {} # Store the results in advance to avoid redundant computations.
    
    for node in non_fixed_nodes:
        # Use an external function to compute the next value of each node.
        next_val = node_state_update_using_pyboolnet(dynamics_pyboolnet, network_state_value_dict_form, node)
        
        # Add a node to the candidate list 
        # only if its next value differs from its current value.
        if next_val != network_state_value_dict_form[node]:
            changing_candidates.append(node)
            updated_values[node] = next_val

    # 3. If there are candidate nodes that would cause a state change:
    if changing_candidates:
        # Randomly select one of them.
        chosen_node = random.choice(changing_candidates)
        
        # Create and return a new state dictionary reflecting the update.
        network_state_value_dict_form[chosen_node] = updated_values[chosen_node]
        
    # 4. If no non-fixed nodes caused a state change (Stg stagnation)
    # return the existing state.
    network_state_value_next = Network_state(dynamics_pyboolnet)
    network_state_value_next.put_state_dict_form(network_state_value_dict_form)
    return network_state_value_next

def all_possible_asynchro_next_states(dynamics_pyboolnet, 
                                    network_state_value_now:Network_state, 
                                    node_fixed_state_map={}):
    """Returns a list of all possible next states after an asynchronous update.
    if the state is in point attractor, it returns empty list"""
    possible_next_states = []
    
    network_state_value_dict_form = network_state_value_now.get_state_dict_form()
    network_state_value_dict_form = network_state_value_now.apply_perturbation_to_dict_form_state(network_state_value_dict_form, node_fixed_state_map)
    # `network_state_value_now` is overridden by `node_fixed_state_map
    
    # 1. Filter only the non-fixed nodes.
    non_fixed_nodes = [node for node in network_state_value_dict_form.keys() if node not in node_fixed_state_map]
    
    # 2. Identify nodes that would actually cause a state change if updated (candidates).
    changing_candidates = []
    updated_values = {} # Store the results in advance to avoid redundant computations.
    
    for node in non_fixed_nodes:
        # Use an external function to compute the next value of each node.
        next_val = node_state_update_using_pyboolnet(dynamics_pyboolnet, network_state_value_dict_form, node)
        
        # Add a node to the candidate list 
        # only if its next value differs from its current value.
        if next_val != network_state_value_dict_form[node]:
            changing_candidates.append(node)
            updated_values[node] = next_val
    
    for node, next_val in updated_values.items():
        next_state_dict_form = network_state_value_dict_form.copy()
        next_state_dict_form[node] = next_val
        
        next_state_obj = Network_state(dynamics_pyboolnet)
        next_state_obj.put_state_dict_form(next_state_dict_form)
        possible_next_states.append(next_state_obj)
    
    return possible_next_states

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
    """if index == len(list_to_put), 
    then `value` is appended to `list_to_put`."""
    if index == len(list_to_put):
        list_to_put.append(value)
    else:
        list_to_put.insert(index, value)

def _put_int_values_in_list(values:iter, list_to_put):
    """Insert the `values` into `list_to_put`, 
    which is already sorted in ascending order,
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
