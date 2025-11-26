class Network_state:
    """network state is saved as an int format.
    conversion from int format to other format is done by this class methods."""
    def __init__(self, dynamics_pyboolnet):
        self.dynamics_pyboolnet = dynamics_pyboolnet
        self.state_int_form = None
    
    def __repr__(self):
        return "model state {}".format(self.get_state_dict_form())
    


    #
    # These functions are used to insert and extract state values from this state object.
    # Since the state object stores network states by converting them into an integer,
    # these functions are required when storing a state in a familiar format
    # or retrieving the state information from the object in a familiar format.
    #
    def get_state_dict_form(self):
        return self._convert_int_form_to_dict_form(self.dynamics_pyboolnet, self.state_int_form)
    
    def put_state_int_form(self, int_form):
        """ when this method is used, the 'inf_form' should be achieved by methods in this class"""
        self.state_int_form = int_form
    
    def put_state_list_form(self, list_form):
        """The list_form represents a list containing the state values of the nodes.
        These node state values must be ordered according to the node sequence
        defined in `self.dynamics_pyboolnet.get_node_names()`.""" 
        self.state_int_form = self._convert_list_form_to_int_form(list_form)

    def put_state_dict_form(self, dict_form):
        self.state_int_form = self._convert_dict_form_to_int_form(self.dynamics_pyboolnet, dict_form)

    #
    # These functions apply perturbations to the state object,
    # modifying the state of specific nodes to the desired values.
    # The network state is altered by the perturbation, but it is not fixed permanently.
    #
    def apply_perturbation(self, perturbation:{}, make_new_obj=False):
        """if 'make_new_obj' is True, new network state object with perturbed state is made and returned"""
        dict_form_state = self.get_state_dict_form()
        dict_form_after_perturbation = self.apply_perturbation_to_dict_form_state(dict_form_state, perturbation)
        if make_new_obj:
            perturbed_state = Network_state(self.dynamics_pyboolnet)
            perturbed_state.put_state_dict_form(dict_form_after_perturbation)
            return perturbed_state
        else:
            self.put_state_dict_form(dict_form_after_perturbation)

    @classmethod
    def apply_perturbation_to_dict_form_state(cls, dict_form, perturbation:{}):
        """dict_form is copied, and overrided by perturbation, then returned."""
        dict_form_with_perturbation = dict_form.copy()
        for node, state in perturbation.items():
            dict_form_with_perturbation[node] = state
        return dict_form_with_perturbation
    

    #
    # First, each state is stored as its integer representation. 
    # When needed, you can convert it to a dict form or list form using the convert functions here.
    # The dict form is structured as {node_name: state_value, ...}.
    # The list form is [state1, state2, ...], where the i-th element corresponds
    # to the state of the i-th node returned by dynamics_pyboolnet.get_node_names().
    #
    @classmethod
    def _convert_int_form_to_dict_form(cls, dynamics_pyboolnet, int_form):
        list_form = cls._convert_int_form_to_list_form(dynamics_pyboolnet, int_form)
        return dict(zip(dynamics_pyboolnet.get_node_names(), list_form))
    
    @classmethod
    def _convert_int_form_to_list_form(cls, dynamics_pyboolnet, int_form):
        list_form = [0]*len(dynamics_pyboolnet.get_node_names())
        i = 0
        while int_form:
            list_form[i] = int_form%2
            int_form = int_form >>1
            i += 1
        
        return list_form
    
    @classmethod
    def _convert_list_form_to_int_form(cls, list_form):
        int_form = 0
        for i, state in enumerate(list_form):
            int_form += state*pow(2,i)
        
        return int_form
    
    @classmethod
    def _convert_dict_form_to_int_form(cls, dynamics_pyboolnet, dict_form):
        list_form =[]
        for node in dynamics_pyboolnet.get_node_names():
            list_form.append(dict_form[node])
        return cls._convert_list_form_to_int_form(list_form)
    
    #
    # These magic methods are used to compare the integer-converted values
    # between state objects. This allows us to store states in a list using
    # the bisect function, making it easy to check whether a newly created
    # state object already exists in the list of discovered states.
    #
    def __eq__(self, int_form_or_obj):
        if isinstance(int_form_or_obj,Network_state):
            return self.state_int_form == int_form_or_obj.state_int_form
        elif isinstance(int_form_or_obj,int):
            return self.state_int_form == int_form_or_obj
    
    def __lt__(self, int_form_or_obj):
        if isinstance(int_form_or_obj,Network_state):
            return self.state_int_form < int_form_or_obj.state_int_form
        elif isinstance(int_form_or_obj,int):
            return self.state_int_form < int_form_or_obj
        
    def __le__(self, int_form_or_obj):
        if isinstance(int_form_or_obj,Network_state):
            return self.state_int_form <= int_form_or_obj.state_int_form
        elif isinstance(int_form_or_obj,int):
            return self.state_int_form <= int_form_or_obj
        
    def __gt__(self, int_form_or_obj):
        if isinstance(int_form_or_obj,Network_state):
            return self.state_int_form > int_form_or_obj.state_int_form
        elif isinstance(int_form_or_obj,int):
            return self.state_int_form > int_form_or_obj
    
    def __ge__(self, int_form_or_obj):
        if isinstance(int_form_or_obj,Network_state):
            return self.state_int_form >= int_form_or_obj.state_int_form
        elif isinstance(int_form_or_obj,int):
            return self.state_int_form >= int_form_or_obj
        

class Attractor:
    """object of this class save the information of attractor.
    in here, the attractor emerged from synchronous update model."""
    def __init__(self, dynamics_pyboolnet):
        self.dynamics_pyboolnet = dynamics_pyboolnet
        self.attractor_states = [] # elements of this list become Network_state objects
        self.perturbation = {}

    #
    # methods to get information from this attractor object
    #
    def is_point_attractor(self):
        return len(self.attractor_states) == 1
    
    def show_states(self):
        for network_state in self.attractor_states:
            print(str(network_state))
    
    def get_attractor_states(self):
        return self.attractor_states
    
    def get_average_state(self):
        """For each node in the network model, compute the average of its state values
        across all network states in the attractor. Return the result as a dict."""
        node_statesum = {}
        for network_state in self.attractor_states:
            dict_form = network_state.get_state_dict_form()
            for node, state in dict_form.items():
                node_statesum[node] = node_statesum.setdefault(node, 0) + state
        
        node_averagestate = {node: statesum/len(self.attractor_states) for node, statesum in node_statesum.items()}
        return node_averagestate
    
    #
    # These methods are used to insert the attractor information
    #
    def put_attractor_states_in_synchro_using_network_state_forms(self, network_states_of_attractor=[], 
                                                              perturbation={}):
        """insert attractor information obtained from a synchronous update.
        Each network state must be a Network_state object, stored in a list
        in the exact transition order within the attractor."""
        self.update_type = "synchronous"
        self.perturbation = perturbation.copy()

        first_state_in_cycle = network_states_of_attractor.index(min(network_states_of_attractor))
        self.attractor_states = network_states_of_attractor[first_state_in_cycle:] + network_states_of_attractor[:first_state_in_cycle]
        # A cyclic attractor stores multiple network state objects in self.attractor_states.
        # Depending on which network state is stored first, the same cyclic attractor
        # may result in different self.attractor_states representations.
        # To make identical cyclic attractors have the same representation,
        # we place the network state with the smallest integer-converted value
        # at the beginning of the list.

    def put_attractor_states_in_synchro_using_dict_state_forms(self, dict_form_network_states_of_attractor:[], perturbation={}):
        """insert attractor information obtained from a synchronous update.
        Each network state must be in dict form, stored in a list
        in the exact transition order within the attractor."""
        network_states = []
        for state_dict_form in dict_form_network_states_of_attractor:
            network_state_obj = Network_state(self.dynamics_pyboolnet)
            network_state_obj.put_state_dict_form(state_dict_form)
            network_states.append(network_state_obj)
        #dictionary describing each network state is converted to Network_state object
        self.put_attractor_states_in_synchro_using_network_state_forms(network_states, perturbation)


    #
    # These methods are used to compare attractor objects.
    #    
    def __eq__(self, att_obj):
        return self.attractor_states == att_obj.attractor_states
    
    def extract_non_cyclic_states(self):
        """Identify non-fluctuating nodes in the attractor and return them with their non-fluctuating value in dict form.
        After computing the average state, if its value is either 0 or 1,
        the node is classified as non-fluctuating."""
        node_averagestate_map = self.get_average_state()
        node_nonfluctuatingstate_map = {}

        for node, averagestate in node_averagestate_map.items():
            if averagestate in (0,1):
                node_nonfluctuatingstate_map[node] = int(averagestate)
        
        return node_nonfluctuatingstate_map