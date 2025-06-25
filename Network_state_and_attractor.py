class Network_state:
    """기본적으로 state를 int form으로 저장하지만, 
    그것을 본래의 state form으로 바꾸는 것을 이 class의 method가 전부 처리하도록.
    
    network state가 가지는 구체적인 값은 network state value라고 부른다."""
    def __init__(self, dynamics_pyboolnet):
        self.dynamics_pyboolnet = dynamics_pyboolnet
        self.state_int_form = None
    
    def __repr__(self):
        return "model state {}".format(self.get_state_dict_form())
    

    #
    # 이 함수들은 state 값을 이 state 객체에 넣고 빼는 용도이다.
    # state 객체는 state를 int 값으로 변환해서 저장하기 때문에,
    # 우리에게 익숙한 형태의 state를 객체에 저장하거나,
    # 객체의 state 정보를 익숙한 형태로 얻어내려면 이 함수들을 사용한다.
    #
    def get_state_dict_form(self):
        return self._convert_int_form_to_dict_form(self.dynamics_pyboolnet, self.state_int_form)
    
    def put_state_int_form(self, int_form):
        """이 함수를 통해 값을 입력할 경우, 이 int form은 이 class의 cls method를 통해 구할 것."""
        self.state_int_form = int_form
    
    def put_state_list_form(self, list_form):
        """state는 self.dynamics_pyboolnet.get_node_names()의 node 순서에 대응된다.""" 
        self.state_int_form = self._convert_list_form_to_int_form(list_form)

    def put_state_dict_form(self, dict_form):
        self.state_int_form = self._convert_dict_form_to_int_form(self.dynamics_pyboolnet, dict_form)

    #
    # state 객체에, 그 state에 perturbation (temporal perturbation)을 가하여
    # 특정 node의 state를 원하는 값으로 바꾸어 놓는데 사용하는 함수들이다.
    #
    def apply_perturbation(self, perturbation:{}, make_new_obj=False):
        """make_new_obj이 True면, 현재 state에 perturbation을 반영한 새로운 obj을 형성하여 return """
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
        """주어진 dict_form에 fixed_node_state_map 을 덧씌워서 return"""
        dict_form_with_perturbation = dict_form.copy()
        for node, state in perturbation.items():
            dict_form_with_perturbation[node] = state
        return dict_form_with_perturbation
    

    #
    # 일단 각 state를 int 형태로 변환하여 저장해놓는다. 그러나 필요할 때, 
    # 여기의 convert 함수들을 사용하여 dict form이나 list form으로 변환하여 쓸 수 있다.
    # dict form의 경우는 {node명: state 값}의 형태로 이루어짐
    # list form의 경우는 [state값1, state값2,,,] 형태이며, 
    # 이 list form의 i번째 state는 dynamics_pyboolnet의 get_node_names method의 결과물의 i번째 node의 state를 의미한다.
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
    # 이 magic methods들은 state 객체들 사이의 int 변환 값들의 대소를 비교하기 위한 methods이다.
    # 이는, 뒤에 이 states를 list에 bisect 함수를 사용하여 저장할 수 있게 하는것으로
    # 새로운 state객체가 이미 찾아진 states list 안에 이미 있는지 확인하는 것을 쉽게 해준다.
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
    """주어진 model의 attractor 정보를 저장하는 객체
    synchronous update를 가정함"""
    def __init__(self, dynamics_pyboolnet):
        self.dynamics_pyboolnet = dynamics_pyboolnet
        self.attractor_states = []
        #여기 들어가는 것은 위의 Network_state 객체들
        self.perturbation = {}

    #
    # attractor에 대한 정보를 얻는데 사용하는 함수들을 모아놓음.
    #
    def is_point_attractor(self):
        return len(self.attractor_states) == 1
    
    def show_states(self):
        for network_state in self.attractor_states:
            print(str(network_state))
    
    def get_attractor_states(self):
        return self.attractor_states
    
    def get_average_state(self):
        """각 node에 대해서, state의 평균을 구한 뒤, dict로 return"""
        node_statesum = {}
        for network_state in self.attractor_states:
            dict_form = network_state.get_state_dict_form()
            for node, state in dict_form.items():
                node_statesum[node] = node_statesum.setdefault(node, 0) + state
        
        node_averagestate = {node: statesum/len(self.attractor_states) for node, statesum in node_statesum.items()}
        return node_averagestate
    
    #
    # 이 attractor 객체에 정보를 넣는데 사용하는 methods 모음.
    #
    def put_attractor_states_in_synchro_using_network_state_forms(self, network_states_of_attractor=[], 
                                                              perturbation={}):
        """synchronous update로 구한 attractor의 정보를 입력한다.
        이 때 각 state는 Network_state 객체로 되어 있으며, 그것이 attractor 내 transition 순서에 맞게 
        list 안에 들어있도록 할 것."""
        self.update_type = "synchronous"
        self.perturbation = perturbation.copy()

        first_state_in_cycle = network_states_of_attractor.index(min(network_states_of_attractor))
        self.attractor_states = network_states_of_attractor[first_state_in_cycle:] + network_states_of_attractor[:first_state_in_cycle]
        #같은 cyclic attractor이면 list의 시작 state가 같게 만드려고 함.
        #state 중, int form으로 변환했을 때 그 값이 가장 작은 state기 list의 가장 처음에 오도록 한다.

    def put_attractor_states_in_synchro_using_dict_state_forms(self, dict_form_network_states_of_attractor:[], perturbation={}):
        """synchronous update로 구한 attractor의 정보를 입력한다.
        각각의 state는 dict로 되어 있으며, 그것이 attractor내 transition 순서에 맞게
        list안에 들어있도록 할 것."""
        network_states = [Network_state(self.dynamics_pyboolnet).put_state_dict_form(state_dict_form) 
                                  for state_dict_form in dict_form_network_states_of_attractor]
        #일단 dict form state를 state 객체로 변환한다.
        self.put_attractor_states_in_synchro_using_state_obj_forms(network_states, perturbation)


    #
    # attractor 객체 사이의 비교와 관련된 methods 들을 모아놓은 부분
    #    
    def __eq__(self, att_obj):
        return self.attractor_states == att_obj.attractor_states
    
    def extract_non_cyclic_states(self):
        """attractor에서 non fluctuating nodes를 찾아서 dict form으로 return
        
        averagestate를 계산한 뒤, 그것이 0이거나 1이면 non fluctuating state라고 판단."""
        node_averagestate_map = self.get_average_state()
        node_nonfluctuatingstate_map = {}

        for node, averagestate in node_averagestate_map.items():
            if averagestate in (0,1):
                node_nonfluctuatingstate_map[node] = int(averagestate)
        
        return node_nonfluctuatingstate_map