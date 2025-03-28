import os

import Model_read_using_pyboolnet
from iATG_module import iATG

# 설치해야 하는 패키지
# numpy
# pyBoolnet
# networkx


address_rEMT_bnet = os.path.join('.','example network', 'rEMT_model.bnet')
IC_basal = {"TGFb":0}
IC_transition = {"TGFb":1}
fixed_node_state_map = {"RAS":1}
phenotype_nodes = ["Ecadherin", "ZEB1"]

# random initial states를 사용해서 attractor landscape를 구할 때 필요한 parameter
waiting_num = 1000
# waiting_num이 클 수록 더 많은 initial states를 써서 계산한다.
difference_threshold = 0.00001
# 이 값이 작을수록 attractor landscape를 더 정밀하게 계산한다.

if __name__ == "__main__":
    dynamics_rEMT_network = Model_read_using_pyboolnet.read_pyboolnet_file(address_rEMT_bnet)

    # rEMT에 대한 iATG를 구축하는 부분
    iatg_of_rEMT = iATG(dynamics_rEMT_network, IC_basal, IC_transition, fixed_node_state_map)
    iatg_of_rEMT.set_phenotype_nodes(phenotype_nodes)
    # iatg_of_toy.calculate_attractor_landscapes_for_each_IC(use_all_initials=True)
    iatg_of_rEMT.calculate_attractor_landscapes_for_each_IC(use_all_initials=False,
                                                    waiting_num=waiting_num, 
                                                    difference_threshold=difference_threshold,
                                                    verbose=False)
    iatg_of_rEMT.get_attractor_transitions_induced_by_IC_change_and_calculate_TPs()


    # iATG에 대해, iCAs를 계산한 뒤, major iCAs를 골라내는 단계.
    iatg_of_rEMT.find_iCAs_and_calculate_iCA_sizes()
    print("\n\n계산된 iCAs를 size 별로 보여준다.")
    iatg_of_rEMT.show_iCAs_and_their_sizes()
        

    threshold_for_major_iCA_selection = float(input("major iCA selection을 위한 threshold를 입력하시오 (0 초과, 1 이하):"))
    iatg_of_rEMT.select_major_iCAs_according_to_threshold(threshold_for_major_iCA_selection)
    # 그러면 iatg_of_rEMT.major_iCAs 에 조건에 맞는 iCAs가 들어간다.

    print("selected major iCAs are as follows")
    for ica in iatg_of_rEMT.major_iCAs:
        print(ica)
    print('\n\n')

    # 선택된 major iCAs 에 대해 ITPs를 계산한다.
    ITPs_to_analyze = []
    for i, ica in enumerate(iatg_of_rEMT.major_iCAs):
        print("calculate ITP for",ica)
        ica.search_ITPs()
        
        # 그리고 각 iCA 에서 어떤 ITP를 irreversibility kernel 분석 대상으로 삼을 것인지 정한다.
        for j, itp in enumerate(ica.ITPs):
            print(j,'th ITP of iCA ',ica, ' has phenotype such as')
            itp.get_phenotype_score()
        print("select index of ITP to analyze. 아무것도 선택하지 않으려면 -1을 입력")
        itp_index_to_select = int(input(":"))
        if itp_index_to_select >=0 and itp_index_to_select < len(ica.ITPs):
            selected_itp = ica.ITPs[itp_index_to_select]
            ITPs_to_analyze.append(selected_itp)
            print(itp_index_to_select,"th ITP having following phenotype is selected")
            selected_itp.get_phenotype_score()
    
    # 선택된 ITPs에 대해 irreversibility kernel 을 탐색한다.
    for itp in ITPs_to_analyze:
        itp.find_irreversibility_kernel()
        print("\nirreversibility krenel of ",itp)
        for i, irreversibility_motif in enumerate(itp.irreversibility_motifs):
            coherency_condition = itp.coherency_conditions[i]

            print(i,'th kernel')
            print("\tcoherency condition: ", coherency_condition)
            print("\tirreversibility motif: ", irreversibility_motif)

        print('\nreverse controls in this ITP is')
        reverse_controls = itp.find_reverse_controls()
        print(reverse_controls)





    
