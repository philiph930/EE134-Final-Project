import numpy as np
import random

# binary spliting
def binary_splitting_round(s):
    # s: np.array the infectious status & test status
    num = 0
    flag = sum(s[:,0])>0
    assert flag
    stages = 0
    if len(s[:,0])==1:
        s[0,1] = s[0,0]
        return num,s,stages
    
    B1, B2 = np.array_split(s.copy(), 2,axis=0)
    flag = sum(B1[:,0])>0
    num+=1
    stages += 1
    
    if flag:
        n,stmp,stage = binary_splitting_round(B1)
        s[:len(B1),1] = stmp[:,1]
    else:
        s[:len(B1),1] = 0
        n,stmp,stage = binary_splitting_round(B2)
        s[len(B1):,1] = stmp[:,1]
    num += n
    stages += stage
    return num,s,stages 

def binary_splitting(s):
    # modified bs
    # s: 1-d array the infectious status
    st = np.zeros((len(s),2))
    st[:,0] = s
    st[:,1] = np.nan
    nums = 0
    count = sum(np.isnan(st[:,1]))
    stages = 0
    # the undetermined people
    while count!=0:
        mask = np.isnan(st[:,1])
        flag = sum(st[mask,0]>0)>0
        nums += 1
        stages+=1
        if not flag:
            st[mask,1] = 0
        else:
            n,stmp,stage = binary_splitting_round(st[mask,:])
            st[mask,1] = stmp[:,1]
            nums += n
            stages += stage
        count = sum(np.isnan(st[:,1]))
        
    assert sum(st[:,0]!=st[:,1])==0
    return nums,stages, st[:,1]

def _binary_splitting_parallel(st):
    if st.shape[0] == 0:
        return 0,0
    
    nums = 0
    stages = 0
    # flag = sum(st[:,0])>0
    if len(st[:,0])==1:
        st[0,1] = st[0,0]
        return nums,stages
    
    nums += 2
    stages += 1
    B1, B2 = np.array_split(st, 2,axis=0)
    flag1 = sum(B1[:,0])>0
    flag2 = sum(B2[:,0])>0
    if flag1 and flag2:
        n1,s1 = _binary_splitting_parallel(B1)
        n2,s2 = _binary_splitting_parallel(B2)
        nums += n1+n2
        stages += max(s1,s2)
    elif flag1:
        n1,s1 = _binary_splitting_parallel(B1)
        nums += n1
        stages += s1
        B2[:,1] = 0
    elif flag2:
        n2,s2 = _binary_splitting_parallel(B2)
        nums += n2
        stages += s2
        B1[:,1] = 0
    else:
        B1[:,1] = 0
        B2[:,1] = 0

    return nums,stages

def binary_splitting_parallel(s):
    # modified bs
    # s: 1-d array the infectious status
    st = np.zeros((len(s),2))
    st[:,0] = s
    st[:,1] = np.nan
    
    nums,stages = _binary_splitting_parallel(st)
    stages += 1
    assert sum(st[:,0]!=st[:,1])==0
    return nums,stages, st[:,1]

# diag
def diagalg_iter(s):
    # s(np.array): binary string of infection status
    k = int(np.log2(len(s)))
    l = int(2**(k-1))
    lp = 0
    p = np.zeros(k+1)
    group = dict()
    num = np.ones(k+1,dtype=np.int32)
    for i in range(k):
        p[i] = sum(s[lp:lp+l])>0
        group[i] = s[lp:lp+l]
        num[i] = l
        lp+=l
        l = l//2

    p[-1] = s[-1]
    group[k] = np.array([s[-1]])
    # p(array): pattern
    # group(dict): indicate the group information
    # num(array): the group size
    return p.astype(np.int32), group,num


def diag_splitting(s):
    # s(np.array): binary string of infection status
    num_tests = 0
    stages = 0
    pattern, group, nums = diagalg_iter(s)
    stages +=1
    num_tests += len(pattern)
    indices = np.where(pattern == 1)[0]
    flag = 0
    for i in indices:
        if nums[i]>1:
            num_test,stage = diag_splitting(group[i])
            num_tests += num_test
            if not flag:
                stages+=stage
                flag = 1
    return num_tests,stages

#T1 Testing

def binary_split_T1(st, num_infected):
    num_tests = 0
    num_stages = 0

    B1, B2 = np.array_split(st, 2, axis=0)
    B1_sum = sum(B1[:,0])
    B2_sum = sum(B2[:,0])
    
    num_tests += 2
    num_stages += 1
    
    flag1 = B1_sum>0
    flag2 = B2_sum>0
    if flag1 and flag2:
        n1,s1 = splitting(B1, B1_sum, num_infected)
        n2,s2 = splitting(B2, B2_sum, num_infected)
        num_tests += n1+n2
        num_stages += max(s1,s2)
    elif flag1:
        n1,s1 = splitting(B1, B1_sum, num_infected)
        num_tests += n1
        num_stages += s1
        B2[:,1] = 0
    elif flag2:
        n2,s2 = splitting(B2, B2_sum, num_infected)
        num_tests += n2
        num_stages += s2
        B1[:,1] = 0
    else:
        B1[:,1] = 0
        B2[:,1] = 0

    return num_tests, num_stages

def splitting(st, num_infected, prev_infected):
    num_tests = 0
    num_stages = 0
    
    if st.shape[0] == 0:
        return 0,0
    
    if len(st[:,0])==1:
        if prev_infected == -1:
            num_tests += 1
        st[0,1] = st[0,0]
        return num_tests, num_stages
    
    elif num_infected == 0 and prev_infected == -1:
        num_tests += 1
        st[:,1] = st[:,0]
    
    elif prev_infected == -1 or num_infected >= 16:
        k = min(4, int(np.log2(num_infected)))
        k = min(len(st) // 4, k)
        k = max(1, k)

        num_tests += 2 ** (k + 1)
        num_stages += 1
        
        B = np.array_split(st, 2 ** (k + 1), axis=0)

        group_stages = []

        index = 0
        for group in B:
            group_sum = sum(group[:,0])
            
            if np.all(group_sum==0):
                B[index][:,1] = 0
                index += 1
                continue
            else:
                tests, stages = splitting(group, group_sum, num_infected)
                num_tests += tests
                group_stages.append(stages)
            
            index += 1

        num_stages += max(group_stages)
    
    else:
        tests, stages = binary_split_T1(st, num_infected)
        num_tests += tests
        num_stages += stages

    return num_tests, num_stages

def Qtesting1(s):

    if s.shape[0] == 0:
        return 0,0
    
    num_tests = 0
    num_stages = 0
    
    num_stages += 1

    st = np.zeros((len(s),2))
    st[:,0] = s
    st[:,1] = np.nan
    
    num_infected = sum(st[:,0])

    tests, stages = splitting(st, num_infected, -1)
    num_tests += tests
    num_stages += stages
    
    
    return num_tests, num_stages, st[:,1]

#T2 Testing

def range_sum(sum):
    if sum == 0:
        return 0
    
    elif sum == 1:
        return 1
    
    elif sum == 2 or sum == 3:
        return 2
    
    elif sum >= 4 and sum < 8:
        return 3
    
    else:
        return 4

def binary_split_T2(st):
    num_tests = 0
    num_stages = 0

    if st.shape[0] == 0:
        return 0, 0

    if st.shape[0] == 1:
        st[:, 1] = st[:, 0]
        return 1, 1

    B1, B2 = np.array_split(st, 2, axis=0)
    B1_range = range_sum(sum(B1[:,0]))
    B2_range = range_sum(sum(B2[:,0]))
    
    num_tests += 2
    num_stages += 1
    
    flag1 = B1_range > 0
    flag2 = B2_range > 0

    if flag1 and flag2:
        n1, s1 = binary_split_T2(B1)
        n2, s2 = binary_split_T2(B2)
        num_tests += n1 + n2
        num_stages += max(s1, s2)
    elif flag1:
        n1, s1 = binary_split_T2(B1)
        num_tests += n1
        num_stages += s1
        B2[:, 1] = 0
    elif flag2:
        n2, s2 = binary_split_T2(B2)
        num_tests += n2
        num_stages += s2
        B1[:, 1] = 0
    else:
        B1[:, 1] = 0
        B2[:, 1] = 0

    return num_tests, num_stages
    
def range_testing(st, num_infected, first_iteration):
    num_tests = 0
    num_stages = 0

    if st.shape[0] == 0:
        return 0, 0
    
    if st.shape[0] == 1:
        if first_iteration == 1:
            num_tests += 1
        st[0, 1] = st[0, 0]
        return num_tests, num_stages

    interval = range_sum(num_infected)
    
    if interval == 0:
        st[:, 1] = 0
    
    elif interval == 1:
        num_stages += 1
        B = np.array_split(st, 8, axis=0)
        
        for index, group in enumerate(B):
            
            if group.shape[0] == 0:
                continue
            
            num_tests += 1
            group_interval = range_sum(sum(group[:, 0]))
            
            if group_interval == 1:
                tests, stages = binary_split_T2(group)
                num_tests += tests
                num_stages = max(num_stages, stages + 1)
                
                while index + 1 < len(B):
                    index += 1
                    if B[index].shape[0] != 0:
                        B[index][:, 1] = 0
                break
                
            else:
                B[index][:, 1] = 0
            
    elif interval == 2: 
        infected_index = []
        B = np.array_split(st, 8, axis=0)

        num_stages += 1
        for index, group in enumerate(B):
            if group.shape[0] == 0:
                continue

            num_tests += 1
            num_infected = sum(group[:, 0])
            group_interval = range_sum(num_infected)

            if group_interval == 1 or group_interval == 2:
                infected_index.append(index)
            else:
                B[index][:, 1] = 0

        for indice in infected_index:
            tests, stages = range_testing(B[indice], sum(B[indice][:, 0]), 0)
            num_tests += tests
            num_stages = max(num_stages, stages + 1)

    elif interval == 3:
        infected_index = []
        B = np.array_split(st, 16, axis=0)

        num_stages += 1
        for index, group in enumerate(B):
            if group.shape[0] == 0:
                continue

            num_tests += 1
            group_interval = range_sum(sum(group[:, 0]))

            if group_interval > 0:
                infected_index.append(index)
            else:
                B[index][:, 1] = 0

        for indice in infected_index:
            tests, stages = range_testing(B[indice], sum(B[indice][:, 0]), 0)
            num_tests += tests
            num_stages = max(num_stages, stages + 1)

    else:
        num_stages += 1
        B = np.array_split(st, 2, axis=0)
        
        infected_index = []
        for index, group in enumerate(B):
            num_tests += 1
            group_interval = range_sum(sum(group[:, 0]))

            if group_interval > 0:
                infected_index.append(index)
            else:
                B[index][:, 1] = 0
            
        for indice in infected_index:
            tests, stages = range_testing(B[indice], sum(B[indice][:, 0]), 0)
            num_tests += tests
            num_stages = max(num_stages, stages + 1)
        
    return num_tests, num_stages
  
def Qtesting2(s):
    '''
    s(np.array): binary string of infection status
    '''

    st = np.zeros((len(s), 2))
    st[:, 0] = s
    st[:, 1] = np.nan
    
    if st.shape[0] == 0:
        return 0, 0, st[:, 1]
    
    num_tests = 0
    num_stages = 0

    num_tests += 1
    num_stages += 1

    num_infected = sum(st[:, 0])
    tests, stages = range_testing(st, num_infected, 1)
    
    num_tests += tests
    num_stages = max(num_stages, stages + 1)
    return num_tests, num_stages, st[:, 1]

    ###################################################

def identify_comm_infected1(st, infected_comm, portion):
    num_tests = 0
    num_stages = 0

    # If the community is small enough, test each individual
    if len(infected_comm) <= 6:
        num_stages += 1
        for individual in infected_comm:
            num_tests += 1
            st[individual, 1] = st[individual, 0]
        
        return num_tests, num_stages
    
    # Hierarchical testing for larger communities
    
    if portion > 0.5:
        stage_list = []
        infected_subgroups = np.array_split(infected_comm, 8)

        index = 0
        for group in infected_subgroups:
            
            if len(group) == 0:
                index += 1
                continue
                
            individual_index = 0
            comm_status = np.zeros(len(group))
            for individual in group:
                comm_status[individual_index] = st[individual, 0]
                individual_index += 1
            index += 1

            tests, stages, arr = Qtesting1(comm_status)
            num_tests += tests
            stage_list.append(stages)

            arr_index = 0
            if len(arr) != 0:
                for individual in group:
                    st[individual, 1] = arr[arr_index]
                    arr_index += 1

        num_stages += max(stage_list)
        
    elif portion > 0.2 and portion < 0.5:
        stage_list = []
        infected_subgroups = np.array_split(infected_comm, 6)

        index = 0
        for group in infected_subgroups:
            
            if len(group) == 0:
                index += 1
                continue
                
            individual_index = 0
            comm_status = np.zeros(len(group))
            for individual in group:
                comm_status[individual_index] = st[individual, 0]
                individual_index += 1
            index += 1

            tests, stages, arr = Qtesting1(comm_status)
            num_tests += tests
            stage_list.append(stages)

            arr_index = 0
            if len(arr) != 0:
                for individual in group:
                    st[individual, 1] = arr[arr_index]
                    arr_index += 1

        num_stages += max(stage_list)
        
    else:
        stage_list = []
        infected_subgroups = np.array_split(infected_comm, 4)

        index = 0
        for group in infected_subgroups:
            
            if len(group) == 0:
                index += 1
                continue
                
            individual_index = 0
            comm_status = np.zeros(len(group))
            for individual in group:
                comm_status[individual_index] = st[individual, 0]
                individual_index += 1
            index += 1

            tests, stages, arr = Qtesting1(comm_status)
            num_tests += tests
            stage_list.append(stages)

            arr_index = 0
            if len(arr) != 0:
                for individual in group:
                    st[individual, 1] = arr[arr_index]
                    arr_index += 1

        num_stages += max(stage_list) 
        
    return tests, stages

def Qtesting1_comm_aware(s,communities):
    '''
    s(np.array): binary string of infection status
    communities(list): the community information
    '''
    num_tests = 0
    num_stages = 0

    st = np.zeros((len(s), 2))
    st[:, 0] = s
    st[:, 1] = np.nan
    
    if st.shape[0] == 0:
        return 0, 0, st[:, 1]
    
    num_stages += 1
    infected_communities = []
    infected_portion_list = []
    infected_portion = 0
    
    for community in communities:
        infected = 0
        num_tests += 1
        for individual in community:
            if st[individual, 0] == 1:
                infected_portion += 1
                infected = 1
                break
                
        infected_portion = infected_portion / len(community)
        
        if infected == 1:
            infected_communities.append(community)
            infected_portion_list.append(infected_portion)
        else:
            for individual in community:
                st[individual, 1] = 0

    stages_list = []
    portion_index = 0
    for infected_community in infected_communities:
        tests, stages = identify_comm_infected1(st, infected_community, infected_portion_list[portion_index])
        num_tests += tests
        portion_index += 1
        stages_list.append(stages) 
                
    if stages_list:
        num_stages += max(stages_list)

    return num_tests, num_stages, st[:,1]

def identify_comm_infected2(st, infected_comm, range):
    num_tests = 0
    num_stages = 0

    # If the community is small enough, test each individual
    if len(infected_comm) <= 6:
        num_stages += 1
        for individual in infected_comm:
            num_tests += 1
            st[individual, 1] = st[individual, 0]
        
        return num_tests, num_stages
    
    # Hierarchical testing for larger communities
    
    if len(infected_comm) < (2**(range + 1)) and range != 4:
        stage_list = []
        infected_subgroups = np.array_split(infected_comm, 8)

        index = 0
        for group in infected_subgroups:
            
            if len(group) == 0:
                index += 1
                continue
                
            individual_index = 0
            comm_status = np.zeros(len(group))
            for individual in group:
                comm_status[individual_index] = st[individual, 0]
                individual_index += 1
            index += 1

            tests, stages, arr = Qtesting2(comm_status)
            num_tests += tests
            stage_list.append(stages)

            arr_index = 0
            if len(arr) != 0:
                for individual in group:
                    st[individual, 1] = arr[arr_index]
                    arr_index += 1

        num_stages += max(stage_list)
        
    elif len(infected_comm) < (2**(range + 2)) and range != 4:
        stage_list = []
        infected_subgroups = np.array_split(infected_comm, 6)

        index = 0
        for group in infected_subgroups:
            
            if len(group) == 0:
                index += 1
                continue
                
            individual_index = 0
            comm_status = np.zeros(len(group))
            for individual in group:
                comm_status[individual_index] = st[individual, 0]
                individual_index += 1
            index += 1

            tests, stages, arr = Qtesting2(comm_status)
            num_tests += tests
            stage_list.append(stages)

            arr_index = 0
            if len(arr) != 0:
                for individual in group:
                    st[individual, 1] = arr[arr_index]
                    arr_index += 1

        num_stages += max(stage_list)
        
    else:
        stage_list = []
        infected_subgroups = np.array_split(infected_comm, 4)

        index = 0
        for group in infected_subgroups:
            
            if len(group) == 0:
                index += 1
                continue
                
            individual_index = 0
            comm_status = np.zeros(len(group))
            for individual in group:
                comm_status[individual_index] = st[individual, 0]
                individual_index += 1
            index += 1

            tests, stages, arr = Qtesting2(comm_status)
            num_tests += tests
            stage_list.append(stages)

            arr_index = 0
            if len(arr) != 0:
                for individual in group:
                    st[individual, 1] = arr[arr_index]
                    arr_index += 1

        num_stages += max(stage_list) 
        
    return tests, stages

def Qtesting2_comm_aware(s, communities):
    num_tests = 0
    num_stages = 0

    st = np.zeros((len(s), 2))
    st[:, 0] = s
    st[:, 1] = np.nan
    
    if st.shape[0] == 0:
        return 0, 0, st[:, 1]
    
    num_stages += 1
    infected_communities = []
    infected_range_list = []
    infected_range = 0
    
    for community in communities:
        infected = 0
        num_tests += 1
        for individual in community:
            if st[individual, 0] == 1:
                infected_range += 1
                infected = 1
                break
                
        infected_range = range_sum(infected_range)
        
        if infected == 1:
            infected_communities.append(community)
            infected_range_list.append(infected_range)
        else:
            for individual in community:
                st[individual, 1] = 0

    stages_list = []
    range_index = 0
    for infected_community in infected_communities:
        tests, stages = identify_comm_infected2(st, infected_community, infected_range_list[range_index])
        num_tests += tests
        range_index += 1
        stages_list.append(stages) 
                
    if stages_list:
        num_stages += max(stages_list)

    return num_tests, num_stages, st[:,1]