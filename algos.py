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
    
    elif prev_infected == -1:
        k = min(4, int(np.log2(num_infected)))
        k = min(len(st) // 4, k)
        k = max(1, k)

        num_tests += 2 ** k
        num_stages += 1
        
        B = np.array_split(st, 2 ** k, axis=0)

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

def Qtesting1_comm_aware(s,communities):
    '''
    s(np.array): binary string of infection status
    communities(list): the community information
    '''
    
def Qtesting1_comm_aware(s, communities):
    num_tests = 0
    stages = 0
    for community in communities:
        num_tests_c, stages_c = binary_splitting(s[community])[:2]
        num_tests += num_tests_c
        stages = max(stages, stages_c)
    return num_tests, stages

def Qtesting2_comm_aware(s, communities):
    num_tests = 0
    stages = 0
    for community in communities:
        num_tests_c, stages_c = diag_splitting(s[community])
        num_tests += num_tests_c
        stages = max(stages, stages_c)
    return num_tests, stages