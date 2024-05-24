import numpy as np
import networkx as nx

def SBM(N,M,q0,q1):
    '''This function is designed to generate the Stochastic Block Model.
    input params (consistent with the project description):
    sizes (list of int): the sizes of the communities
    q0 (float): probability of inter-subgroup connection
    q1 (float): probability of intra-subgroup connection

    output:
    G (N*N): adjacency matrix of the generated graph
    '''
    #################################################
    
    blocks = []
    approx_block_size = N // M
    
    for i in range(M):
        blocks.append(approx_block_size)
    
    total_approx_size = approx_block_size * M
    while total_approx_size < N:
        index = 0
        blocks[index] = blocks[index] + 1
        index = index + 1
        total_approx_size = total_approx_size + 1
        
    G = np.zeros((N, N))
    
    for i in range(M):
        for j in range(M):
            if i == j:
                q = q0
            else:
                q = q1
                
            start_i = sum(blocks[:i])
            end_i = start_i + blocks[i]
            start_j = sum(blocks[:j])
            end_j = start_j + blocks[j]
            
            for row in range(start_i, end_i):
                for col in range(start_j, end_j):
                    if row == col:
                        G[row, col] = 0
                        continue
                        
                    if G[row, col] == 1:
                        continue
                    
                    if np.random.rand() < q:
                        G[row, col] = 1
                        G[col, row] = 1
                    else:
                        G[row, col] = 0

    #################################################

    return G
