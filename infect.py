import numpy as np
import random

def infect_step(G,p1,individuals,N):
    '''The function serves as the infection model for each day.
    input params (consistent with the project description):
    G (ndarray N*N): the adjacency matrix.
    p1: the probability each individual infects neighbours.
    '''

    ###################################################
    
    individuals_updated = individuals

    original_infected = []
    for i in range(N):
        if individuals[i] == 1:
            original_infected.append(i)

    for individual in original_infected:
        neighbors = []
        node_num = 0
        for edge in G[individual,:]:
            if edge == 1:
                neighbors.append(node_num)
            node_num = node_num + 1

        for neighbor in neighbors:
            if individuals_updated[neighbor] == 1:
                continue

            if np.random.rand() < p1:
                individuals_updated[neighbor] = 1
                        
    ###################################################
    return individuals_updated




def infect(G,p0,p1,time_steps):
    '''The function serves as the infection model for each day.
    input params (consistent with the project description):
    G (ndarray N*N): the adjacency matrix.
    p0: the infection probability for initial status.
    p1: the probability each individual infects neighbours.
    time_steps: log N
    '''
    N = G.shape[0]
    individuals = np.zeros(N)
    ###################################################
    for i in range(N):
        if np.random.rand() < p0:
            individuals[i] = 1

    for i in range(time_steps):
        individuals = infect_step(G, p1, individuals, N)

    ###################################################

    return individuals