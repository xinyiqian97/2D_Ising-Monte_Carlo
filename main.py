import numpy as np
import plots

"""
Initialize N*N lattice with different setting
1 - all empty
2 - all occupied
3 - ordered (0.5 occupancy)
4 - random
"""
def initialize(N, Config_start, seed = 0):
    if Config_start == 1:
        config = - np.ones((N, N))
    if Config_start == 2:
        config = np.ones((N, N))
    if Config_start == 3:
        config = np.ones((N, N))
        for i in range(1, N, 1):
            offset = 1 + i % 2
            for j in range(offset, N, 2):
                config[i][j] = -1
    if Config_start == 4:
        np.random.seed(seed)
        config = 2 * np.random.randint(2, size=(N, N)) - np.ones((N,N))
        
    return config
"""
def periodic(config, N):
    config[0,:] = config[N,:]
    config[N + 1,:] = config[1,:]
    config[:,0] = config[:,N]
    config[:,N + 1] = config[:,1]
    config[0][0] = config[N][N]
    config[0][N + 1] = config[1][N]
    config[N + 1][0] = config[N][1]
    config[N + 1][N + 1] = config[1][1]
    return
"""

"""
Energy for whole system
E = \sum V1 * NN + \sum V2 * NNN + \sum mu * occupancy
"""
def energy(V1, V2, mu, config):
    E_nn = np.roll(config, 1, axis = 1) + np.roll(config, -1, axis = 1)
    E_nnn = np.roll(E_nn, 1, axis = 0) + np.roll(E_nn, -1, axis = 0)
    E_nn = E_nn + np.roll(config, 1, axis = 0) + np.roll(config, -1, axis = 0)
    E = V1 * np.multiply(E_nn, config) + V2 * np.multiply(E_nnn, config) + mu * config_start
    E_tot = np.sum(E)
    return E_tot

"""
Energy increment when flipping (i,j)
"""
def diff_energy(V1, V2, mu, i, j, contig, l = -1, r = 1, u = -1, d = 1):
    dE_nn = contig[i + u][j] + contig[i + d][j] + contig[i][j + l] + contig[i][j + r]
    dE_nnn = contig[i + u][j + l] + contig[i + d][j + r] + contig[i + d][j + l] + contig[i + u][j + r]
    dE = V1 * dE_nn + V2 * dE_nnn + mu
    dE = - 2 * contig[i][j] * dE
    return dE                

if __name__ == '__main__':
    """
    Initial Parameters
    """
    N = 32 # number of lattices - N*N grid
    V1 = 1.0 #param for nearest neighbor interaction
    V2 = -2.0 #param for second nearest neighbor interaction
    num_eq = 4000 # equilibrium steps
    interval = 100 # interval for save
    num_interval = 400 #total number of interval
    #temperature at start/final/increment
    t_start = 7.0 
    t_final = 7.0 
    t_inc = -1
    #mu at start/final/increment
    mu_start = -4.0
    mu_final = 0.0
    mu_inc = 0.1
    #flag at start 1 (all -1) or 2 (all +1) or 3 (ordering) or 4 (random)
    Config_start = 4
    Rand_seed = 30
    #flag for whether including temperature inside loop
    T_inside_loop = False
    #flag for whther print all data
    All_print = False
    #flag for print last config
    LastConfig_print = False

    """
    Initialization
    """
    if Config_start == 4:
        config_start = initialize(N, Config_start, Rand_seed)
    else:
        config_start = initialize(N, Config_start)
    energy_start = energy(V1, V2, mu_start, config_start)
    # apply periodic condition
    #if Config_start == 3 or Config_start == 4:
    #    periodic(config_start, N)

    """
    Main Loop
    total number of interval * num_interval sweeps 
    Metropolis Algorithm
    """
    config = config_start
    energy = energy_start
    #Metro
    for m in range(num_interval):
        if All_print:
            print(config)
            print(energy)
        for n in range(interval):
            #apply metropolis
            a = 1
        

    if (LastConfig_print == True):
        print(config)
        print(energy)
        
