import numpy as np
import output

"""
Initialize N*N lattice with different setting
1 - all empty
2 - all occupied
3 - ordered (0.5 occupancy)
4 - random
"""
def initialize(N, Config_init, seed = 0):
    if Config_init == 1:
        config = - np.ones((N, N))
    if Config_init == 2:
        config = np.ones((N, N))
    if Config_init == 3:
        config = np.ones((N, N))
        for i in range(0, N, 1):
            offset = i % 2
            for j in range(offset, N, 2):
                config[i][j] = -1
    if Config_init == 4:
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

def coverage(config, N):
    theta = np.sum(config) / (N * N)
    theta = (1 + theta) / 2
    return theta

"""
Energy for whole system
E = \sum V1 * NN + \sum V2 * NNN + \sum mu * occupancy
"""
def H(V1, V2, mu, config):
    E_nn = np.roll(config, 1, axis = 1) + np.roll(config, -1, axis = 1)
    E_nnn = np.roll(E_nn, 1, axis = 0) + np.roll(E_nn, -1, axis = 0)
    E_nn = E_nn + np.roll(config, 1, axis = 0) + np.roll(config, -1, axis = 0)
    E = (1/2) * (V1 * np.multiply(E_nn, config) + V2 * np.multiply(E_nnn, config)) + mu * config_start
    E_tot = np.sum(E)
    return E_tot

"""
Energy increment when flipping (i,j)
"""
def diff_energy(V1, V2, mu, i, j, config, l = -1, r = 1, u = -1, d = 1):
    dE_nn = config[i + u][j] + config[i + d][j] + config[i][j + l] + config[i][j + r]
    dE_nnn = config[i + u][j + l] + config[i + d][j + r] + config[i + d][j + l] + config[i + u][j + r]
    dE = (V1 * dE_nn + V2 * dE_nnn) / 2 + mu
    dE = - 2 * config[i][j] * dE
    return dE

"""
Sweep sample
"""
def sweep_sample(N):
    index_1 = np.random.randint(N, size=(N))
    index_2 = np.random.randint(N, size=(N))
    return index_1, index_2

"""
Sweep
"""
def sweep(V1, V2, mu, N, config_current, T):
    site_i, site_j = sweep_sample(N)
    config = config_current
    for i,j in zip(site_i, site_j):
        if i == 0 or i == N - 1 or j == 0 or j == N - 1:
            l = ((j + N - 1) % N) - j
            r = ((j + 1) % N) - j
            u = ((i + N - 1) % N) - i
            d = ((i + 1) % N) - i
            dE = diff_energy(V1, V2, mu, i, j, config, l, r, u, d)
        else:
            dE = diff_energy(V1, V2, mu, i, j, config)
        if dE <= 0:
            config[i][j] = - config[i][j]
        else:
            crit = np.random.random()
            if T != 0 :
                if np.exp(- dE / T) > crit:
                    config[i][j] = - config[i][j]
    return config

if __name__ == '__main__':
    """
    Initial Parameters
    """
    N = 64 # number of lattices - N*N grid
    V1 = 1.0 #param for nearest neighbor interaction
    V2 = -2.0 #param for second nearest neighbor interaction
    num_eq = 4000 # equilibrium steps
    interval = 50 # interval for save
    num_interval = 400 #total number of interval
    #temperature at start/final/increment
    T_start = 0
    T_final = 5 
    T_inc = 0.5
    #mu at start/final/increment
    mu_start = 0.0
    mu_final = 4.0
    mu_inc = 1
    #flag at start 1 (all -1) or 2 (all +1) or 3 (ordering) or 4 (random)
    Config_init = 4
    Rand_seed = 20
    #flag for whether print all data
    All_print = True
    freq = 20
    #flag for print last config
    LastConfig_print = True
    #run type
    # C certain mu & T
    # P varying mu & T
    # M iterate over mu range
    runtype = 'C'
    
    #grid for plot
    X = np.linspace(1, N, N)
    Y = np.linspace(1, N, N)
    X, Y = np.meshgrid(X, Y)

    """
    Initialization
    """
    if Config_init == 4:
        config_start = initialize(N, Config_init, Rand_seed)
    else:
        config_start = initialize(N, Config_init)
    energy_start = H(V1, V2, mu_start, config_start)
    theta_start = coverage(config_start, N)
    

    """
    Main Loop
    total number of interval * num_interval sweeps 
    Metropolis Algorithm
    """
    config = config_start
    # single run for fixed temperature and mu
    if runtype == 'C':
        energy = [energy_start]
        theta = [theta_start]
        T = T_start
        mu = mu_start
        output.config_save(X, Y, config, 0, interval)
        #Metro
        for m in range(num_interval):
            for n in range(interval):
                config = sweep(V1, V2, mu, N, config, T)
                if All_print:
                    energy.append(H(V1, V2, mu, config))
                    theta.append(coverage(config, N))
                    if m % freq == 0:
                        output.config_save(X, Y, config, m + 1, interval)
        output.energy_save(energy, interval, T)
        output.theta_save(theta, interval, T)
    
        if (LastConfig_print == True):
            output.config_save(X, Y, config, 'last', interval)
    
    # Temperature iterate
    """
    if runtype == 'T':
        E = []
        for T in range(T_start, T_final + T_inc, T_inc):
            for m in range(num_interval - 1):
                for n in range(interval):
                    config = sweep(V1, V2, mu_start, N, config, T)
            energy = [H(V1, V2, mu_start, config)]
            for n in range(interval):
                config = sweep(V1, V2, mu_start, N, config, T)
                energy.append(H(V1, V2, mu_start, config))
            avg_energy = sum(energy) / len(energy)
            E.append(avg_energy)
        output.T_save(E, T_start, T_final, T_inc, mu_start)
    """
    if runtype == 'P':
        num_T_inc = round((T_final - T_start) / T_inc) + 1
        num_mu_inc = round((mu_final - mu_start) / mu_inc) + 1
        T_iter = np.linspace(T_start, T_final, num_T_inc)
        mu_iter = np.linspace(mu_start, mu_final, num_mu_inc)
        for T in T_iter:
            for mu in mu_iter:
                if Config_init == 4:
                    config = initialize(N, Config_init, Rand_seed)
                else:
                    config = initialize(N, Config_init)
                #output.config_save(X, Y, config_start, mu, interval)
                for m in range(num_interval):
                    for n in range(interval):
                        config = sweep(V1, V2, mu, N, config, T)
                theta = coverage(config, N)
                output.P_save(X, Y, config, T, mu, Config_init, theta)            

    # mu iterate
    if runtype == 'M':
        Th = []
        T = [0.2, 0.4, 0.6, 0.8, 1]
        num_inc = round((mu_final - mu_start) / mu_inc) + 1
        M = np.linspace(mu_start, mu_final, num_inc)
        for i in range(num_inc):
            mu = M[i]
            if Config_init == 4:
                config = initialize(N, Config_init, Rand_seed)
            else:
                config = initialize(N, Config_init)
            for m in range(num_interval - 3):
                for n in range(interval):
                    config = sweep(V1, V2, mu, N, config, T_start)
            theta = [coverage(config, N)]
            for n in range(3 * interval):
                config = sweep(V1, V2, mu, N, config, T_start)
                theta.append(coverage(config, N))
            avg_theta = sum(theta) / len(theta)
            Th.append(avg_theta)
        output.M_save(Th, (1 / T_start) * M, T_start)
        
