import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
from matplotlib import ticker

def config_save(X, Y, S, name, interval):
    plt.rc('font', family = 'serif', size = 16)
    plt.rc('figure', figsize = (10.0,10.0))
    fig, ax_scatter = plt.subplots()
    ax_scatter.scatter(X, Y, c = S, s = 80 * np.ones(X.shape), cmap = 'summer')
    #ax_scatter.set_xlabel(r'this is x label 2')
    #ax_scatter.set_ylabel(r'this is y label 2')
    ax_scatter.axis('equal')
    title = "config at " + str(name) + " * " + str(interval) +"sweeps"
    ax_scatter.set_title(title, fontsize = 32, fontweight='bold')
    img = "config" + str(name) + ".jpg"
    plt.savefig(img)
    plt.cla()
    plt.close(fig)


def energy_save(energy, interval, T):
    t = np.linspace(0, interval * (len(energy) - 1), len(energy))
    plt.rc('font', family = 'serif', size = 16)
    plt.rc('figure', figsize = (10.0,10.0))
    fig, ax = plt.subplots()
    ax.plot(t, energy, c = 'tomato')
    ax.set_xlabel(r'steps')
    ax.set_ylabel(r'regularized energy')
    title = "Energy at " + "kBT = " + str(T)
    ax.set_title(title, fontsize = 32, fontweight='bold')
    img = "energy" + str(T) + "K.jpg"
    plt.savefig(img)
    plt.cla()

def theta_save(theta, interval, T):
    t = np.linspace(0, interval * (len(theta) - 1), len(theta))
    plt.rc('font', family = 'serif', size = 16)
    plt.rc('figure', figsize = (10.0,10.0))
    fig, ax = plt.subplots()
    ax.plot(t, theta, c = 'royalblue')
    ax.set_xlabel(r'steps')
    ax.set_ylabel(r'coverage')
    ax.set_ylim(0, 1)
    title = "Coverage at " + "kBT = " + str(T)
    ax.set_title(title, fontsize = 32, fontweight='bold')
    img = "coverage" + str(T) + "K.jpg"
    plt.savefig(img)
    plt.cla()
    plt.close(fig)

def P_save(X, Y, S, T, mu, Config_init, theta):
    plt.rc('font', family = 'serif', size = 16)
    plt.rc('figure', figsize = (10.0,10.0))
    fig, ax_scatter = plt.subplots()
    ax_scatter.scatter(X, Y, c = S, s = 100 * np.ones(X.shape), cmap = 'summer')
    label = "theta = " + str(theta)
    ax_scatter.set_xlabel(label)
    #ax_scatter.set_ylabel(r'this is y label 2')
    ax_scatter.axis('equal')
    title = "kBT = " + str(T) + ", mu = " + str(mu)
    ax_scatter.set_title(title, fontsize = 32, fontweight='bold')
    img = "config-c" + str(Config_init) + "-T" + str(T) + "-mu" + str(mu + 4) + ".jpg"
    plt.savefig(img)
    plt.cla()
    plt.close(fig)

def M_save(Th, M, T):
    plt.rc('font', family = 'serif', size = 16)
    plt.rc('figure', figsize = (10.0,10.0))
    fig, ax = plt.subplots()
    ax.plot(M, Th, c = 'forestgreen')
    ax.set_xlabel(r'mu/kBT')
    ax.set_ylabel(r'theta')
    title = "Coverage at kBT = " + str(T) 
    ax.set_title(title, fontsize = 32, fontweight='bold')
    img = "mu" + str(T) + ".jpg"
    plt.savefig(img)
    plt.cla()
    plt.close(fig)
