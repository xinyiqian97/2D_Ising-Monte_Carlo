import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
from matplotlib import ticker

def config_save(X, Y, S, name, interval):
    plt.rc('font', family = 'serif', size = 16)
    plt.rc('figure', figsize = (10.0,10.0))
    fig, ax_scatter = plt.subplots()
    ax_scatter.scatter(X, Y, c = S, s = 32 * np.ones(X.shape), cmap = 'bwr')
    #ax_scatter.set_xlabel(r'this is x label 2')
    #ax_scatter.set_ylabel(r'this is y label 2')
    ax_scatter.axis('equal')
    title = "config at " + str(name) + " * " + str(interval) +"sweeps"
    ax_scatter.set_title(title, fontsize = 32, fontweight='bold')
    img = "config" + str(name) + ".jpg"
    plt.savefig(img)


#def energy_save(energy):
