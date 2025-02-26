import numpy as np
import sys
from polypy import LatticePolymer, MonteCarloFactory, MonteCarlo
import visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os
import time

# Measuring the computing time
start_time = time.time()

# Plot params
mpl.rcParams['figure.figsize'] = (8,7)
mpl.rcParams['font.size'] = 15
mpl.rcParams['axes.linewidth'] = 3
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']


def make_axis(axe):
    if not log:
        axe.xaxis.set_minor_locator(AutoMinorLocator())
    axe.yaxis.set_minor_locator(AutoMinorLocator())
    axe.xaxis.set_tick_params(which='major', direction='in', size=10, width=2, top='on', bottom='on')
    axe.xaxis.set_tick_params(which='minor', direction='in', size=6, width=1, top='on', bottom='on')
    axe.yaxis.set_tick_params(which='major', direction='in', size=10, width=2, right='on', left='on')
    axe.yaxis.set_tick_params(which='minor', direction='in', size=6, width=1, right='on', left='on')
    pass

# Run params
arg = sys.argv[1]
if arg not in ['sarw', 'isaw', 'bisaw']:
    raise NotImplementedError("Please provide a run type in ['sarw', 'isarw', 'bisaw']")

# Params
N = 100                 # Number of monomers
poly_per_run = 10       
runs = 1
c_m = 0.05              # PERM lower threshold's coefficient
c_p = 6.5               # PERM higher threshold's coefficient
log = True              # True if SARW plotted in logarithmic scale


# SARW simulation
if arg == 'sarw':
    # Generating a group of polymers with PERM method

    if not os.path.exists('%s_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (arg, runs, N, poly_per_run, c_m, c_p)):
        mcgroup = MonteCarloFactory(N=N)
        mcgroup.multiple_PERM(runs=runs, poly_per_run=poly_per_run, c_m=c_m, c_p=c_p, \
                            save='%s_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (arg, runs, N, poly_per_run, c_m, c_p))
    

    # Universal law calculation
    def r(L):
        nu = 0.586
        return L**(2*nu)
    
    # Creating values of N at which we evaluate the observable
    if log :
        ks=np.logspace(1, np.log10(N), num=15)  
        ks = [int(k) for k in ks]
    else :
        ks= [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1]
        ks = [int(N*k) for k in ks]

    x = np.linspace(ks[0], ks[-1], 100)
    ts = []
    es = []

    # Loading already saved simulation if lines 47-49 are commented
    mcgroup = MonteCarloFactory(load='sarw_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (runs, N, poly_per_run, c_m, c_p))

    # Creating results data
    for k in ks:
        t = mcgroup.compute_observable(LatticePolymer.gyration, k)
        e = mcgroup.error(LatticePolymer.gyration, k)
        if log :
            ts.append(t/r(k-1)) 
            es.append(e/r(k-1))
        else :
            ts.append(t)
            es.append(e)
        print("re(%d) = %f +/- %f" % (k, t, e))
        print('r_theo(%d)=' % k, r(k-1))

    plt.errorbar(ks, ts, yerr=es, fmt='ro', mfc='red', mec='k', ms=6.0, capsize=3, lw=0.9, label='PERM')
    if not log : plt.plot(x, r(x), 'b-', lw=2, label='Approximate law')
    if log : 
        plt.xscale('log')
    ax = plt.gca()
    make_axis(ax)
    plt.legend(loc=2, frameon=False)
    if not log:
        plt.ylabel(r'$\langle r_g^2\rangle$')
    else: 
        plt.ylabel(r'$\frac{\langle r_g^2\rangle}{N^{2\nu}}$')
    plt.xlabel(r'$N$')
    plt.savefig('Re_%druns_%dmonom_%dpoly_%.2f_%.2f.jpg' % (runs, N, poly_per_run, c_m, c_p), dpi=300)

    end_time = time.time()
    plt.show()


# ISAW simulation
elif arg == 'isaw':
    # Chosen interaction strengths
    energy = [1.3, 1.305, 1.3087, 1.310, 1.315]

    # Creating values of N at which we evaluate the observable
    ks = np.logspace(1, np.log10(N), num=15)
    ks = [int(k) for k in ks]
    fmts = ['r+-', 'gs-', 'bx-', 'm^-', 'cD-']
    for i, en in enumerate(energy):

        
        if not os.path.exists('%s_%.3fq_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                    (arg, en, runs, N, poly_per_run, c_m, c_p)):

            # Generating a group of polymers with PERM method
            mcgroup = MonteCarloFactory(N=N, boltzmann_energy=en)
            mcgroup.multiple_PERM(runs=runs, poly_per_run=poly_per_run, c_m=c_m, c_p=c_p, \
                                save='%s_%.3fq_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                    (arg, en, runs, N, poly_per_run, c_m, c_p))
        
        ts = np.empty(shape=(len(energy), len(ks)))
        es = np.empty(shape=(len(energy), len(ks)))

        # Loading saved simulations
        mcgroup = MonteCarloFactory(load='%s_%.3fq_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                (arg, en, runs, N, poly_per_run, c_m, c_p))
        
        # Creating results data
        for j, k in enumerate(ks):
            t = mcgroup.compute_observable(LatticePolymer.length, k)
            e = mcgroup.error(LatticePolymer.length, k)
            ts[i,j] = (t/k)
            es[i,j] = (e/k)

        plt.plot(ks, ts[i, :], fmts[i], mfc='none', ms=8.0, lw=0.5, label='q=%.4f' % en)

    # Plotting the data  
    ax = plt.gca()
    make_axis(ax)
    plt.xscale('log')
    plt.legend(loc=2)
    plt.ylabel(r'$\frac{1}{N}\langle r_e^2\rangle$')
    plt.xlabel(r'$N$')
    plt.savefig('Re_%s_%druns_%dmonom_%dpoly_%.2f_%.2f.jpg' % (arg, runs, N, poly_per_run, c_m, c_p), dpi=300)

    end_time = time.time()
    plt.show()

# BISAW simulation
else:
    # Setting the energy
    energy = 1.5

    # Chosen forces
    force = [1, 1.45, 1.55, 1.6, 1.65, 1.70]

    # Creating values of N at which we evaluate the observable
    ks = np.linspace(2, N, num=50)
    ks = [int(k) for k in ks]
    fmts = ['r-', 'g.', 'c--', 'm-', 'y.', 'b-']
    for i, f in enumerate(force):
        if not os.path.exists('%s_%.2fb_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                    (arg, f, runs, N, poly_per_run, c_m, c_p)):
            # Generating a group of polymers with PERM method
            mcgroup = MonteCarloFactory(N=N, boltzmann_energy=energy, boltzmann_force=f)
            mcgroup.multiple_PERM(runs=runs, poly_per_run=poly_per_run, c_m=c_m, c_p=c_p, \
                                save='%s_%.2fb_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                    (arg, f, runs, N, poly_per_run, c_m, c_p))
        
        # Loading saved simulations
        mcgroup = MonteCarloFactory(load='%s_%.2fb_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                (arg, f, runs, N, poly_per_run, c_m, c_p))
        
        # Creating results data
        ts = np.empty(shape=(len(force), len(ks)))
        es = np.empty(shape=(len(force), len(ks)))
        for j, k in enumerate(ks):
           
            t = mcgroup.compute_observable(LatticePolymer.extension, k)
            e = mcgroup.error(LatticePolymer.extension, k)
            ts[i,j] = t
            es[i,j] = e

        plt.plot(ks, ts[i, :], fmts[i], lw=0.5, label='b=%.3f' % f)
        
    # Plotting the data
    ax = plt.gca()
    make_axis(ax)
    plt.legend(loc='upper left')
    plt.ylabel(r'$\langle x\rangle$')
    plt.xlabel(r'$N$')
    plt.savefig('X_%s_%druns_%dmonom_%dpoly_%.2f_%.2f.jpg' % (arg, runs, N, poly_per_run, c_m, c_p), dpi=300)
    end_time = time.time()
    plt.show()



# Heatmap
visualisation.polyCloud3D(mcgroup, N)
    
# Chopped heatmap
# visualisation.polyCloud3Dchop(mcgroup,N)

elapsed_time = end_time - start_time

print("Elapsed time:", elapsed_time, "seconds")