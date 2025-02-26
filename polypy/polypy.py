import numpy as np
from random import choice, choices, uniform
from copy import deepcopy as copy
from tqdm import tqdm
from colorama import Fore, Style
import pickle 
# from itertools import zip_longest

class LatticePolymer:
    def __init__(self, N, boltzmann_energy=1, boltzmann_force=1):
        '''
        Initializes a polymer chain that is to be simulated by a self-avoiding 
        random walk. 
        Initializes the weight according to the choice between an interacting 
        and not interacting random walk.

        Parameters
        ----------
        N : int
            Polymer length
        boltzmann_energy : float
            Boltzmann energy as in exp(-beta*epsilon)
        boltzmann_force : float
            Boltzmann force as in exp(beta*force)
        '''
        # Setting Interaction and constraints
        self.q = boltzmann_energy
        self.b = boltzmann_force

        if self.q == 1:
            self.interacting = False
        else:
            self.interacting = True
        if self.b != 1:
            self.forced = True
        else:
            self.forced = False

        self.N = N
        
        # Number of clones and number of pruned steps
        # self.n_c = 0
        # self.n_p = 0

    def gen_walk(self, start=1, perm=False, c_m=0.2, c_p=2):
        '''
        Generates a chain of random steps to simulate the polymer. 
        It starts at the center of the grid.

        Parameters
        ----------
        start : int
            Number of monomers already present in the polymer
        perm : bool
            Tells if pruning/enriching is applied during sampling
        c_m : float
            Pruning strength (as in lower threshold = c_m * current estimator of Z)
        c_p : float
            Cloning strength (as in upper threshold = c_p * current estimator of Z)
        '''
        # Positioning the initial monomer
        if start == 1:
            self.pos = [[0, 0, 0]]
            self.weight = 0
        
        # self.heatup = 0
        # heatup_thres = max(10, self.N // 500)
        # if self.tours.count(self.tour) >= self.relaxation and self.origin <= int(0.2*self.N):
        #     print('%sRelaxing PERM...%s' % (Fore.YELLOW, Style.RESET_ALL))

        # Number of 
        # self.c0 = len(self.clones)
        # self.k_ = 0
        # self.c_ = 0

        # Looping on the walk
        try:
            for step in tqdm(range(start, self.N)):
                # Stopping the walk when it reaches a closed-loop of neighbors
                if self.number_neighbors() == 0:
                    break
                
                # Updating new weight and partition function
                self.update_weight(step)
            
                # self.heatup+=1
                # if start == 1 and not self.clones:
                #     heatup_thres = self.start_heatup
                # else:
                #     heatup_thres = max(5, self.N//500)
                if perm: # and self.heatup >= heatup_thres:
                        # input('%sAdding FORCING to polymer.%s' % (Fore.YELLOW, Style.RESET_ALL))
                        # c_m /= 10
                    # Applying the pruning/enriching algorithm
                    self.control_weight(step, c_m, c_p)
                
                # Mean number of killed and cloned polymers
                # if len(self.clones) == self.c0:
                #     self.c_ += 1
                # self.k_ += 1
                    
        # If control_weight kills a polymer
        except BreakException:
            pass
        
        # Converting list of positions as numpy.array
        self.pos = np.array(self.pos)

    def number_neighbors(self):
        '''
        This function computes the number of free neighboring sites of the last visited site in a history of steps. 
        It is relevant to use in the case where a polymer chain can no longer be extended 
        (the final step is surrounded by 6 neighbors) and to calculate the weight.
        '''
        neighbors = self.neighborhood(self.pos[-1])
        c = 0
        for neighbor in neighbors:
            if neighbor not in self.pos:
                c += 1
        return c
    
    def random_step(self):
        '''
        This function generates a random step starting from the last visited site.
        '''
        if not self.forced:
            # in this case, we choose uniformly
            x, y, z = choice(self.neighborhood(self.pos[-1]))

        else:
            # in this case, we choose with a bias
            weights = [np.sqrt(self.b), 1/np.sqrt(self.b), 1, 1, 1, 1]
            s = np.sum(weights)
            weights = [w/s for w in weights]
            x, y, z = choices(self.neighborhood(self.pos[-1]), weights=weights)[0]

        return x, y, z
    
    def control_weight(self, step, c_m, c_p):
        '''
        This function applies the pruning/enriching algorithm to the Rosenbluth sampling.

        Parameters
        ----------
        step : int
            Current step in the generation
        c_m : float
            Pruning strength (as in lower threshold = c_m * current estimator of Z)
        c_p : float
            Cloning strength (as in upper threshold = c_p * current estimator of Z)
        '''
        # Current weight thresholds
        W_m = np.log10(c_m)+self.logZ[step]
        W_p = np.log10(c_p)+self.logZ[step]

        # Relaxation
        # if self.tours.count(self.tour) >= self.relaxation and self.origin <= int(0.2*self.N):
        #     W_m = 0
        #     W_p = 10**(10000)

        # self.heatup=0 

        # Pruning
        if self.weight < W_m:

            # Number of clones and pruned steps
            # self.n_p += 1
            # self.n_p += 1

            # Mean number of killed polymers
            # self.k.append(self.k_)
            # self.k_ = 0
            
            if uniform(0, 1) < 0.5:
                print('%sPolymer has been KILLED!%s' % (Fore.RED, Style.RESET_ALL))
                raise BreakException()
            else:
                print('%sPolymer has SURVIVED!%s' % (Fore.GREEN, Style.RESET_ALL))
                self.weight += np.log10(2)

        elif self.weight > W_p and step != self.N-1: # and step <= int(0.95*self.N):

            # Number of cloned steps
            # self.n_c += 1
            
            self.weight -= np.log10(2)
            self.clones.append(self.checkpoint())
            print('%sPolymer has been CLONED!%s' % (Fore.MAGENTA, Style.RESET_ALL))

            # Mean number of cloned polymers 
            # self.c.append(self.c_)
            # self.c0 = len(self.clones)

    
    def checkpoint(self):
        '''
        This function saves the key properties of a polymer at any step of its growth.
        '''
        return {'weight': copy(self.weight), 'pos': copy(self.pos)}
    
    def reset(self, w, p):
        '''
        This function resets the weight and the positions of a polymer to w and p.
        '''
        self.weight = w
        self.pos = p

    def update_weight(self, step):
        '''
        Updates weight according to the chosen random walk pattern.

        Parameters
        ----------
        step : int
            Current step in the generation
        '''
        # Number of free neighbor nearest sites
        numb_neigh = self.number_neighbors()

        # Generating a new direction
        x, y, z = self.random_step()
        while [x, y, z] in self.pos:
            # Generating new step if the step is already present in the history of steps
            x, y, z = self.random_step()
        self.pos.append([x,y,z])

        # Updating the weight
        if not self.forced:
            if not self.interacting:
                self.weight += np.log10(numb_neigh)
            else: 
                # Number of closest non-bound occupied sites
                numb_pairs = 5 - self.number_neighbors()
                self.weight += np.log10(numb_neigh*self.q**(numb_pairs))

        else:
            numb_pairs = 5 - self.number_neighbors()
            # x-travel between last two monomers
            delta_x = self.pos[-1][0] - self.pos[-2][0]
            if delta_x == 0:
                self.weight += np.log10(self.q**(numb_pairs))
            elif delta_x == 1:
                self.weight += np.log10(self.q**(numb_pairs)*np.sqrt(self.b))
            else:
                self.weight += np.log10(self.q**(numb_pairs)/np.sqrt(self.b))

        # Calculation of Z (for Monte Carlo purposes)
        # Finding the max weight
        self.weights[step].append(self.weight)
        w_max = np.array(self.weights[step]).max()
        self.diff_weights = [x-w_max for x in self.weights[step]]

        # a is the list of [a_i] in our project report
        self.a = np.sum(np.power(10, self.diff_weights))
        trials = len(self.weights[step])
        self.logZ[step] = np.log10(1/(trials)) + np.log10(self.a) + w_max

    '''
    In the next three following staticmethods, pos is a numpy.array of shape(x, 3) with x varying.
    '''
    @staticmethod
    def length(pos):
        '''
        Computes the squared end-to-end length of a set of vectors.
        '''
        return np.linalg.norm(pos[-1]-pos[0], 2)**2
    
    @staticmethod
    def gyration(pos):
        '''
        Computes the radius of gyration (how "curled up" the polymer is) of a set of vectors.
        '''
        rCM = np.average(pos, axis=0)
        return np.average(np.linalg.norm(pos - rCM, ord=2, axis=1)**2)
    
    @staticmethod
    def extension(pos):
        '''
        Computes the extension of a set of vector in the x direction.
        '''
        return pos[-1][0] - pos[0][0]
    
    @staticmethod
    def neighborhood(r):
        '''
        Checks the neighboring sites of a given vector r.
        '''
        x, y, z = r
        neighbors = [[(x+1), y, z], [(x-1), y, z], [x, (y+1), z], [x, (y-1), z], [x, y, (z+1)], [x, y, (z-1)]]
        return neighbors

def MonteCarloFactory(load=None, *args, **kwargs):
    '''
    This function serves as a factory to load or initiate a new MonteCarlo instance.

    Parameters
    ----------
    load : str
        path to saved MonteCarlo object
    *args : arguments to pass in MonteCarlo
    **kwargs : keyword arguments to pass in MonteCarlo
    '''
    if load:
        with open(load, "rb") as f:
            return pickle.load(f)
    else:
        return MonteCarlo(*args, **kwargs)
    
class MonteCarlo(LatticePolymer):
    '''
    Generates collection of polymers.
    Computes thermodynamic observables.
    '''
    def __init__(self, N=100, boltzmann_energy=1, boltzmann_force=1):
        '''
        Parameters
        ----------
        Init parameters of LatticePolymer
        '''
        LatticePolymer.__init__(self, N, boltzmann_energy, boltzmann_force)
        self.weights = [[] for _ in range(self.N)]
        self.logZ = np.zeros(shape=self.N)
        self.history = {'weight': [], 'pos': [], 'origin': []}

    def rosenbluth(self, perm=False, save=None, **kwargs):
        '''
        Fills the history with the polymers simulated by a random walk with a Rosenbluth 
        sampling strategy.

        Parameters
        ----------
        perm : bool
            Tells if pruning/enriching is applied during sampling
        save : str
            Path to save MonteCarlo object
        c_m : float (keyword)
            Pruning strength (as in lower threshold = c_m * current estimator of Z)
        c_p : float (keyword)
            Cloning strength (as in upper threshold = c_p * current estimator of Z)
        '''
        # Run variables
        self.clones = []
        # self.clones will eventually take the form of a list of clones properties

        # Run globals
        self.perm = perm
        c_m = kwargs.get('c_m', 0.2)     # lower threshold
        c_p = kwargs.get('c_p', 10*c_m)  # upper threshold

        # Relaxation
        # self.relaxation = kwargs.get('relaxation', max(250, self.n//5))

        # pruning/enriching is only applied after start trials
        start = 2  
        
        # Run iterators
        # self.start_heatup = max(5, self.N//500)
        # self.tours = []    
        self.trial = 0                           # total number of polymers created
        self.desired_trials = 0                  # number of polymers which were of size self.N
        self.cloning_freeze = 0                  # number of trials with no clones generated
        self.tour = 0                            # number of tours                     

        # Mean number of cloned and killed polymers
        # self.c = []
        # self.k = []

        while self.desired_trials < self.n:
            print('Simulating Polymer %d / Trial %d / Tour %d' % (self.desired_trials, self.trial, self.tour))
            
            self.origin = 0   # starting position of polymer growth

            # if no clones are generated during 7500 trials, we start a new run
            if self.cloning_freeze >= 7500:
               break

            # For the first trials, no PERM is applied
            if (self.trial < start and self.r==1) or not self.perm:
                self.gen_walk(perm = False)

            else:
                # Cheking if a clone has been generated for this trial
                if self.clones:
                    clone = self.clones[-1]
                    # Number of monomers already present in current polymer
                    m = len(clone['pos'])                                
                    self.origin = m
                    self.reset(clone['weight'], clone['pos'])

                    # Processing polymer growth on top of the clone
                    self.gen_walk(m, perm=True, c_m=c_m, c_p=c_p)        
                    self.clones.remove(clone)
                    self.cloning_freeze = 0

                # Else generating polymer from scratch
                else:    
                    self.gen_walk(perm=True, c_m=c_m, c_p=c_p)
                    self.cloning_freeze += 1

            # Filling the history
            self.history['weight'].append(self.weight)
            self.history['pos'].append(self.pos) 
            self.history['origin'].append(self.origin)

            # if self.cloning_freeze >= 40:
            #     self.start_heatup += self.N//500
            #     input('%sPress SPACE to increase PERM cooldown to %d steps%s' % (Fore.CYAN, self.start_heatup, Style.RESET_ALL))

            if not self.clones:
                self.tour += 1
            self.trial += 1
            if self.pos.shape[0] == self.N:
                self.desired_trials += 1

            # self.tours.append(self.tour)
                
        if save != None:
            self.save(save)

    def multiple_PERM(self, runs=50, poly_per_run=100, c_m=(1/3), c_p=3, save=None):
        '''
        This function runs a collection of PERM simulations to achieve tour-decorrelation.

        Parameters
        ----------
        runs : int
            Number of decorrelated runs
        poly_per_run : int
            Number of desired polymers generated in each run
        c_m : float 
            Pruning strength (as in lower threshold = c_m * current estimator of Z)
        c_p : float
            Cloning strength (as in upper threshold = c_p * current estimator of Z)
        save : str
            Path to save MonteCarlo object
        '''
        self.n = poly_per_run

        # Variables
        # self.weights_stack = []
        # self.history_stack = []

        # Iterators
        self.r = 1                    # current run
        self.all_tours = 0            # total number of tours generated
        self.all_trials = 0           # total number of trials generated
        self.all_desired_trials = 0   # total number of desired trials generated

        # Sequential decorrelated runs
        while self.r <= runs:
            print('%sStarting run %d/%d%s' % (Fore.YELLOW, self.r, runs, Style.RESET_ALL))
            self.rosenbluth(perm=True, c_m=c_m, c_p=c_p, save=None)

            # self.weights_stack.append([list(i) for i in zip_longest(*self.weights)])
            # self.history_stack.append(self.history)

            self.all_tours += self.tour
            self.all_desired_trials += self.desired_trials
            self.all_trials += self.trial
            self.r += 1
            print('%sSimulated so far %d polymers / %d trials / %d tours%s' % \
                  (Fore.YELLOW, self.all_desired_trials, self.all_trials, \
                   self.all_tours, Style.RESET_ALL))
            

        # Rebuilding stacked MC attributes

        # all_origins = [x['origin'] for x in self.history_stack]
        # all_positions = [x['pos'] for x in self.history_stack]

        # self.history['origin'] = sum(all_origins, [])
        # self.history['pos'] = sum(all_positions, [])
        # # self.weights = list(np.concatenate(self.weights_stack).T)
        # transposed_weights = sum(self.weights_stack, [])
        # self.weights = [list(i) for i in zip_longest(*transposed_weights)]

        # # Computing whole Z
        # self.Z = np.zeros(shape=self.N)

        # # Finding the max weight
        # for step in range(1, self.N):
        #     logweights = [x for x in self.weights[step] if x is not None]
        #     w_max = np.array(logweights).max()
        #     self.y = [x-w_max for x in logweights]

        #     self.zfactor = np.sum(np.power(10, self.y))
        #     trials = len(self.weights[step])
        #     self.Z[step] = np.log10(1/(trials)) + np.log10(self.zfactor) + w_max

        if save != None:
            self.save(save)

    def compute_observable(self, obs, N):
        '''
        Computes an observable average at a given number of monomers.

        Parameters
        ----------
        obs : func
            function taking in argument a set of vectors and returning a desired observable
        N : int
            number of desired monomers
        '''
        logweights = [x for x in self.weights[N-1]] # if x is not None]
        trials = len(logweights)
        # positions of all the polymers in the history of length at least N
        positions = [pos[:N] for i, pos in enumerate(self.history['pos']) \
                     if pos.shape[0] >= N and self.history['origin'][i] < N]

        observables = [obs(pos) for pos in positions]
        # normalized sampling weights (list of [W_i^N] in our project report)
        weights = np.power(10, [w-np.log10(trials)-self.logZ[N-1] for w in logweights])
        
        return np.average(observables, weights=weights)

    def error(self, obs, N):
        '''
        Estimate the Monte Carlo error using approximate analytical formula and boostrapping

        Parameters
        ----------
        obs : func
            function taking in argument a set of vectors and returning a desired observable
        N : int
            number of desired monomers
        '''
        # Approximate analytical formula
        # err_a = np.sqrt((n/(n-1))*(np.sum([self.history[trial].weight**2*(self.history[trial].length()**2- np.sum([self.history[j].length() for j in range(self.n)]))**2) for trial in range(self.n)]/(w**2))
        logweights = [x for x in self.weights[N-1]] # if x is not None]
        trials = len(logweights)
        # positions of all the polymers in the history of length at least N
        positions = [pos[:N] for i, pos in enumerate(self.history['pos']) \
                     if pos.shape[0] >= N and self.history['origin'][i] < N]

        observables = np.array([obs(pos) for pos in positions])
        # normalized sampling weights (list of [W_i^N] in our project report)
        weights = np.power(10, [w-np.log10(trials)-self.logZ[N-1] for w in logweights])

        # Boostrapping
        num_bootstrap_samples = trials//10

        # array of each new sampled observable average
        bootstrap_obs = np.zeros(num_bootstrap_samples)
    
        for i in range(num_bootstrap_samples):
            # Generating bootstrapped samples
            bootstrap_indices = np.random.choice(trials, size=trials, replace=True)
            bootstrap_obs[i] = np.average(observables[bootstrap_indices], weights=weights[bootstrap_indices])

        error = np.std(bootstrap_obs)

        return error
    
    def save(self, file):
        '''
        Saving MonteCarlo object as .pkl file

        Parameters
        ----------
        file : str
            Path to save MonteCarlo object
        '''
        with open(file, "wb") as f:
            pickle.dump(self, f)

class BreakException(Exception):
    pass