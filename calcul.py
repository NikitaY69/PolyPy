import numpy as np
from random import choice
from copy import copy

class LatticePolymer:
    def __init__(self, N=100, constraint='force', beta_eps=0):
        '''
        Initializes a polymer chain that is to be simulated by a self-avoiding 
        random walk. 
        Initializes the weight according to the choice between an interacting and not interacting random walk.
        Parameters
        ----------
        N : int
            Polymer length
        constraint : string
            Thermo-mecanical constraint applied to the chain. Only 'force' implemented for the moment.
        beta_eps : float
            strength of interacting energy compared to inverse temperature. 
            If beta_eps = 0 there is no closest-non-paired neighbor interaction.
        '''
        # Setting Interaction 
        self.beta_eps = beta_eps
        if self.beta_eps == 0:
            self.interacting = False
        else:
            self.interacting = True
        
        self.N = N
        self.constraint = constraint
        if self.constraint not in ['force', 'length']:
            raise NotImplementedError('Please select constraint in ["force", "length"].')
        
        self.weight = 1
        if self.interacting:
            # This initialization is necessary because the first iteration will wrongly count -1 occupied neighboring sites,
            # The weight has to be balanced in accordance.
            self.weight = np.exp(-self.beta_eps)
        self.weights = [1]
        # there is probably a more clever way to keep track of the sequential weights
        # probably just dividing the global weight by some configurational weight

    def gen_walk(self):
        '''
        Generates a chain of random steps to simulate the polymer. It starts at the center of the grid.
        '''
        # Positioning the initial monomer
        self.pos = [[0, 0, 0]]

        # Looping on the walk
        for step in range(1, self.N):
            self.update_weight()
            if self.number_neighbors() == 0:
                # Stoping the walk when it reaches a closed-loop of neighbors
                break
            # Generating a new direction
            x, y, z = self.random_step()
            while [x, y, z] in self.pos:
                # Generating new step if the step is already present in the history of steps
                x, y, z = self.random_step()

            self.pos.append([x,y,z])

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
    
    def number_pairs(self):
        '''
        This function computes the number of occupied neighbors at a given site in a history of steps. 
        It is relevent to use in the case where a polymer chain can no longer be extended (the final step is surrounded by 6 neighbors),
        or to calculate the weight.
        '''
        neighbors = self.neighborhood(self.pos[-1])
        # Counting the number of occupied neighbors
        c = -1  
        # Since we will certainly count the occupied neighbor from which we just moved, we start the count at -1. 
        # An adjustment was added in the __init__ function to account for the difference in treatment for the very first iteration,
        # where there is no occupied neighbor that we shouldn't count.
        for neighbor in neighbors:
            if neighbor in self.pos:
                c += 1
        return c
    
    def random_step(self):
        '''
        This function generates a random step starting from last visited site.
        '''
        x, y, z = choice(self.neighborhood(self.pos[-1]))
        return x, y, z
    
    def update_weight(self):
        '''
        Updates weight according to the chosen random walk pattern
        '''
        if not self.interacting:
            self.weight *= self.number_neighbors()
            self.weights.append(self.weight)
        else: 
            self.weight *= self.number_neighbors()*np.exp(-self.beta_eps*self.number_pairs())
            self.weights.append(self.weight)

    def length(self):
        '''
        Computes the squared length of a polymer (squared norm between beginning and end of said polymer).
        '''
        return np.linalg.norm(self.pos[-1]-self.pos[0], 2)**2
    
    def gyration(self):
        '''
        Computes the radius of gyration (how "curled up" the polymer is).
        '''
        N = self.pos.shape[0]                       # this in case the polymer is stuck before reaching self.N monomers
        rCM = np.sum(self.pos, axis=0)/N
        return np.sum(np.linalg.norm(self.pos - rCM, ord=2, axis=1)**2)/N
    
    @staticmethod
    def neighborhood(r):
        '''
        Checks the neighboring sites of a given vector r.
        '''
        x, y, z = r
        neighbors = [[(x+1), y, z], [(x-1), y, z], [x, (y+1), z], [x, (y-1), z], [x, y, (z+1)], [x, y, (z-1)]]
        return neighbors


class MonteCarlo(LatticePolymer):
    '''
    Generates collection of polymers.
    Returns thermodynamic observables.
    '''
    def __init__(self, n=10, N=100, constraint='force', beta_eps=0):
        '''
        Parameters
        ----------
        n : int
          Number of monte carlo steps (number of generated polymers)
        '''
        self.n = n
        LatticePolymer.__init__(self, N, constraint, beta_eps)
        self.history = np.empty(shape=self.n, dtype=MonteCarlo)

    def rosenbluth(self):
        '''
        Fills the history with the polymers simulated by a random walk with the normal Rosenbluth method.
        '''
        for trial in range(self.n):
            poly = copy(self)
            poly.gen_walk()
            self.history[trial] = poly
    
    def compute_re(self):
        '''
        Computes the weighted average squared norm of a group of polymers
        '''
        # Weighted averaged length
        return np.average([self.history[trial].length() for trial in range(self.n)], \
                          weights=[self.history[trial].weight for trial in range(self.n)])
    
    def estimate_Z(self, L):
        '''
        This function estimates a partition function for sized-L polymers.
        '''
        W = [self.history[trial].weights[L] for trial in range(self.n)] # weights of all the generated sized-L polymers
        return np.average(W)