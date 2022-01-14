import matlab.engine
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## MATLAB model arguments -----
holding_p = -70
holding_t = 4.5 
P1 = 50
P1_t = 29.5
P2 = 50
P2_t = P1_t
eng = matlab.engine.start_matlab()


## custom functions -----
def init_pop(N, X0, step_size=1.0):
    """generate the initial population

    Arguments:
        N {int} -- [population size]
        X0 {list} -- [initial parameter values]

    Keyword Arguments:
        step_size {float} -- [sigma for normal noise] (default: {1.0})

    Returns:
        population {list} -- [N chromosomes with normal noise]
    """

    population = list()
    for _ in range(N):
        chrom = list()

        # add normal noise
        for x0 in X0:
            chrom.append(np.random.normal(loc=x0, scale=step_size))

        population.append(chrom)        
    
    return population


def fitness(y, idx, chroms):
    """evalueate chromosomes with SSE

    Arguments:
        y {numpy array} -- [experimental values]
        idx {int} -- index
        chroms {list} -- chromosomes, each chromosome represents a parameter set

    Returns:
        fits {list} -- evaluation values
    """
    fits = list()
    for i in range(len(chroms)):
        X = matlab.double(chroms[i])
        outputs = eng.IKslow1(holding_p, holding_t, P1, P1_t, P2, P2_t, X, nargout=4)
        A = np.array(outputs[2]._data).reshape(outputs[2].size[::-1]).T
        y_hat = np.max(A[:, idx])
        
        sse = np.sum((y-y_hat)**2)            
        fits.append(sse)

    return fits


def elite_idx(fits, num_elites):
    return np.argpartition(fits, -num_elites)[-num_elites:]


def cross_over(cross_over_num, chroms, fits, rate=0.5):
    # number of parameters 
    num_vars = len(chroms[0])

    # selection probabilities
    probs = fits/np.nansum(fits)
    probs = np.nan_to_num(probs)

    new_gen = list()
    for _ in range(cross_over_num):
        # select parents
        parents_idx = np.random.choice(len(chroms), size=2, replace=False, p=probs)
        parents = np.take(chroms, parents_idx, axis=0)
        
        # set cross-over points
        cross_over_pts = [u >= rate for u in np.random.uniform(0, 1, num_vars)]

        # cross over
        chrom1 = list()
        chrom2 = list()
        for i in range(len(cross_over_pts)):
            if cross_over_pts[i] == True:
                chrom1.append(parents[1][i])
                chrom2.append(parents[0][i])
            else:
                chrom1.append(parents[0][i])
                chrom2.append(parents[1][i])

        new_gen.append(chrom1)
        new_gen.append(chrom2)

    return new_gen


def mutate(chroms, rate, step_size=1):
    mutate_idx = np.random.choice([True, False], size=np.shape(chroms), p=[rate, 1-rate])
    mutant = chroms
    
    if True in mutate_idx:
        chrom_idx, var_idx = np.where(mutate_idx)
        for i in range(len(chrom_idx)):
            mutant[chrom_idx[i]][var_idx[i]] += np.random.normal(scale=step_size)
    else:
        pass

    return mutant


## main function -----
# initialization
data = pd.read_excel('./potassium-KO.xlsx')
y = np.array(data['A2 FF'])
X0 = [22.5000, 2.05800, 45.2000, 1200.00, 45.2000, 0.493000, 170.000]
idx = 64  # IKslow1
burnin = 30
iters = 70
pop_size = 100
cross_over_num = int((pop_size-2)/2)
errs = list()

chroms = init_pop(pop_size, X0, step_size=1.0)
# run GA
for i in range(iters):
    begin_time = time.time()

    fits = fitness(y, idx, chroms)
    errs.append(np.min(fits))
    
    elites = np.take(chroms, elite_idx(fits, num_elites=2), axis=0)
    cross_over_gen = cross_over(cross_over_num, chroms, fits)
    mutate_gen = mutate(cross_over_gen, 0.01, step_size=1)
    chroms = np.vstack((elites, mutate_gen)).tolist()
    
    end_time = time.time()
    print(f"### Iter {i+1}/{iters} | Running time {end_time-begin_time}")

fits = fitness(y, idx, chroms)


outputs3 = eng.IKslow1(holding_p, holding_t, P1, P1_t, P2, P2_t, matlab.double(chroms[0]), nargout=4)
A3 = np.array(outputs3[2]._data).reshape(outputs3[2].size[::-1]).T

testX = chroms[53]
testX[5] = -testX[5]
testX[6] = -testX[6]
testX[6] = 1
outputs4 = eng.IKslow1(holding_p, holding_t, P1, P1_t, P2, P2_t, matlab.double(testX), nargout=4)
A4 = np.array(outputs4[2]._data).reshape(outputs4[2].size[::-1]).T
plt.plot(A4[:,idx])



