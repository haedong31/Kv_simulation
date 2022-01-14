import matlab.engine
import numpy as np
import time
import matplotlib.pyplot as plt


## global variables -----
holding_p = -70
holding_t = 4.5 
P1 = 50
P1_t = 29.5
P2 = 50
P2_t = P1_t
X0 = [22.5000, 2.05800, 45.2000, 1200.00, 45.2000, 0.493000, 170.000]
idx = 64

eng = matlab.engine.start_matlab()
outputs = eng.IKslow1(holding_p, holding_t, P1, P1_t, P2, P2_t, matlab.double(X0), nargout=4)
A = np.array(outputs[2]._data).reshape(outputs[2].size[::-1]).T
y0 = np.max(A[:, idx])
# plt.plot(A[:, idx])

iters = 100
pop_size = 100
num_elites = 2
chroms = [None] * (iters+1)
effs = [None] * (iters+1)


## custom functions -----
def fitness(y0, chroms, idx):
    # |y_hat(i) - y0|
    
    fits = list()
    for i in range(len(chroms)):
        X = matlab.double(chroms[i])
        
        outputs = eng.IKslow1(holding_p, holding_t, P1, P1_t, P2, P2_t, X, nargout=4)
        A = np.array(outputs[2]._data).reshape(outputs[2].size[::-1]).T
        lnA = np.log(A[:, idx])
        y_hat = np.max(lnA)
        
        num_sig = sum([1 for elt in chroms[i] if elt == 0.0 or elt == 1.0])
        if num_sig == 0:
            num_sig = 1
            
        fits.append(np.abs(y_hat - y0)/num_sig)

    return fits


def init_pop(N, add_x0, mul_x0, drop_rate, step_size=1.0):
    population = list()

    for _ in range(N):
        # randomly select parameters to drop off
        drop_add = np.random.choice([0, 1], size=len(add_x0), p=[drop_rate, 1-drop_rate])
        drop_mul = np.random.choice([0, 1], size=len(mul_x0), p=[drop_rate, 1-drop_rate])

        # addictive terms
        add_chrom = list()
        for i, a in enumerate(drop_add):
            if a == 0:
                add_chrom.append(0)
            else:
                add_chrom.append(add_x0[i])

        # multiplicative terms
        mul_chrom = list()
        for j, m in enumerate(drop_mul):
            if m == 0:
                mul_chrom.append(1)
            else:
                mul_chrom.append(mul_x0[j])
        
        population.append(add_chrom + mul_chrom)
    
    return population


def elite_idx(fits, num_elites):
    return np.argpartition(fits, -num_elites)[-num_elites:]


def cross_over(N, num_elites, candidates, fits, rate=0.5):
    # number of parameters 
    num_vars = len(candidates[0])

    # selection probabilities
    probs = fits/sum(fits)

    cross_over_num = int((N-num_elites)/2)
    new_gen = list()
    for _ in range(cross_over_num):
        # select parents
        parents_idx = np.random.choice(len(candidates), size=2, replace=False, p=probs)
        parents = np.take(candidates, parents_idx, axis=0)
        
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


def mutate(candidates, rate, x0, pt):
    # pt: break point of additive and multiplicative parameters (index of the last additive term)
    
    mutate_idx = np.random.choice([True, False], size=np.shape(candidates), p=[rate, 1-rate])
    mutant = candidates

    if True in mutate_idx:
        chrom_idx, var_idx = np.where(mutate_idx)
        for i in range(len(chrom_idx)):
            if var_idx[i] <= pt:
                if mutant[chrom_idx[i]][var_idx[i]] == 0:
                    mutant[chrom_idx[i]][var_idx[i]] = x0[var_idx[i]]
                else:
                    mutant[chrom_idx[i]][var_idx[i]] = 0
            else:
                if mutant[chrom_idx[i]][var_idx[i]] == 1:
                    mutant[chrom_idx[i]][var_idx[i]] = x0[var_idx[i]]
                else:
                    mutant[chrom_idx[i]][var_idx[i]] = 1
    else:
        pass

    return mutant


## main function -----
# initialization
chroms = init_pop(pop_size, X0[0:5], X0[5:7], 0.3)

# run GA
for i in range(iters):
    begin_time = time.time()
    
    if (i+1)%10 == 0:
        print(f"### Iter {i+1}/{iters}")
    
    fits = fitness(y0, chroms, idx)
    effs[i] = np.max(fits)

    elites = np.take(chroms, elite_idx(fits, num_elites=2), axis=0)
    cross_over_gen = cross_over(pop_size, num_elites, chroms, fits)
    mutate_gen = mutate(cross_over_gen, 0.01, X0, 4)
    chroms = np.vstack((elites, mutate_gen)).tolist()
    
    end_time = time.time()
    print(end_time - begin_time)

fits = fitness(y0, chroms, idx)
