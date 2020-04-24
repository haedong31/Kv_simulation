import numpy as np
import matlab.engine
import matplotlib.pyplot as plt


## custom functions -----
def init_pop(N, add_x0, mul_x0, drop_rate, step_size=1.0):
    population = list()

    for _ in range(N):
        # randomly select parameters to drop off
        drop_add = np.ranom.choice([0, 1], size=len(add_x0), p=[drop_rate, 1-drop_rate])
        drop_mul = np.random.choice([0, 1], size=len(mul_x0), p=[drop_rate, 1-drop_rate])

        add_chrom = list()
        for n, a in enumerate(drop_add):
            if a == 0:
                add_chrom.append(0)
            else:
                add_chrom.append(add_x0[n] + np.random.normal(scale=step_size))

        mul_chrom = list()
        for m in enumerate(drop_mul):
            if m == 0:
                mul_chrom.append(1)
            else:
                mul_chrom.append(mul_x0[n] + np.random.normal(scale=step_size))
        
        population.append(add_chrom + mul_chrom)
    
    return population


def eval_fn(Y_hat, Y0, chroms):
    if len(X) != len(Y):
        EnvironmentError('lengths of X and Y are not matching')
    
    fits = list()
    for i in range(len(chroms)):
        fit_i = list()

        for j in range(len(X)):
            fit = (Y[j] - boltzmann(X[j], chroms[i][0], chroms[i][1]))**2
            fit_i.append(fit)

        fits.append(sum(fit_i))

    return fits


def evolve(chroms, fits, N1, N2):
    next_gens = list()
    
    # find superior chromosomes
    super_idx = np.argpartition(fits, N1)[:N1]
    super_chroms = [chroms[idx] for idx in super_idx]
    
    # mutation
    for chrom in super_chroms:
        # clone of parent
        next_gens.append(chrom)
        
        # mutated children
        for _ in range(N2):
            mutant = chrom + np.random.normal(size=2)
            next_gens.append(mutant)
    
    return next_gens


## import data -----
holding_p = -70
holding_t = 4.5 
P1 = 50
P1_t = 29.5 
P2 = 50
P2_t = 29.5 
X0 = matlab.double([0.00000481333, 26.5, 26.5, 0.0000953333, 26.5])

eng = matlab.engine.start_matlab()
outputs = eng.IKslow2(holding_p, holding_t, P1, P1_t, P2, P2_t, X0, nargout=4)
t = np.array(outputs[0]._data)
A = np.array(outputs[2]._data).reshape(outputs[2].size[::-1]).T
IKslow20 = A[:,63]


## main function -----
# hyper parameters
iters = 10000
N0 = 100
N1= 10
N2 = 9
chroms = [None] * (iters+1)
errs = [None] * (iters+1)

# initialization
chroms0 = init_pop(N0)
chroms[0] = chroms0

# run GA
for i in range(iters):
    fits = eval_fn(voltages, ssa, chroms[i])
    errs[i] = min(fits)
    next_gen = evolve(chroms[i], fits, N1, N2)
    chroms[i+1] = next_gen
    
    print('### run %d' % i)

fits = eval_fn(voltages, ssa, chroms[iters])
errs[iters] = min(fits)

last = chroms[iters][np.argmin(fits)]
fitted_ssa = [boltzmann(v, last[0], last[1]) for v in range(-100, -20, 1)]

plt.figure(0)
plt.plot(range(iters+1), errs)
plt.xlabel('Iterations')
plt.ylabel('Sum of squres error')

plt.figure(1)
plt.plot(range(-100, -20, 1), fitted_ssa)
plt.plot(voltages, ssa, 'bo')
plt.xlabel('mV')
plt.ylabel('S-S Activation')
