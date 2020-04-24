import csv
import numpy as np
import matplotlib.pyplot as plt


## custom functions -----
def boltzmann(v, para1, para2):
    return 1/(1 + np.exp(-(v-para1)/para2))


def init_pop(N):
    population = list()
    for _ in range(N):
        population.append(np.random.normal(size=2))

    return population


def eval_fn(X, Y, chroms):
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
ssa = list()
voltages = list()
num_pts = 0
with open('./Nav_WT_ggmax.csv', mode='rt') as csv_file:
    csv_reader = csv.reader(csv_file)
    next(csv_reader)  # skip the header

    for nth, row in enumerate(csv_reader):
        voltages.append(float(row[0]))
        ssa.append(float(row[1]))
        num_pts += 1
        print('### reading data: %dth data point: %s' % (num_pts, row))


## main function -----
# code parameters
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
