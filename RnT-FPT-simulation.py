import numpy as np

pj = 0.1      # probability of jump
pl = 0.8      # probability of laziness
nt = 10**6    # time step (set large to allow time to hit the target)
nh = 4*10**5  # No. of stochastic realizations
xtgt = 4      # target site
x0 = 0        # initial condition
pr = 1. - pj  # probability of remain

file1 = open("RnT-FPT-nT4-n00-pj0-1-pl0-8-nt6-nh4e5.dat","w",1)

# ----------------------------------
# RTP first passage from simulation
# ----------------------------------
xarr = np.zeros(nh,dtype=int)
tarr = np.zeros(nh,dtype=int)
statearr = np.zeros(nh,dtype=int)
np.random.seed(131)
for ih in range(nh):
    x = x0
    state = 1
    for it in range(1,nt):
        if (state > 0):
            coin = np.random.uniform(0.,1.,1)  
            if (1.-pj*pl) >= coin:
                x += 1
            else:
                x -= 1
                state = -1
        elif (state < 0):
            coin = np.random.uniform(0.,1.,1)  
            if pr*pl >= coin:
                x -= 1
            else:
                x += 1
                state = 1
        if x == xtgt:
            xarr[ih] = x
            tarr[ih] = it
            statearr[ih] = state
            break 
        else:
            continue

# Check how many times target not achieved
file1.write('# Total no. trials = '+str(nh)+'\n')
file1.write('# Max waiting time = '+str(nt)+'\n')
file1.write('# Trgt not reached = '+str(np.count_nonzero(xarr != xtgt))+'\n')

# Remove those cases from consideration
itemindex, = np.where(xarr != xtgt)
xarr = np.delete(xarr, itemindex)
tarr = np.delete(tarr, itemindex)
statearr = np.delete(statearr, itemindex)

# Histogram of fpt data for achieved targets
tarr_dim = nh-len(itemindex)
tmax = int(np.max(tarr))
count = np.zeros(tmax)
for i in range(tmax):
    count[i] = np.count_nonzero(tarr == i+1)
    
file1.write('# FPT Norma const. = '+str(np.sum(count)/tarr_dim)+'\n')
file1.write('# Maximum  of  FPT = '+str(tmax)+'\n')
file1.write('# Average  of  FPT = '+str(np.average(tarr))+'\n')
file1.write('# Variance of  FPT = '+str(np.std(tarr)**2.)+'\n')

for i in range(min(tmax,10**5)):
    file1.write('{: 7d} {:.8e}\n'.format(i+1,count[i]/tarr_dim) )

file1.close()