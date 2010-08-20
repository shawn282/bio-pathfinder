from pylab import *

replication_time = 30 # the number of time units it takes a ribosome to dublicate
R0 = 10 # the amount of free ribosomes (non-replicating) at time = 0
total_time = 5000
P_reps = 0.5 ** array(range(1, 20)) # the probability for that a ribosome will start replicating at each time point
R_table = zeros((total_time, len(P_reps)))
growth_rate = zeros((len(P_reps), 1))

for i in range(len(P_reps)):
	replication_queue = [0] * replication_time
	R = R0
	P_rep = P_reps[i]
	for t in range(total_time):
		R_finished = replication_queue.pop(0) * 2 # finished replication
		R_starting = R * P_rep # starting replication
		replication_queue.append(R_starting)
		R = R + R_finished - R_starting
		R_table[t, i] = R
	growth_rate[i, 0] = log(R_table[total_time-1, i] / R_table[total_time-1-replication_time, i])
	
#semilogy(R_table)
loglog(P_reps, growth_rate)
xlabel("P[replication]")
ylabel("Growth rate")
show()
