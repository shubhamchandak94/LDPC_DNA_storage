import numpy as np
import math
import scipy.stats

def main():
	N = 100000
	num_iter = 200
	k = 13
	l = 3
	alpha = 1.0*l/(k-l)
	print(alpha)
	eps = 0.005
	for cov in [2.20,2.22,2.24,2.26]:
		print(cov)
		lamb = cov/(1+alpha)
		P_b = sampled_DE(lamb, num_iter, N, eps, k, l)
		print(P_b)
	
def sampled_DE(lamb, num_iter, N, eps, k, l):
	# returns num_iter array containing estimated bit error rates after t iterations, t = 1,...,num_iter
	P_hat = np.zeros(N)
	P = np.zeros(N)	
	
	P_b = np.zeros(num_iter)

	for t in range(num_iter):
		# var to check nodes
		total_counts = scipy.stats.poisson.rvs(lamb, size=N)
		zero_counts = scipy.stats.binom.rvs(total_counts, 1-eps)
		l_Y = (2*zero_counts-total_counts)*np.log((1-eps)/eps)
		nbr = np.random.randint(N,size=(N,l-1),dtype=int)
		P = l_Y + np.sum(P_hat[nbr],axis=1)

		#check to var nodes
		nbr = np.random.randint(N,size=(N,k-1),dtype=int)
		P_hat = 2*np.arctanh(np.prod(np.tanh(P[nbr]/2),axis=1))
		
		
		#probability error
		total_counts = scipy.stats.poisson.rvs(lamb, size=N)
		zero_counts = scipy.stats.binom.rvs(total_counts, 1-eps)
		l_Y = (2*zero_counts-total_counts)*np.log((1-eps)/eps)
		nbr = np.random.randint(N,size=(N,l),dtype=int)
		P = l_Y + np.sum(P_hat[nbr],axis=1)

		P_b[t] = np.mean((P<=0))

	return P_b   

if __name__ == '__main__':
	main()
