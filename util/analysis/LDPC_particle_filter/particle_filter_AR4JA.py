import numpy as np
import math
import scipy.stats

def main():
	N = 10000
	num_iter = 100
	alpha = 1.0
	eps = 0.00001
	for cov in [1.5,1.6,1.7,1.8]:#[1.6,1.62,1.64,1.66,1.68,1.7]:
		print(cov)
		lamb = cov/(1+alpha)
		P_b = sampled_DE(lamb, num_iter, N, eps)
		print(P_b)
	
def sampled_DE(lamb, num_iter, N, eps):
	# returns num_iter array containing estimated bit error rates after t iterations, t = 1,...,num_iter
	
	P_hat_a_1 = np.zeros(N)
	P_hat_a_5 = np.zeros(N)
	P_hat_b_2 = np.zeros(N)
	P_hat_b_3 = np.zeros(N)
	P_hat_b_4 = np.zeros(N)
	P_hat_b_5 = np.zeros(N)
	P_hat_c_2 = np.zeros(N)
	P_hat_c_3 = np.zeros(N)
	P_hat_c_4 = np.zeros(N)
	P_hat_c_5 = np.zeros(N)

	P_1_a = np.zeros(N)
	P_2_b = np.zeros(N)
	P_2_c = np.zeros(N)
	P_5_a = np.zeros(N)
	P_5_b = np.zeros(N)
	P_5_c = np.zeros(N)
	P_3_b = np.zeros(N)
	P_3_c = np.zeros(N)
	P_4_b = np.zeros(N)
	P_4_c = np.zeros(N)
	
	P_b = np.zeros(num_iter)

	for t in range(num_iter):
		# var to check nodes
		l_Y_5 = np.zeros(N) # punctured variable node
		total_counts = scipy.stats.poisson.rvs(lamb, size=4*N)
		zero_counts = scipy.stats.binom.rvs(total_counts, 1-eps)
		l_Y_rest = (2*zero_counts-total_counts)*np.log((1-eps)/eps)
		l_Y_1 = l_Y_rest[0:N]
		l_Y_2 = l_Y_rest[N:2*N]
		l_Y_3 = l_Y_rest[2*N:3*N]
		l_Y_4 = l_Y_rest[3*N:4*N]
		
	
		# 1
		P_1_a = l_Y_1

		# 2
		nbr_b = np.random.randint(N,size=(N,2),dtype=int)
		nbr_c = np.random.randint(N,size=(N,1),dtype=int)
		sum_b_1 = P_hat_b_2[nbr_b[:,0]]
		sum_b = np.sum(P_hat_b_2[nbr_b],axis=1)
		sum_c = np.sum(P_hat_c_2[nbr_c],axis=1)	
		P_2_b = l_Y_2 + sum_b_1 + sum_c 
		P_2_c = l_Y_2 + sum_b

		# 3
		nbr_b = np.random.randint(N,size=(N,2),dtype=int)
		nbr_c = np.random.randint(N,size=(N,1),dtype=int)
		sum_b_1 = P_hat_b_3[nbr_b[:,0]]
		sum_b = np.sum(P_hat_b_3[nbr_b],axis=1)
		sum_c = np.sum(P_hat_c_3[nbr_c],axis=1)	
		P_3_b = l_Y_3 + sum_b_1 + sum_c
		P_3_c = l_Y_3 + sum_b
		
		# 4
		nbr_b = np.random.randint(N,size=(N,1),dtype=int)
		nbr_c = np.random.randint(N,size=(N,1),dtype=int)
		sum_b = np.sum(P_hat_b_4[nbr_b],axis=1)
		sum_c = np.sum(P_hat_c_4[nbr_c],axis=1)	
		P_4_b = l_Y_4 + sum_c
		P_4_c = l_Y_4 + sum_b
	
		# 5
		nbr_a = np.random.randint(N,size=(N,2),dtype=int)
		nbr_b = np.random.randint(N,size=(N,1),dtype=int)
		nbr_c = np.random.randint(N,size=(N,3),dtype=int)
		sum_a_1 = P_hat_a_5[nbr_a[:,0]]
		sum_a = np.sum(P_hat_a_5[nbr_a],axis=1)
		sum_b = np.sum(P_hat_b_5[nbr_b],axis=1)	
		sum_c_2 = np.sum(P_hat_c_5[nbr_c[:,:2]],axis=1)
		sum_c = np.sum(P_hat_c_5[nbr_c],axis=1)
		P_5_a = l_Y_5 + sum_a_1 + sum_b + sum_c
		P_5_b = l_Y_5 + sum_a + sum_c
		P_5_c = l_Y_5 + sum_a + sum_b + sum_c_2
		#check to var nodes

		# a
		nbr_1 = np.random.randint(N,size=(N,1),dtype=int)
		nbr_5 = np.random.randint(N,size=(N,2),dtype=int)
		prod_1 = np.prod(np.tanh(P_1_a[nbr_1]/2),axis=1)
		prod_5_1 = np.tanh(P_5_a[nbr_5[:,0]]/2)
		prod_5 = np.prod(np.tanh(P_5_a[nbr_5]/2),axis=1)
		
		P_hat_a_1 = 2*np.arctanh(prod_5)
		P_hat_a_5 = 2*np.arctanh(prod_1*prod_5_1)

		# b
		nbr_2 = np.random.randint(N,size=(N,2),dtype=int)
		nbr_3 = np.random.randint(N,size=(N,2),dtype=int)	
		nbr_4 = np.random.randint(N,size=(N,1),dtype=int)
		nbr_5 = np.random.randint(N,size=(N,1),dtype=int)
		prod_2_1 = np.tanh(P_2_b[nbr_2[:,0]]/2)
		prod_2 = np.prod(np.tanh(P_2_b[nbr_2]/2),axis=1)
		prod_3_1 = np.tanh(P_3_b[nbr_3[:,0]]/2)
		prod_3 = np.prod(np.tanh(P_3_b[nbr_3]/2),axis=1)
		prod_4 = np.prod(np.tanh(P_4_b[nbr_4]/2),axis=1)
		prod_5 = np.prod(np.tanh(P_5_b[nbr_5]/2),axis=1)
	
		P_hat_b_2 = 2*np.arctanh(prod_2_1*prod_3*prod_4*prod_5)
		P_hat_b_3 = 2*np.arctanh(prod_2*prod_3_1*prod_4*prod_5)
		P_hat_b_4 = 2*np.arctanh(prod_2*prod_3*prod_5)
		P_hat_b_5 = 2*np.arctanh(prod_2*prod_3*prod_4)

		# c
		nbr_2 = np.random.randint(N,size=(N,1),dtype=int)
		nbr_3 = np.random.randint(N,size=(N,1),dtype=int)	
		nbr_4 = np.random.randint(N,size=(N,1),dtype=int)
		nbr_5 = np.random.randint(N,size=(N,3),dtype=int)
		prod_2 = np.prod(np.tanh(P_2_c[nbr_2]/2),axis=1)
		prod_3 = np.prod(np.tanh(P_3_c[nbr_3]/2),axis=1)
		prod_4 = np.prod(np.tanh(P_4_c[nbr_4]/2),axis=1)
		prod_5_2 = np.prod(np.tanh(P_5_c[nbr_5[:,:2]]/2),axis=1)
		prod_5 = np.prod(np.tanh(P_5_c[nbr_5]/2),axis=1)
	
		P_hat_c_2 = 2*np.arctanh(prod_3*prod_4*prod_5)
		P_hat_c_3 = 2*np.arctanh(prod_2*prod_4*prod_5)
		P_hat_c_4 = 2*np.arctanh(prod_2*prod_3*prod_5)
		P_hat_c_5 = 2*np.arctanh(prod_2*prod_3*prod_4*prod_5_2)
		
		#probability error
		l_Y_5 = np.zeros(N) # punctured variable node
		total_counts = scipy.stats.poisson.rvs(lamb, size=4*N)
		zero_counts = scipy.stats.binom.rvs(total_counts, 1-eps)
		l_Y_rest = (2*zero_counts-total_counts)*np.log((1-eps)/eps)
		l_Y_1 = l_Y_rest[0:N]
		l_Y_2 = l_Y_rest[N:2*N]
		l_Y_3 = l_Y_rest[2*N:3*N]
		l_Y_4 = l_Y_rest[3*N:4*N]
		
	
		# 1
		nbr_a = np.random.randint(N,size=(N,1),dtype=int)
		P_1 = l_Y_1 + np.sum(P_hat_a_1[nbr_a],axis=1)

		# 2
		nbr_b = np.random.randint(N,size=(N,2),dtype=int)
		nbr_c = np.random.randint(N,size=(N,1),dtype=int)
		sum_b = np.sum(P_hat_b_2[nbr_b],axis=1)
		sum_c = np.sum(P_hat_c_2[nbr_c],axis=1)	
		P_2 = l_Y_2 + sum_b + sum_c

		# 3
		nbr_b = np.random.randint(N,size=(N,2),dtype=int)
		nbr_c = np.random.randint(N,size=(N,1),dtype=int)
		sum_b = np.sum(P_hat_b_3[nbr_b],axis=1)
		sum_c = np.sum(P_hat_c_3[nbr_c],axis=1)	
		P_3 = l_Y_3 + sum_b + sum_c
		
		# 4
		nbr_b = np.random.randint(N,size=(N,1),dtype=int)
		nbr_c = np.random.randint(N,size=(N,1),dtype=int)
		sum_b = np.sum(P_hat_b_4[nbr_b],axis=1)
		sum_c = np.sum(P_hat_c_4[nbr_c],axis=1)	
		P_4 = l_Y_4 + sum_b + sum_c
	
		# 5
		nbr_a = np.random.randint(N,size=(N,2),dtype=int)
		nbr_b = np.random.randint(N,size=(N,1),dtype=int)
		nbr_c = np.random.randint(N,size=(N,3),dtype=int)
		sum_a = np.sum(P_hat_a_5[nbr_a],axis=1)
		sum_b = np.sum(P_hat_b_5[nbr_b],axis=1)	
		sum_c = np.sum(P_hat_c_5[nbr_c],axis=1)
		P_5 = l_Y_5 + sum_a + sum_b + sum_c
		P_b[t] = np.mean((P_1<=0) + (P_2<=0) + (P_3<=0) + (P_4<=0))/4
		# not including 5 which is not transmitted
		# should not matter 

	return P_b   

if __name__ == '__main__':
	main()
