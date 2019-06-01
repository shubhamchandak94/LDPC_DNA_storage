import numpy as np
import math
import scipy.stats

def main():
	N = 100000
	num_iter = 200
	k = 33
	eps_bar = 0.09	
	lambda_list = find_poly_lambda(eps_bar,k)
	rate = 1 - 1.0/(k*np.sum(lambda_list/np.arange(2,np.size(lambda_list)+2)))
	print('rate',rate)
	alpha = 1.0/rate - 1
	print('alpha',alpha)
	eps = 0.005
	for cov in [2.8,2.9,3.0]:
		print(cov)
		lamb = cov/(1+alpha)
		P_b = sampled_DE(lamb, num_iter, N, eps, k, lambda_list)
		print(P_b)
	
def sampled_DE(lamb, num_iter, N, eps, k, lambda_list):
	# returns num_iter array containing estimated bit error rates after t iterations, t = 1,...,num_iter
	P_hat = np.zeros(N)
	P = np.zeros(N)	
	
	l_values = np.arange(2,np.size(lambda_list)+2)
	for t in range(num_iter):
		# var to check nodes
		total_counts = scipy.stats.poisson.rvs(lamb, size=N)
		zero_counts = scipy.stats.binom.rvs(total_counts, 1-eps)
		l_Y = (2*zero_counts-total_counts)*np.log((1-eps)/eps)
		var_degrees = np.random.choice(l_values,N,p=lambda_list)
		nbrs = np.random.randint(N,size=np.sum(var_degrees)-N,dtype=int)
		cur_pos = 0
		for i in range(N):
			P[i] = l_Y[i] + np.sum(P_hat[nbrs[cur_pos:cur_pos+var_degrees[i] - 1]])
			cur_pos += var_degrees[i] - 1

		#check to var nodes
		nbr = np.random.randint(N,size=(N,k-1),dtype=int)
		P_hat = 2*np.arctanh(np.prod(np.tanh(P[nbr]/2),axis=1))
		
		
	#probability error
	total_counts = scipy.stats.poisson.rvs(lamb, size=N)
	zero_counts = scipy.stats.binom.rvs(total_counts, 1-eps)
	l_Y = (2*zero_counts-total_counts)*np.log((1-eps)/eps)
	var_degrees = np.random.choice(l_values,N,p=lambda_list)
	nbrs = np.random.randint(N,size=np.sum(var_degrees),dtype=int)
	cur_pos = 0
	for i in range(N):
		P[i] = l_Y[i] + np.sum(P_hat[nbrs[cur_pos:cur_pos+var_degrees[i]]])
		cur_pos += var_degrees[i]

	P_b = np.mean((P<=0))

	return P_b   

def find_poly_lambda(eps_bar,k):
    lambda_list = []
    l = 2
    while(True):
        #lambda_list.append((-1)**l/eps_bar*math.gamma(1+1/(k-1))\
        #                   /(math.gamma(l)*math.gamma(-l+2+1/(k-1))))
        lambda_list.append(abs(1.0/eps_bar*math.exp(math.lgamma(1+1.0/(k-1))-math.lgamma(l)\
                                                       -math.lgamma(-l+2+1.0/(k-1)))))
        if(sum(lambda_list) >= 1):
            break
        l+=1
    lambda_list = np.array(lambda_list)
    lambda_list = lambda_list/np.sum(lambda_list)
    return lambda_list

def variable_node_degree_distribution(lambda_list):
    L_list = np.zeros(np.size(lambda_list))
    for i in range(np.size(lambda_list)):
        L_list[i] = lambda_list[i]/(i+2)
    L_list = L_list/np.sum(L_list)   
    return(L_list)

if __name__ == '__main__':
	main()
