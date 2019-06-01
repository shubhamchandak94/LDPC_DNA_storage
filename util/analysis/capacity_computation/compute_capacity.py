import numpy as np
import random
import argparse
import json
import scipy.stats
import scipy.special

def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--eps', type=float, default=0.0)
    parser.add_argument('--outfile', type=str,required=True)
    return parser

def compute_capacity(eps,outfile):
	alpha = np.linspace(0.01,1.0,100)	
	if eps == 0:
		coverage = (1+alpha)*np.log(1+1/alpha)
		f_out = open(outfile,'w')
		f_out.write("c_w,c_r\n")
		for i in range(np.size(alpha)):
			f_out.write(str("%.2f" % round(1+alpha[i],2)))
			f_out.write(',')
			f_out.write(str("%.2f" % round(coverage[i],2)))
			f_out.write('\n')
	else:
		precision = 0.0001
		# first compute the capacities for a range of lambdas
		# smallest lambda needed is 1/1+alpha_max
		# and largest lambda needed is one for which capacity exceeds
		# 1/1+alpha_max
		alpha_min = np.min(alpha)
		alpha_max = np.max(alpha)
		#k_array is the range of k used in poisson summations
		k_array = np.array(range(1,1001))
		exp_log_term = np.zeros(np.shape(k_array))
		for k in k_array:
			k_0 = np.array(range(0,k+1))
			binom_pdf = scipy.stats.binom.pmf(k_0,k,eps)
			exp_log_term[k-1] = \
			np.sum(binom_pdf*np.logaddexp2(0,(2*k_0-k)*np.log2(1/eps-1)))
		
		lamb = np.floor(1/(1+alpha_max)*1/precision)*precision
		# round down to resolution
		lambda_counts = []
		capacity_counts = []
		while True:
			poisson_pdf = scipy.stats.poisson.pmf(k_array, lamb)
			capacity = 1-np.exp(-lamb)-np.sum(poisson_pdf*exp_log_term)
			lambda_counts.append(lamb)
			capacity_counts.append(capacity)
			if capacity > 1/(1+alpha_min):
				break
			lamb += precision
		lambda_counts = np.array(lambda_counts)
		capacity_counts = np.array(capacity_counts)
		coverage_counts = np.zeros(np.size(alpha))
		for i in range(np.size(alpha)):
			idx = np.nonzero(capacity_counts > 1/(1+alpha[i]))[0][0]
			coverage_counts[i] = lambda_counts[idx]*(1+alpha[i])
		
		bin_ent_eps = np.zeros(np.shape(k_array))
		for k in k_array:
			if k%2 == 0: #special handling for k/2
				error_prob_k = scipy.stats.binom.cdf(k/2-1,k,1-eps)
				error_prob_k += 0.5*scipy.stats.binom.pmf(k/2,k,1-eps)
			else:
				error_prob_k = scipy.stats.binom.cdf(k/2,k,1-eps)
			bin_ent_eps[k-1] = scipy.special.entr(error_prob_k)/np.log(2)\
					+ scipy.special.entr(1-error_prob_k)/np.log(2)
		
		lamb = np.floor(1/(1+alpha_max)*1/precision)*precision
		# round down to resolution
		lambda_maj = []
		capacity_maj = []
		while True:
			poisson_pdf = scipy.stats.poisson.pmf(k_array, lamb)
			capacity = 1-np.exp(-lamb)-np.sum(poisson_pdf*bin_ent_eps)
			lambda_maj.append(lamb)
			capacity_maj.append(capacity)
			if capacity > 1/(1+alpha_min):
				break
			lamb += precision
		lambda_maj = np.array(lambda_maj)
		capacity_maj = np.array(capacity_maj)
		coverage_maj = np.zeros(np.size(alpha))
		for i in range(np.size(alpha)):
			idx = np.nonzero(capacity_maj > 1/(1+alpha[i]))[0][0]
			coverage_maj[i] = lambda_maj[idx]*(1+alpha[i])
		
		f_out = open(outfile,'w')
		f_out.write("c_w,c_r\n")
		for i in range(np.size(alpha)):
			f_out.write(str("%.2f" % round(1+alpha[i],2)))
			f_out.write(',')
			f_out.write(str("%.2f" % round(coverage_counts[i],2)))
			f_out.write('\n')

def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    compute_capacity(config.eps, config.outfile)	

if __name__ == '__main__':
    main()
