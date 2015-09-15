##############################################################################
### To use the functions in this lib simply import this python module using
### import myStats
### Then you'll able able to call functions with the proper arguments using
### returnVal=myStats.func1(arg1,arg2)
##############################################################################
##############################################################################
import numpy as np

################### FUNC benjamini_hochberg_correction  #####################
#### Given an array of p-values (not necessarily sorted) and the number of total 
### tests that were performed to gather these p-values, this function performs
###  the multiple hypothesis testing correction described by Benjamini-Hochberg.
###
### 
### If the number of tests are much more compared to p-value array and
### the omitted p-values are all supposed to be zero you should use it like:
### q_array=benjamini_hochberg_correction([0.03,0.4,0.7,0.01],10)
### 
### If the number of tests are equal to the ones in p-values array then:
### p_array=[0.03,0.4,0.7,0.01]
### q_array=benjamini_hochberg_correction(p_array,len(p_array))
##############################################################################
def benjamini_hochberg_correction(p_values, num_total_tests):
	# assumes that p_values vector might not be sorted
	pvalsArray=np.array(p_values)
	order=pvalsArray.argsort()
	sorted_pvals=np.take(p_values,order)
	q_values=[1.0 for i in range(len(p_values))]
	prev_bh_value = 0
	for i, p_value in enumerate(sorted_pvals):
		bh_value = p_value * num_total_tests / (i + 1)
		# Sometimes this correction can give values greater than 1,
		# so we set those values at 1
		bh_value = min(bh_value, 1)

		# To preserve monotonicity in the values, we take the
		# maximum of the previous value or this one, so that we
		# don't yield a value less than the previous.
		bh_value = max(bh_value, prev_bh_value)
		prev_bh_value = bh_value
		qIndex=order[i]
		q_values[qIndex] =bh_value
	#END for
	return q_values

################### FUNC meanAndVariance  ##################################
#### Given an array this function calculates the mean and variance 
##############################################################################
def meanAndVariance(a):
	#calculates variance using E(x^2) - (Ex)^2
	sumxsquared=0
	sumx=0
	for x in a:
		sumxsquared += (x*x)
		sumx += x
	mean=sumx/float(len(a))
	squaremean=sumxsquared/float(len(a))
	variance=squaremean - (mean*mean)
	return (mean, variance)


