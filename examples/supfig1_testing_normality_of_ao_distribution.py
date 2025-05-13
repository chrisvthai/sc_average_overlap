import matplotlib.pyplot as plt
import numpy as np
import random
import math
from typing import List, Optional, Union

def compute_overlap(l1: Union[List, np.ndarray], l2: Union[List, np.ndarray]):
	assert len(l1) == len(l2), "lists must be equal in length"
	return len(set(l1).intersection(l2)) / len(l1)

def average_overlap(l1: Union[List, np.ndarray], l2: Union[List, np.ndarray]):
	assert type(l1) in [list, np.ndarray], "first array must be np.ndarray or List"
	assert type(l2) in [list, np.ndarray], "second array must be np.ndarray or List"
	assert len(l1) == len(l2), "lists must be equal in length"

	total_depth = len(l1)

	# First generate the ranked lists, which are different from the vectors themselves
	# argsort gives indices of the largest elements. this is in increasing order, so list reversal is needed
	#l1_ranks = np.argsort(l1)[::-1]
	#l2_ranks = np.argsort(l2)[::-1]
	l1_ranks = l1
	l2_ranks = l2 

	ao = 0
	overlaps_array = []

	for d in range(1, total_depth+1):
		o = compute_overlap(l1_ranks[:d], l2_ranks[:d])
		overlaps_array.append(o)
		ao += o

	ao /= total_depth
	return ao, overlaps_array


list_lengths = [100,  500,] # 1000, 2000]
base_list = list(range(0, 2000))

# 1. generate distributions for each list length to be tested
n_samples = 1000
ao_score_dist_list = {}
for l in list_lengths:
	l1 = base_list[:l]
	l2 = base_list[:l]

	ao_score_dist = []
	for i in range(0, n_samples):
		random.shuffle(l1)
		random.shuffle(l2)

		ao, _ = average_overlap(l1, l2)
		ao_score_dist.append(ao)

	ao_score_dist_list[l] = ao_score_dist

# Fit normal distributions and calculate standard deviations
for l in list_lengths:

	plt.hist(ao_score_dist_list[l], bins=20)
	plt.show()
