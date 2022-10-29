import numpy as np
import pandas as pd
import networkx as nx
import scipy as sc
from scipy import io
from matplotlib import pyplot as plt
from signet.cluster import Cluster
from signet.block_models import SSBM
from sklearn.metrics import adjusted_rand_score




mat = io.mmread('soc-sign-epinions.mtx')
data = sc.sparse.csc_matrix(mat)



A_p = np.zeros(shape = data.shape)
A_n = np.zeros(shape = data.shape)
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        if data[i, j] == -1:
            A_n[i][j] = 1;
        if data[i, j] == 1:
            A_p[i][j] = 1;