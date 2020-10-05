#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to estimate error using the monte carlo technique
# Created on Sat Jul  4 12:40:59 2020

import pandas_montecarlo as pm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# terms in equation
a = pd.Series(np.random.normal(size = 1000,scale = 2)) + 18
b = pd.Series(np.random.normal(size = 1000, scale = 0.1)) + 19
# errors
error_a = 0.5
error_b = 1
# calculation / equation
c = a * b
# calculate data_1 and data_2 incorporating 
# randomly select a value for a and b wihtin error range

# a simulated
a_sim = a.montecarlo(sims=1000, bust=-0.5, goal = 1)


mc = c.montecarlo(sims=1000, bust=-0.1, goal=1)



