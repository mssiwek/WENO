#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 18:21:54 2019

@author: magda
"""

import numpy as np
from numpy.random import rand

sigma            = 0.2
def gaussian(x, sigma=sigma):
    return np.exp(-x**2/sigma)

def simple_step(x, cd = 0., left = 1., right = 0.5):
    if x <= cd:
        u = left
    if x > cd:
        u = right
    return u

def rand_step(x, cd = 0., left = 2., right = 0.):
    if x <= cd:
        u = rand(1)+left
    if x > cd:
        u = rand(1)+right
    return u

def square_well(x, cd1 = -0.5, cd2 = 0.5):
    if cd1 < x < cd2:
        u = -1.
    else:
        u = 0.
    return u

def quadratic(x, b, a):
    return a*x**2 + b

def quadratic_well(x, cd1 = -0.5, cd2 = 0.5, depth=-10, ref_ht = 0.):
    if cd1 < x < cd2:
        u = quadratic(x, depth, a = (ref_ht - depth)/(cd1**2))
    else:
        u = ref_ht
    return u

def trig_disc(x, cd = 0., f = 100.):
    if x <= cd:
        u = 3.*np.sin(x)
    if x > cd:
        u = 1.
    return u
