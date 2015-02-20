#! /usr/bin/python
"""
model_comparison.py

Functions to calculate AIC, AICc, and BIC

lnl = log-likelihood
p = number of parameters in model
n = sample size

Copyright 2015 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

def aic(lnl, p):
    return round( ( -2 * lnl ) - ( 2 * p ), 2 )

def aicc(lnl, p, n):
    return round( ( -2 * lnl ) + ( ( 2 * p ) * ( n / (n - p - 1) ) ), 2 )

def bic(lnl, p, n):
    import math
    return round( ( -2 * lnl ) + ( p * math.log(n) ), 2 )

def print_ic_values(lnl, p, n):
    print("AIC:\t" + str( aic(lnl, p) ) )
    print("AICc:\t" + str( aicc(lnl, p, n) ) )
    print("BIC:\t" + str( bic(lnl, p, n) ) )
