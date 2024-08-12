# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 10:49:45 2024

@author: charlesr
"""

import numpy as np
from scipy.optimize import fsolve


def get_bond_price_from_cdf(coupon, tenor, recovery_rate, riskfree_rate, cdf):
    
    # Get times of coupon payments
    coupon_times = np.arange(tenor, 0, -0.5)
    
    # Partiton interval [0, T] uniformly 
    time_partition = np.linspace(0, tenor, 2 * len(coupon_times) + 1)
    
    # Initialize bond price
    bond_price = 0
    
    # Loop over subintervals
    for t0, t1 in zip(time_partition[:-1], time_partition[1:]):
        
        # Cash flow if default times probability of default in [t0, t1]
        bond_price += recovery_rate * np.exp(-t1 * riskfree_rate) * (cdf(t1) - cdf(t0))
        
        # Get payments in interval
        payment_times = coupon_times[(t0 < coupon_times) & (coupon_times <= t1)]
        
        # PV of payments times probabilities of making payments. Sum results
        bond_price += coupon/2 * np.sum([(1 - cdf(t)) * np.exp(-riskfree_rate * t) 
                                         for t in payment_times])
        
    # Add principal payment times probability of payment
    bond_price += 100 * np.exp(-riskfree_rate * tenor) * (1 - cdf(tenor))
    
    return bond_price


def get_bond_price_from_spread(coupon, tenor, recovery_rate, riskfree_rate, spread):
    
    # Get times of coupon payments
    coupon_times = np.arange(tenor, 0, -0.5)

    # Get cash flows
    cash_flows = np.repeat(coupon/2, len(coupon_times)) 
    
    # Add principal payment to final cash flow
    cash_flows[0] += 100
    
    # Calculate discount factors
    discount_factors = np.exp(-(riskfree_rate + spread) * coupon_times)
    
    # Get PV
    bond_price = np.sum(cash_flows * discount_factors)
    
    return bond_price


def get_spread(coupon, tenor, recovery_rate, riskfree_rate, cdf, spread0 = 0.05):
    
    # Get bond price using cdf
    bond_price = get_bond_price_from_cdf(coupon, tenor, recovery_rate, 
                                             riskfree_rate, cdf)
    
    try:
        
        # Get objective function
        obj = lambda spread: bond_price - get_bond_price_from_spread(coupon, 
                                tenor, recovery_rate, riskfree_rate, spread)
        
        # Solve for spread
        spread = 100**2 * fsolve(obj, x0 = spread0)[0]
        
        return spread
    
    except:
        
        return np.nan
    
        
# =============================================================================
# tenor = 10
# coupon = 5
# recovery_rate = 40
# riskfree_rate = 0.04
# 
# cdf = lambda t: 1 - np.exp(-0.10 * t)
# =============================================================================
        
    
    

