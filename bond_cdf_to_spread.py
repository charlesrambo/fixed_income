# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 10:49:45 2024

@author: charlesr
"""

import numpy as np
from scipy.optimize import fsolve


def get_bond_price_from_cdf(coupon, tenor, recovery_rate, riskfree_rate, cdf):
    
    # Get times of coupon payments
    t_coupons = np.arange(tenor, 0, -0.5)
    
    # Partiton interval [0, tenor] uniformly 
    t_partition = np.linspace(0, tenor, 2 * len(t_coupons) + 1)
    
    # Add in coupons in cases where there's no default before time t
    bond = coupon/2 * np.sum([(1 - cdf(t)) * np.exp(-riskfree_rate * t) 
                                          for t in t_coupons])    
        
    # Add principal payment in cases where there's no default during the life of bond
    bond += 100 * np.exp(-riskfree_rate * tenor) * (1 - cdf(tenor))
    
    # Add in value in cases where there's a default
    bond += recovery_rate * np.sum([np.exp(-riskfree_rate * t1) * (cdf(t1) - cdf(t0)) 
                        for t0, t1 in zip(t_partition[:-1], t_partition[1:])])
    
    return bond


def get_bond_price_from_spread(coupon, tenor, recovery_rate, riskfree_rate, spread):
    
    # Get times of coupon payments
    t_coupons = np.arange(tenor, 0, -0.5)

    # Get cash flows
    cash_flows = np.repeat(coupon/2, len(t_coupons)) 
    
    # Add principal payment to final cash flow
    cash_flows[0] += 100
    
    # Calculate discount factors
    discount_factors = np.exp(-(riskfree_rate + spread) * t_coupons)
    
    # Get PV
    bond = np.sum(cash_flows * discount_factors)
    
    return bond


def get_spread(coupon, tenor, recovery_rate, riskfree_rate, cdf, spread0 = 0.05):
    
    # Get bond price using cdf
    bond = get_bond_price_from_cdf(coupon, tenor, recovery_rate, 
                                             riskfree_rate, cdf)
    
    try:
        
        # Get objective function
        obj = lambda spread: bond - get_bond_price_from_spread(coupon, tenor, 
                                        recovery_rate, riskfree_rate, spread)
        
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
# cdf = lambda t: 1 - np.exp(-0.03 * t)
# =============================================================================
        
    
    

