# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 13:54:13 2025

@author: USER
"""

import numpy as np

def KPM_Soln(a,b,V0,eps_range):
    
    h_bar = 1.054*1e-34
    m = 9.109*1e-31
    alpha_0 = (2*m*V0*1.602*1e-19/h_bar**2)**(1/2)
    
    def KPM_p(eps):
        reutrn (1-2*eps)/(2*(eps*(eps-1))**(1/2))*np.sin(alpha_0*a*1e-10*eps)