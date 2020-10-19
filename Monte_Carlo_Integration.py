from scipy import random
import numpy as np

"""MONTE-CARLO INTEGRATION (SAMPLING METHOD)-----------------------------------"""

# Requires lower and upper limits of integration and the number of randomly
# generated sample values between these limits, and the FUNCTION you wish to 
# integrate.

def MC_int(func,lwr,upp,N):
    
    # Generate list of random numbers with values between integration limits
    rng_list= random.uniform(lwr, upp, N)
    
    # start value of our integration
    integral = 0.0 
    sumsq = 0.0                     

    # Perform evaluation of function for all random numbers listed in rng_list,
    # and do the same for the square of the function
    for i in range(N):                 
        x = (rng_list[i])               
        integral += func(x)
        sumsq += func(x)**2
    
    # Calculation below evaluates the integration and associated error
    answer = ((upp-lwr)/N)*integral
    error = (upp-lwr)*np.sqrt(((sumsq/N)-(integral/N)**2)/N)
      
    return answer, error
