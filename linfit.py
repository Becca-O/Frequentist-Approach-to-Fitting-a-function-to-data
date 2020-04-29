#============================================
# program: Library for fitting functions with
# linear coefficients (linfit). If the 
# functions are non-linear, use glinfit.
#============================================
import numpy as np
import scipy.special
#===================================
# function: linfit
# purpose : fits a straight line with
#           parameters a+b*x to data set.
# input   : x   : float vector of length n: "independent" variable (assumed to have no uncertainties)
#           y   : float vector of length n: data points
#           sig : float vector of length n: measurement uncertainties for data points y.
# output  : a    : float number: fit parameter (here: offset)
#           b    : float number: fit parameter (here: slope)
#           siga : float number: uncertainty of a
#           sigb : float number: uncertainty of b
#           chi2 : floating point number. Should be around n-2, i.e. 
#                  the number of data points less the number of parameters ("degrees of freedom")
#           q    : quality of fit estimate (should be between ~0.1 and 10^(-3))
#===================================
def linfit(x,y,sig):
    
    # ??????????????????????????????
    
    S = 0.0
    S_xx = 0.0
    S_y = 0.0
    S_xy = 0.0
    S_x = 0.0
    
    n = x.size # number of data points
    m = 2 # degrees of freedom (number of parameters)
    for i in range(0, n-1):
        S = S + (1/(sig[i]**2))
        S_x = S_x + (x[i]/(sig[i]**2))
        S_y = S_y + (y[i]/(sig[i]**2))
        S_xx = S_xx + ((x[i]**2)/(sig[i]**2))
        S_xy = S_xy + ((x[i]*y[i])/(sig[i]**2))
        
    omega = (S * S_xx) - (S_x**2)
    a = ((S_xx * S_y) - (S_xy * S_x)) / omega # the offset
    b = ((S * S_xy) - (S_x * S_y)) / omega # the slope
    
    siga = np.sqrt(S_xx/omega) # the error in the offset
    sigb = np.sqrt(S/omega) # the error in the slope
    
    # χ2
    chi2 = 0.0
    for i in range(0, n-1):
        expected = a + b*x[i]
        chi2 = chi2 + ((y[i] - expected)/(sig[i]))**2
         
    q = scipy.special.gammaincc(0.5*float(n-m),0.5*chi2)
    # PROBABILITY: how likely it would be to find another data set with a χ2
    # as poor as the one just found, assuming that the model is correct

    return a, b, siga, sigb, chi2, q

    # ??????????????????????????????


#===================================
# function: glinfit
# purpose : fits a set of basis functions with linear 
#           coefficients to a data set.
#           Measurement uncertainties need to be provided.
# input   : x   : float vector of length n: "independent" variable (assumed to have no uncertainties)
#           y   : float vector of length n: data points
#           sig : float vector of length n: measurement uncertainties for data points y.
#           m   : integer: number of parameters
#           fMOD: function pointer vector of length m. (f(x))[j] must return
#                the j'th basis function, j=0...m-1. For example, if
#                we have a function basis [1,x,x^2] (parabolic fit),
#                then (f(x0))[0] = 1, (f(x0))[1] = x0, (f(x0))[2] = x0^2
# output  : a    : float vector of length m: fit parameters (here: offset and slope)
#           siga : float vector of length n: uncertainties of a ("errors")
#           chi2 : floating point number. Should be around n-m, i.e. 
#                  the number of data points less the number of parameters ("degrees of freedom")
#           q    : quality of fit estimate (should be between ~0.1 and 10^(-3))
#============================================
def glinfit(x,y,sig,m,fMOD):
    # ??????????????????????????????
    
    n = len(x)
    A = np.zeros([n,m])
    b_vector = np.zeros([n,1])
    siga = np.zeros([n])
    for i in range(0, n-1):
        A[i,0] = 1/sig[i]
        A[i,1] = x[i]/sig[i]
        A[i,2] = np.sin(x[i])/sig[i]
        b_vector[i] = y[i]/sig[i]
    A_T = np.matrix.transpose(A)
    C = np.linalg.inv(np.matmul(A_T, A))
    a_vector = np.matmul(C, np.matmul(A_T, b_vector))
    a = a_vector[0,0]
    b = a_vector[1,0]
    c = a_vector[2,0]
    
    siga = np.zeros(m)
    for i in range(m):
        siga[i] = np.sqrt(C[i,i])

    # χ2
    chi2 = 0.0
    for i in range(n):
        expected = a + (b*x[i]) + (c*np.sin(x[i]))
        chi2 = chi2 + ((y[i] - expected)/(sig[i]))**2

    q = scipy.special.gammaincc(0.5*float(n-m),0.5*chi2)
    # PROBABILITY: how likely it would be to find another data set with a χ2
    # as poor as the one just found, assuming that the model is correct
    
    return a_vector, siga, chi2, q
    
    # ??????????????????????????????