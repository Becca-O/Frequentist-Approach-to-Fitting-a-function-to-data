#============================================
# program: fitdata.py
# purpose: calls linfit and glinfit, on a
#  set of data provided by user. Prints out
#  parameters, uncertainties, chi^2, and Q.
#============================================
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import linfit

#============================================
# reads table of measured data. Expects
# three columns with (x, y, sig)
# Returns arrays x,y,sig.
def readdata(name):
    f   = open(name)
    lst = []
    for line in f:
        lst.append([float(d) for d in line.split()])
    x   = np.array([d[0] for d in lst])
    y   = np.array([d[1] for d in lst])
    sig = np.array([d[2] for d in lst])
    return x,y,sig

#============================================
# MODEL FUNCTION
# Must return a vector of the values of the 
# basis functions at position x. The sum over the vector should
# be the value of the model function at position x.
# (Think expansion in terms of
# basis functions), i.e. ymod = np.sum((ftest(x))
def ftest(x):
    if (isinstance(x,np.ndarray)):
        n = x.size
        return np.array([[np.zeros(n)+1.0],[x],[np.sin(x)]])
    else:
        return np.array([1.0,x,np.sin(x)])

#============================================
# main() should read all 5 data sets, perform
# the fits using linfit and glinfit,  plot the 
# results, and print out the fit parameters and
# their uncertainties.
# The routines linfit and glinfit can be addressed via 
# linfit.linfit(...) and linfit.glinfit(...)
def main():

    # ???????????????????????????????????????
    
    # LINFIT
    data0 = np.loadtxt('data0.txt')
    x_0 = data0[0:len(data0), 0]
    y_0 = data0[0:len(data0), 1]
    sig_0 = data0[0:len(data0), 2]
    a_0, b_0, siga_0, sigb_0, chi2_0, q_0 = linfit.linfit(x_0, y_0, sig_0)
    print(a_0, b_0, siga_0, sigb_0, chi2_0, q_0)
    # PLOT 0
    y_fit0 = a_0 + (b_0 * x_0)
    plt.plot(x_0, y_0, 'o', x_0, y_fit0, '-')
    plt.errorbar(x_0, y_0, sig_0, fmt = '|')
    plt.legend(('Data', 'Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 0')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    data1 = np.loadtxt('data1.txt')
    x_1 = data1[0:len(data1), 0]
    y_1 = data1[0:len(data1), 1]
    sig_1 = data1[0:len(data1), 2]
    a_1, b_1, siga_1, sigb_1, chi2_1, q_1 = linfit.linfit(x_1, y_1, sig_1)
    print(a_1, b_1, siga_1, sigb_1, chi2_1, q_1)
    # PLOT 1
    y_fit1 = a_1 + (b_1 * x_1)
    plt.plot(x_1, y_1, 'o', x_1, y_fit1, '-')
    plt.errorbar(x_1, y_1, sig_1, fmt = '|')
    plt.legend(('Data', 'Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 1')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    data2 = np.loadtxt('data2.txt')
    x_2 = data2[0:len(data2), 0]
    y_2 = data2[0:len(data2), 1]
    sig_2 = data2[0:len(data2), 2]
    a_2, b_2, siga_2, sigb_2, chi2_2, q_2 = linfit.linfit(x_2, y_2, sig_2)
    print(a_2, b_2, siga_2, sigb_2, chi2_2, q_2)
    # PLOT 2
    y_fit2 = a_2 + (b_2 * x_2)
    plt.plot(x_2, y_2, 'o', x_2, y_fit2, '-')
    plt.errorbar(x_2, y_2, sig_2, fmt = '|')
    plt.legend(('Data', 'Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    data3 = np.loadtxt('data3.txt')
    x_3 = data3[0:len(data3), 0]
    y_3 = data3[0:len(data3), 1]
    sig_3 = data3[0:len(data3), 2]
    a_3, b_3, siga_3, sigb_3, chi2_3, q_3 = linfit.linfit(x_3, y_3, sig_3)
    print(a_3, b_3, siga_3, sigb_3, chi2_3, q_3)
    # PLOT 3
    y_fit3 = a_3 + (b_3 * x_3)
    plt.plot(x_3, y_3, 'o', x_3, y_fit3, '-')
    plt.errorbar(x_3, y_3, sig_3, fmt = '|')
    plt.legend(('Data', 'Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 3')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    data4 = np.loadtxt('data4.txt')
    x_4 = data4[0:len(data4), 0]
    y_4 = data4[0:len(data4), 1]
    sig_4 = data4[0:len(data4), 2]
    a_4, b_4, siga_4, sigb_4, chi2_4, q_4 = linfit.linfit(x_4, y_4, sig_4)
    print(a_4, b_4, siga_4, sigb_4, chi2_4, q_4)
    # PLOT 4
    y_fit4 = a_4 + (b_4 * x_4)
    plt.plot(x_4, y_4, 'o', x_4, y_fit4, '-')
    plt.errorbar(x_4, y_4, sig_4, fmt = '|')
    plt.legend(('Data', 'Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 4')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()

    # ---------------------------------------------------------------------
    
    # GLINFIT
    m = 3
    
    fMOD = ftest
    a_0, siga_0, chi2_0, q_0 = linfit.glinfit(x_0, y_0, sig_0, m, fMOD)
    print(a_0, siga_0, chi2_0, q_0)
    # PLOT 0
    y_fit0 = a_0[0] + (a_0[1] * x_0) + (a_0[2] * np.sin(x_0))
    plt.plot(x_0, y_0, 'o', x_0, y_fit0, '-')
    plt.errorbar(x_0, y_0, sig_0, fmt = '|')
    plt.legend(('Data', 'General Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 0')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    a_1, siga_1, chi2_1, q_1 = linfit.glinfit(x_1, y_1, sig_1, m, fMOD)
    print(a_1, siga_1, chi2_1, q_1)
    # PLOT 1
    y_fit1 = a_1[0] + (a_1[1] * x_1) + (a_1[2] * np.sin(x_1))
    plt.plot(x_1, y_1, 'o', x_1, y_fit1, '-')
    plt.errorbar(x_1, y_1, sig_1, fmt = '|')
    plt.legend(('Data', 'General Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 1')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()

    a_2, siga_2, chi2_2, q_2 = linfit.glinfit(x_2, y_2, sig_2, m, fMOD)
    print(a_2, siga_2, chi2_2, q_2)
    # PLOT 2
    y_fit0 = a_2[0] + (a_2[1] * x_2) + (a_2[2] * np.sin(x_2))
    plt.plot(x_2, y_2, 'o', x_2, y_fit2, '-')
    plt.errorbar(x_2, y_2, sig_2, fmt = '|')
    plt.legend(('Data', 'General Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    a_3, siga_3, chi2_3, q_3 = linfit.glinfit(x_3, y_3, sig_3, m, fMOD)
    print(a_3, siga_3, chi2_3, q_3)
    # PLOT 3
    y_fit3 = a_3[0] + (a_3[1] * x_3) + (a_3[2] * np.sin(x_3))
    plt.plot(x_3, y_3, 'o', x_3, y_fit3, '-')
    plt.errorbar(x_3, y_3, sig_3, fmt = '|')
    plt.legend(('Data', 'General Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 3')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    a_4, siga_4, chi2_4, q_4 = linfit.glinfit(x_4, y_4, sig_4, m, fMOD)
    print(a_4, siga_4, chi2_4, q_4)
    # PLOT 4
    y_fit4 = a_4[0] + (a_4[1] * x_4) + (a_4[2] * np.sin(x_4))
    plt.plot(x_4, y_4, 'o', x_4, y_fit4, '-')
    plt.errorbar(x_4, y_4, sig_4, fmt = '|')
    plt.legend(('Data', 'General Least-Squares Fit', 'y Error'),loc = 0)
    plt.title('Data Set 4')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    # ????????????????????????????????????
    
#========================================

main()