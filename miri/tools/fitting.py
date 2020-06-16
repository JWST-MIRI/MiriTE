#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""

Module fitting - Contains general purpose fitting functions for
MIRI data, including:

* gaussPlusPoly - Generate a Gaussian plus a polynomial
* gaussian - Generate a Gaussian
* nonLinFit - Nonlinear least squares fit using the Levenberg-Marquardt method
      in the scipy.optimize.curve_fit method
      this should be possible with any function in the format of
      ydata = f(xdata, *params)    
* LinFit - Class for linear fitting

:History:

09 Jan 2013: Created
02 Sep 2013: Suppress unwanted output.
09 Dec 2014: Suppress even more unwanted output.
10 Dec 2014: Replaced y_sigma == None with y_sigma is None.
08 Nov 2017: Moved from LRS pipeline package to general purpose fitting
             package.

@author:  Juergen Schreiber

"""

# Import standard python packages/modules
# e.g. import os

# Import packages/modules from the MIRI development environment
import numpy as np
from scipy import optimize, special
from matplotlib import pyplot


def gaussPlusPoly(x, *p):
    """
    Function of a gaussian peak + a Polynomial 
    
    :Parameters:
    
    x: array
       position where the function is evaluated
    
    p: array 
       * Constant Background          : p[0]
       * Peak height above background : p[1]
       * Central value                : p[2]
       * Standard deviation           : p[3]
       * plus a polynomial with coefficients p[4] + ....: starting with the highest order but without constant offset parameter
    
    """
    gp = (p[0]+p[len(p)-1], p[1], p[2], p[3])
    ga = gaussian(x, *gp)
    pList = []
    for i in range(4, len(p)):
        pList.append(p[i])
    #zero offset constant
    pList.append(0)
    pol = np.polyval(pList, x)
    return pol + ga
    

def gaussian(x, *p):
    """
    function of a gaussian peak
    
    :Parameters:
    
    x: array
       position where the gaussian is evaluated
    p: array
       * Constant Background          : p[0]
       * Peak height above background : p[1]
       * Central value                : p[2]
       * Standard deviation           : p[3]
    
    """
    return p[0]+p[1]*np.exp(-1*(x-p[2])**2/(2*p[3]**2))


def nonLinFit(func, x_data, y_data, y_sigma=None, p_guess=None, plotting=False,
              verbose=True):     
    """
    non linear fit function
    
    :Parameters:
    
    func: callable
        The model function, func(x, ...).  It must take the independent
        variable as the first argument and the parameters to fit as
        separate remaining arguments.
    
    x_data: array 
        An N-length sequence or an (k,N)-shaped array
        for functions with k predictors. The independent variable where the data is measured.
    
    y_data: array 
        N-length sequence, the dependent data: nominally func(xdata, ...)
    
    p_guess: array, optional (default=None) 
        None, scalar, or M-length sequence
        Initial guess for the parameters.  If None, then the initial
        values will all be 1 (if the number of parameters for the 
        function can be determined using introspection, otherwise a 
        ValueError is raised).
    
    y_sigma: array, optional (default = None)
        None or N-length sequence
        If not None, it represents the standard-deviation of y_data.
        This vector, if given, will be used as weights in the
        least-squares problem. If None, they are set to unity.
    
    verbose: boolean, optional (default = True)
        in verbose mode it prints out
            * the optimal values for the parameters so that the sum of the 
              squared error of ``func(xdata, *popt) - ydata`` is minimized
              and their uncertainties
            * the estimated covariance and correlation matrices
              the diagonals provide the variance of the parameter estimate.
            * Chi-Squared
            * DOF: degree of freedom
            * CDF: cumulative distribution function
    
    plotting: boolean, optional (default = False)
        plotting results or not, it overplots:
            * the guess function, 
            * the measured values (inclusive error bars)
            * the final fit 
            * the residual with error bars   
    
    :Returns:
    
    fit: array 
        the computed fit function
    
    popt: array
        the optimal values for the parameters so that the sum of the 
        squared error of ``func(xdata, *popt) - ydata`` is minimized
        and their uncertainties
    
    :Raises:
    
    TypeError
        if the fit is bad
    """
    
    x_func = np.linspace(min(x_data),max(x_data))
    
    initial_plot = func(x_func,*p_guess)
    
    try:
        # Notes: maxfev is the maximum number of func evaluations tried; you
        # can try increasing this value if the fit fails.
        p, cov = optimize.curve_fit(func, x_data, y_data, p0 = p_guess, sigma = y_sigma, maxfev=100*(len(x_data)+1))
    except:
        print("scipy.optimize.curve_fit failed")
        p, cov = p_guess, None       
        
    y_fit = func(x_data,*p)
    y_residual = y_data - y_fit

    perr = np.zeros(len(p))
    # Calculate degrees of freedom of fit
    dof = len(x_data) - len(p) 
 
    ## Output results
     
    if verbose:
        print("\nNumber of Data Points = %7g, Number of Parameters = %1g"\
            %(len(x_data), len(p) ))

        print("Covariance Matrix : \n", cov, "\n")
    if y_sigma is None: 
        y_sigma = 1
    try:
        if verbose:
            print("Correlation Matrix :")
        for i,row in enumerate(cov):
            for j in range(len(p)) :
                if verbose:
                    print("%10f" % (cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])), end=' ')
            if verbose:
                print() 
        # Calculate Chi-squared
        chisq = sum(((y_data-func(x_data,*p))/y_sigma)**2)
        #calculate parameter errors
        for i in range(len(p)):
            perr[i] = cov[i,i]**0.5*max(1,np.sqrt(chisq/dof))
        
        if verbose:
            print("\nEstimated parameters and uncertainties (with initial guesses)")
            for i in range(len(p)) :
                print(("   p[%d] = %10.5f +/- %10.5f      (%10.5f)"\
                       %(i,p[i],cov[i,i]**0.5*max(1,np.sqrt(chisq/dof)), p_guess[i])))

            print("Chi-Squared/dof = %10.5f, CDF = %10.5f%%"\
                %(chisq/dof, 100.*float(special.chdtrc(dof,chisq))))
            if chisq > dof :
                print("\nNOTE:Because Chi-squared > dof, the parameter uncertainty")
                print("      estimates have been scaled up by sqrt(Chi-squared/dof).")
    
    # If cov has not been calculated because of a bad fit, the above block
    #   will cause a python TypeError which is caught by this try-except structure.
    except TypeError:
        print("**** BAD FIT ****")
        print("Parameters were: ",p)
        # Calculate Chi-squared for current parameters
        chisq = sum(((y_data-func(x_data,*p))/y_sigma)**2)
        print("Chi-Squared/dof for these parameter values = %10.5f, CDF = %10.5f%%"\
            %(chisq/dof, 100.*float(special.chdtrc(dof,chisq))))
        print("Uncertainties not calculated.")
        print()
        print("Try a different initial guess for the fit parameters.")
        print("Or if these parameters appear close to a good fit, try giving")
        print("    the fitting program more time by increasing the value of maxfev.")
        chisq = None

    if plotting:
        ## Plot
        # create figure with light gray background
        fig = pyplot.figure(facecolor="0.98")
        # 3 rows, 1 column, subplot 1
        # 3 rows are declared, but there are only 2 plots; this leaves room for text
        #   in the empty 3rd row
        fit = fig.add_subplot(311)
        # remove tick labels from upper plot (for clean look)
        fit.set_xticklabels( () ) 

        # Plot data as red circles, and fitted function as (default) line
        # (The sort is in case the x data are not in sequential order.)
        fit.plot(x_data,y_data,'ro', np.sort(x_data), func(np.sort(x_data),*p))
        # draw starting guess as dashed green line ('r-')
        fit.plot(x_func, initial_plot, 'g-', label="Start", linestyle="--") 
        # Add error bars on data as red crosses.
        fit.errorbar(x_data, y_data, yerr=y_sigma, fmt='r+')

        # separate plot to show residuals
        residuals = fig.add_subplot(312) # 3 rows, 1 column, subplot 2
        residuals.errorbar(x_data, y_residual, yerr=y_sigma, fmt='r+', label="Residuals")
        # make sure residual plot has same x axis as fit plot
        residuals.set_xlim(fit.get_xlim())
        residuals.axhline(y=0) # draw horizontal line at 0 on vertical axis

        # These data look better if 'plain', not scientific, notation is used,
        #   and if the tick labels are not offset by a constant (as is done by default).
        #   Note: This only works for matplotlib version 1.0 and newer, so it is
        #           enclosed in a "try" to avoid errors.
        try:
            pyplot.ticklabel_format(style='plain', useOffset=False, axis='x')
        except: pass

        # print selected information in empty 3rd plot row
        try:
            pyplot.figtext(0.05,0.25,"Converged with ChiSq = " + str(chisq) + ", DOF = " +
                str(dof) + ", CDF = " + str(100*special.chdtrc(dof,chisq))+"%")
            for i, value in enumerate(p):
                pyplot.figtext(0.08,0.16-i*0.03, "p["+str(i)+"]" + " = " +
                           str(p[i]).ljust(18) + " +/- " +
                           str(np.sqrt(cov[i,i])*max(1,np.sqrt(chisq/dof))),
                           fontdict=None)
                # Note: Including family="Monospace" in the above figtext call will
                #       produce nicer looking output, but can cause problems with
                #       some older python installations.
        except TypeError:
            pyplot.figtext(0.05,0.25,"BAD FIT.  Guess again.")


        # Display the plot
        pyplot.show()
        # To print the plot, save it first, and then print the saved image.
        # Closing the plot window will end the python script.
    return y_fit, p, perr


class LinFit():
    """
    Linear Fit class does a
    linear least square fit of x, y arrays
    ``y = a*x + b``
    
    :Synopsis:
    
      * fitter = LinFit(x, y)
      * method fit: fitter.fit()
      * method fitter.calcLinFitError():
        calculates the error on the whole fitted straight line
    
    :Returns: 
    
    returns the coefficients array
       * coeff[0] = a
       * coeff[1] = b
       * coeff[2] = y_err
       * coeff[3] = a_err
       * coeff[4] = b_err
    """

    def __init__(self, x_data, y_data):
        """
        constructor
        
        :Parameters:
        
        x_data: x coordinates array
        y_data: corresponding y coordinates array
        """
        self.x = x_data
        self.y = y_data
    
    def fit(self):
        """
        linear least square fit of x, y arrays
        ``y = a*x + b``
        
        :Returns:
        
        coeff: array
           returns the coefficients (0(a),1(b), the overall error on y (2), and the
           errors on the coefficients (3(a),4(b))
        
        """
        x = self.x
        y = self.y
        n = len(x)
        if n != len(y):
            raise Exception("number of elements of x and y arrays are different")
    
        #method described in 'Walcher: Physikalisches Praktikum, Kapitel 1.2.3 "statistische Fehler"
        sumxy = np.sum(x * y) 
        sumy = np.sum(y)
        sumx = np.sum(x)
        sumsqrx = np.sum(x*x)
    
        a = (n*sumxy - sumx*sumy)/(n*sumsqrx - sumx*sumx)
        b = (sumsqrx*sumy - sumx*sumxy)/(n*sumsqrx - sumx*sumx)
        erry = np.sqrt(np.sum(np.square(y - a*x - b))/(n - 2.))
        erra = erry * np.sqrt(n/(n*sumsqrx - sumx*sumx))
        errb = erry * np.sqrt(sumsqrx/(n*sumsqrx - sumx*sumx))
        self.coeff = np.array([a,b, erry, erra, errb])
        return self.coeff

    def calcLinFitError(self, x_data):
        """
        calculates the error depending on input x data and errors on linear coefficients
        
        :Returns:
        
        err: array
            errors on y data with same length as input x data
        """
        coeffs = self.coeff
        y_max = (coeffs[0]+coeffs[3])*x_data + coeffs[1]-coeffs[4]
        y_min = (coeffs[0]-coeffs[3])*x_data + coeffs[1]+coeffs[4] 
        return np.abs(y_max - y_min)

    def straightLine(self, x_data):
        """
        calculates the straight line ``y = a*x + b`` at array x=x_data
        """
        return self.coeff[0] * x_data + self.coeff[1]
        
           

# A minimal test is run when this file is run as a main program.
if __name__ == '__main__':
    print("Testing the module\n")
    PLOTTING = True     # Set to False to turn off plotting.
    print("test nonlinear fit")
    print("gaussian fit")
    x_data = np.arange(0,100)
    p = (10, 100, 50, 10)
    y_data = gaussian(x_data, *p)
    p_guess = (2, 102, 48, 7)
    fit, presult, perr = nonLinFit(gaussian, x_data, y_data, p_guess = p_guess, plotting = True)
    assert fit is not None
    
    print("gaussian + polynomial fit")
    po = (5, 100, 50, 10, 0.1, 2.)
    newy_data = gaussPlusPoly(x_data, *po)
    newp_guess = (3, 110, 45, 8, 0.05, 1.5)
    new_fit, new_presult, new_perr = nonLinFit(gaussPlusPoly, x_data, newy_data, p_guess = newp_guess, plotting = True)
    assert new_fit is not None
    
    print("linear fit")
    y = np.arange(0,100)
    fitter = LinFit(x_data, y)
    coeff = fitter.fit()
    assert coeff[0] == 1
    
    print("calcLinFitError")
    error = fitter.calcLinFitError(x_data)
    assert np.mean(error) == 0.    
    
    print("straightLine")
    line = fitter.straightLine(x_data)
    assert (np.all(line - y) == 0.)
    
    print("Test finished.")
