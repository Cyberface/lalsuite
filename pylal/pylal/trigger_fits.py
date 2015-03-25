from __future__ import division

# Copyright T. Dent 2015 (thomas.dent@aei.mpg.de)
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.

import numpy
from scipy.stats import kstest

def fit_exponential(vals, thresh):
    '''
    Maximum likelihood fit for the coefficient alpha for a distribution of
    discrete values p(x) = alpha exp(-alpha (x-x_t)) above a threshold x_t.

    vals: sequence of real numbers none of which lies below thresh
    thresh: threshold used in the fitting formula 
    '''
    vals = numpy.array(vals)
    if min(vals) < thresh:
        raise RuntimeError("Values to be fit must all be above threshold!")
    alpha = 1. / (numpy.mean(vals) - thresh)
    # sigma is measurement standard deviation = (-d^2 log L/d alpha^2)^(-1/2)
    sigma_alpha = alpha / len(vals)**0.5
    return alpha, sigma_alpha

def fit_rayleigh(vals, thresh):
    '''
    Maximum likelihood fit for the coefficient alpha for a distribution of
    discrete values p(x) = alpha x exp(-alpha (x**2-x_t**2)/2) above a 
    threshold x_t.

    vals: sequence of real numbers none of which lies below thresh
    thresh: threshold used in the fitting formula
    '''
    vals = numpy.array(vals)
    if min(vals) < thresh:
        raise RuntimeError("Values to be fit must all be above threshold!")
    alpha = 2. / (numpy.mean(vals**2.) - thresh**2.)
    # sigma is measurement standard deviation = (-d^2 log L/d alpha^2)^(-1/2)
    sigma_alpha = alpha / len(vals)**0.5
    return alpha, sigma_alpha

def fit_power(vals, thresh):
    '''
    Maximum likelihood fit for the coefficient alpha for a distribution of
    discrete values p(x) = ((alpha-1)/x_t) (x/x_t)**-alpha above a threshold 
    x_t.

    vals: sequence of real numbers none of which lies below thresh
    thresh: threshold used in the fitting formula
    '''
    vals = numpy.array(vals)
    if min(vals) < thresh:
        raise RuntimeError("Values to be fit must all be above threshold!")
    alpha = numpy.mean(numpy.log(vals/thresh))**-1. + 1.
    # sigma is measurement standard deviation = (-d^2 log L/d alpha^2)^(-1/2)
    sigma_alpha = (alpha - 1.) / len(vals)**0.5
    return alpha, sigma_alpha

def expfit(xvals, alpha, thresh):
    '''
    The fitted exponential function normalized to 1 above threshold

    xvals: the values at which the fit PDF is to be evaluated
    alpha: the fitted exponent factor
    thresh: threshold value for the fitting process 

    To normalize to a given total count multiply by the count.
    '''
    xvals = numpy.array(xvals)
    fit = alpha * numpy.exp(-alpha * (xvals - thresh))
    # set fitted values below threshold to 0
    numpy.putmask(fit, xvals < thresh, 0.)
    return fit

def expfit_cum(xvals, alpha, thresh):
    '''
    The integral of the exponential fit above a given value (reverse CDF)
    normalized to 1 above threshold

    To normalize to a given total count multiply by the count.
    '''
    xvals = numpy.array(xvals)
    cum_fit = numpy.exp(-alpha * (xvals - thresh))
    # set fitted values below threshold to 0
    numpy.putmask(cum_fit, xvals < thresh, 0.)
    return cum_fit

def rayleighfit(xvals, alpha, thresh):
    '''
    The fitted Rayleigh function normalized to 1 above threshold

    xvals: the values at which the fit PDF is to be evaluated
    alpha: the fitted exponent factor
    thresh: threshold value for the fitting process 

    To normalize to a given total count multiply by the count.
    '''
    if thresh < 0:
        raise RuntimeError('Threshold for Rayleigh distribution cannot be negative!')
    xvals = numpy.array(xvals)
    fit = alpha * xvals * numpy.exp(-alpha * (xvals**2. - thresh**2.) / 2.)
    # set fitted values below threshold to 0
    numpy.putmask(fit, xvals < thresh, 0.)
    return fit

def rayleighfit_cum(xvals, alpha, thresh):
    '''
    The integral of the Rayleigh fit above the x-values given (reverse CDF)
    normalized to 1 above threshold

    To normalize to a given total count multiply by the count.
    '''
    if thresh < 0:
        raise RuntimeError('Threshold for Rayleigh distribution cannot be negative!')
    xvals = numpy.array(xvals)
    cum_fit = numpy.exp(-alpha * (xvals**2. - thresh**2.) / 2.)
    # set fitted values below threshold to 0
    numpy.putmask(cum_fit, xvals < thresh, 0.)
    return cum_fit

def powerfit(xvals, alpha, thresh):
    '''
    The fitted power-law function normalized to 1 above threshold

    xvals: the values at which the fit PDF is to be evaluated
    alpha: the fitted power
    thresh: threshold value for the fitting process 

    To normalize to a given total count multiply by the count.
    '''
    xvals = numpy.array(xvals)
    fit = (alpha - 1.) * xvals**(-alpha) * thresh**(alpha - 1.)
    # set fitted values below threshold to 0
    numpy.putmask(fit, xvals < thresh, 0.)
    return fit

def powerfit_cum(xvals, alpha, thresh):
    '''
    The integral of the power-law fit above the x-values given (reverse CDF)
    normalized to 1 above threshold

    To normalize to a given total count multiply by the count.
    '''
    xvals = numpy.array(xvals)
    cum_fit = xvals**(1. - alpha) * thresh**(alpha - 1.)
    # set fitted values below threshold to 0
    numpy.putmask(cum_fit, xvals < thresh, 0.)
    return cum_fit

fitdict = {
    'exponential' : fit_exponential,
    'rayleigh'    : fit_rayleigh,
    'power'       : fit_power
}

fndict = {
    'exponential' : expfit,
    'rayleigh'    : rayleighfit,
    'power'       : powerfit
}

cum_fndict = {
    'exponential' : expfit_cum,
    'rayleigh'    : rayleighfit_cum,
    'power'       : powerfit_cum,
}

def fit_above_thresh(distr, vals, thresh=None):
    '''
    Maximum likelihood fit for the coefficient alpha for a distribution of
    discrete values p(x) = alpha exp(-alpha*x) above a given threshold.
    Values below threshold will be discarded.  If no threshold is specified, 
    the minimum sample value will be used.

    distr: name of distribution, either 'exponential' or 'rayleigh' or 'power'

    vals: sequence of real numbers

    thresh: threshold to apply before fitting - if thresh=None, use the lowest
    value
    '''
    vals = numpy.array(vals)
    if thresh == None:
        thresh = min(vals)
    else:
        vals = vals[vals >= thresh]
    return fitdict[distr](vals, thresh)

def fit_tail(distr, vals, N=1000):
    '''
    Maximum likelihood fit for the coefficient alpha for a distribution of 
    discrete values p(x) = alpha exp(-alpha*x) using only the largest N 
    values.
    '''
    vals = numpy.array(vals)
    if len(vals) < N:
        raise RuntimeError('Not enough input values to do fit!')
    vals.sort()
    vals = vals[-N:]
    return fitdict[distr](vals, min(vals))

def fit_fn(distr, xvals, alpha, thresh):
    '''
    The fitted function normalized to 1 above threshold
    '''
    return fndict[distr](xvals, alpha, thresh)

def cum_fit(distr, xvals, alpha, thresh):
    '''
    The integral of the fitted function above a given value (reverse CDF)
    normalized to 1 above threshold
    '''
    return cum_fndict[distr](xvals, alpha, thresh)

def KS_test(distr, vals, alpha, thresh=None):
    '''
    Perform KS test of the given set of discrete values above a given
    threshold for the fitted distribution function
    ex.: KS_test('exponential', vals, alpha, thresh)

    If no threshold is specified, the minimum sample value will be used.
    '''
    vals = numpy.array(vals)
    if thresh == None:
        thresh = min(vals)
    else:
        vals = vals[vals >= thresh]
    cdf_fn = lambda x: 1 - cum_fndict[distr](x, alpha, thresh)
    return kstest(vals, cdf_fn)

