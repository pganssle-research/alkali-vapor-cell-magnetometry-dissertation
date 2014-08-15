from numpy import *
from scipy import *
from scipy import optimize
from scipy.io import loadmat, savemat
from copy import copy

####################################################################################################
##                                                                                                ##
##                                       Math Utilities                                           ##
##                                                                                                ##
####################################################################################################

def moving_avg(x, navg, mode='box', window_func=None):
    """
    Calculate the moving average of dataset x. This is an in-place transform. The value returned is
    the average of +/- navg points on either side. On the boundaries, the number of averages is 
    reduced to the largest number les than navg. So the first point is unchanged, the second point 
    is the average of points 1-3, etc.
    
    This is a wrapper for a few other methods, moving_avg_linear, moving_avg_cubic, etc. The default
    value is linear, but it may be valuable to use the other methods directly when embedding into a
    large (potentially slow) data-processing application.
    
    @param x    The input - any sort of vector.
    @param navg Number of averages to include in the window.
    @param mode The equation used for the windowing function (used to weight the average). Default
                state is 'box', other options are:
                    'box':          A simple average
                    'triangle':     A triangle window
                    'linear':       Windowing function is 1/abs(i)
                    'cubic':        Windowing function is navg**2/(abs(i)**3)
                    'Hann':         Windowing function is a Hann window.
                    'custom':       Choose any windowing function, it will be spline interpolated
                                    into navg points. This should be passed to window_func.
    
    @return Returns a smoothed curve of the same size as the input.
    @rtype  array
    """
       
    if mode == 'box':
        return moving_avg_box(x, navg);
    else:
        raise UnsupportedError('Mode not yet supported.')
    
def moving_avg_box(x, navg):
    """
    Calculate the moving average of dataset x using a boxtop filter.
    
    @param x    The vector to smooth.
    @param navg The number of points on either side of the window.
    
    @return Removes a smoothed vector.
    @rtype array   
    """
    out = zeros(shape(x))
    np = len(x)
    
    for ii in range(0, navg):
        tna = min([2*(ii+1), np])
        out[ii] = mean(x[0:tna])
        out[-(ii+1)] = mean(x[-tna:])
    
    for ii in range(navg, np-navg):
        out[ii] = mean(x[ii-navg:ii+navg+1])

    return out
    
def moving_avg_norm_window(x, win_fun):
    """
    Calculate the moving average of a dataset using a normalized window function. This will likely 
    be called by various helpers.
    
    @param x        The vector to smooth
    @type  x        Numpy array (numpy.array)
    @param win_fun  The windowing function, must be normalized with an odd number of elements.
    @type win_fun   For now, a normalized numpy array.
  
    @return Returns a smoothed vector.
    """
    navg = len(win_fun/2)
    np = len(x);
    
    for ii in range(0, navg-1):
        tna = min([2*ii, np])
        out[ii] = mean(x[0:tna]*win_fun[navg-ii:navg+ii])
        out[-ii] = mean(x[-tna:]*win_fun[navg-ii:navg+ii])        
        
    for ii in range(navg, np-navg):
        out[ii] = mean(x[ii-navg:ii+navg]*win_fun)
        
    return out

####################################################################################################
##                                                                                                ##
##                                    Searching Utilities                                         ##
##                                                                                                ##
####################################################################################################
def find_where(condition, arr=None, n=None, order='first', ignoreNaN=False):
    """
    Equivalent to a MATLAB-style find(), this returns n arguments where 'condition' is met, either 
    starting from the beginning (order='first') or the end (order='last').
    
    Example - Find the first element of 'x' which equals 3
    x = [0, 1, 4, 6, 3, 9, 7, 3]
    fi = find_where(lambda x: x == 3, x, 1, 'first')
    print 'fi == {}'.format(fi)
    > fi == 4
    
    
    @param condition    A condition imposed on the array. This can either be a function which is
                        executed for each value in the array.
    @type condition     array_like, function
    
    @param arr          An array which is being tested. If condition is an array of logical values, 
                        then arr does not need to be included. Default: None
    @type arr           array_like
    
    @param n            The maximum number of indices to return. Pass None to return all matching 
                        values, otherwise the first n values are chosen. Default: None
    @type n             integer
    
    @param order        The order in which to execute the search. There are two valid states:
                            'first':    Start from index 0 going to the end.
                            'last':     Start from index -1 going to 0.
                            
                        If 'last' is chosen, the vector will be returned in the order it originally
                        appeared in the vector. For example: 
                        
                            find_where(lambda x: x >= 3, [0, 1, 2, 3, 4, 5, 6], n=3, order='last')
                        
                        returns:
                        
                            [4, 5, 6]
    @type   order       string
    
    @param ignoreNaN    Boolean value of whether or not to explicitly ignore NaNs. If condition is
                        a function, this will only check of arr is NaN, not if the result is NaN.
                        If condition is an array, it will be checked for being NaN. Default: False
    @type  ignoreNaN    bool
    
    @return Returns a list of indices or None if no indices are found.
    @rtype  array
    
    @raises ValueError  Raised for invalid inputs.
    """
    
    # If condition is just an array, then the arr parameter is not necessary. In the other case, 
    # we're looking at a function with an argument, so we will convert the condition array probing
    # to a function rather than the other way around (which requires more overhead). Both will be
    # called by the index.
    if not hasattr(condition, '__call__'):
        cond_arr = condition
        if ignoreNaN:
            condition = lambda x: cond_arr[x] and cond_arr[x] is not NaN
        else:
            condition = lambda x: cond_arr[x]
    else:
        cond_fun = condition
        if ignoreNaN:
            condition = lambda x: arr[x] is not NaN and cond_fun(arr[x])
        else:
            condition = lambda x: cond_fun(arr[x])
    
    # Establish local parameters in convenience variables.
    if n == None or n < 1:
        n = len(arr)
 
    isFirst = order == 'first'

    nf = 0
    found = []
    
    # Establish search range
    if isFirst:
        sr = arange(0, len(arr), 1)
    else:
        sr = arange(len(arr), 0, -1) 
        sr = sr-1                    # Without this you'd stop at 1.
        
    # Execute the search
    for i in sr:
        if condition(i):
            nf += 1
            found.append(i)
        
        if nf >= n:
            break
        
    # Return the values as appropriate.
    if nf == 0:
        return None
    elif nf == 1:
        return found[0]
        
    if isFirst:
        return found
    else:
        return [k for k in reversed(found)] # So it appears as it did in the original vector.
        
    
    
####################################################################################################
##                                                                                                ##
##                                     Plotting Utilities                                         ##
##                                                                                                ##
####################################################################################################

def squarest_grid(n, prefer_columns=True):
    """
    Gives the pair of integer factors of n closest to sqrt(n). Designed for automatically finding 
    the optimal subplot to give the squarest grid, i.e. for n=4, return (2, 2), for n=5, return
    (5, 1), for n=6 return (3,2)
    
    @param n    The total number of items in the grid.
    @type n int
    
    @param prefer_columns For non-square grids, setting this to true gives the larger number first,
    setting it to false gives the smaller number first. Default: True
    @type prefer_columns bool
    
    @return Returns a 2-tuple of ints which represent the number of columns and number of rows, 
    respectively, in the grid.
    @rtype (int, int)
    """
    nroot = sqrt(n)
    nrootrnd = int(ceil(nroot))

    [a, b] = [n, 1]     # Initialize to something. Probably not necessary.
    
    for ii in range(nrootrnd, n+1):        # Find the nearest even multiple
        if n % ii == 0:
            a = ii
            b = n/a
            break
            
    if prefer_columns:
        return (a,b)    # a is always the larger number.
    else:
        return (b,a)

def fix_label_ticks(ax, fun=(lambda x: 10**x), fix_x=True):
    # Function for transforming label ticks according to the function fun
    
    if fix_x:
        iter_list = ax.get_xticklabels();
        nlf = lambda x: ax.set_xticklabels(x)
    else:
        iter_list = ax.get_yticklabels();
        nlf = lambda x: ax.set_yticklabels(x)
    
    ff = lambda x: '{{:0.{:d}g}}'.format(max([int(floor(log10(x)))+1, 2]) if x > 0 else 1).format(x)
    
    labels = []
    for item in iter_list:
        citem = (item.get_text()).replace(u'\u2212', '-')
           
        try:
            citem = ff(fun(float(citem)))
        except:
            pass
        
        labels.append(citem)
        
    nlf(labels)

def min_string(x, pos=None):
    # Calculate the minimal string
    return '{{:0.{:d}g}}'.format(max([int(floor(log10(x)))+1, 2]) if x > 0 else 1).format(x)
        
####################################################################################################
##                                                                                                ##
##                                      Fitting Utilities                                         ##
##                                                                                                ##
####################################################################################################
def exp_sq_fit(xdata, params):
    """
    Kernel for square exponential fit without an offset. Returns:
    
    params[1]*exp(-(xdata**2)*params[0]);
    
    So params has the form: [Time Constant, Amplitude]
    
    @param  xdata   The x data values.
    @type   xdata   numpy.array
    
    @param  params  A 2-vector of the form [Time Constant, Amplitude]
    @type   params  numpy.array
    
    @return Returns the ydata described by the parameters.
    @rtype  numpy.array
    """
    
    return params[1]*exp(-(xdata**2)*params[0]);
        
def exp_sq_fit_o(xdata, params):
    """
    Kernel for square exponential fit with an offset. Returns:
    
    params[1]*exp(-(xdata**2)*params[0]) + params[2];
    
    So params has the form: [Time Constant, Amplitude]
    
    @param  xdata   The x data values.
    @type   xdata   numpy.array
    
    @param  params  A 3-vector of the form [Time Constant, Amplitude, Offset]
    @type   params  numpy.array
    
    @return Returns the ydata described by the parameters.
    @rtype  numpy.array
    """
    
    return params[1]*exp(-(xdata**2)*params[0]) + params[2];

def exp_fit(xdata, params):
    """
    Kernel for exponential fit without an offset. Returns:
    
    params[1]*exp(-xdata/params[0]);
    
    So params has the form: [Time Constant, Amplitude]
    
    @param  xdata   The x data values.
    @type   xdata   numpy.array
    
    @param  params  A 2-vector of the form [Time Constant, Amplitude, Offset]
    @type   params  numpy.array
    
    @return Returns the ydata described by the parameters.
    @rtype  numpy.array
    """
    
    return params[1]*exp(-xdata/params[0])
    
def exp_fit_o(xdata, params):
    """
    Kernel for exponential fit with an offset. Returns:
    
    params[1]*exp(-xdata/params[0]) - params[2];
    
    So params has the form: [Time Constant, Amplitude, Offset]
    
    @param  xdata   The x data values.
    @type   xdata   numpy.array
    
    @param  params  A 3-vector of the form [Time Constant, Amplitude, Offset]
    @type   params  numpy.array
    
    @return Returns the ydata described by the parameters.
    @rtype  numpy.array
    """
    
    return params[1]*exp(-xdata/params[0]) - params[2]

def sin_fit_p(xdata, params):
    """
    Kernel for sine function with phase. Returns:
    
    params[2]*sin(params[0]*xdata - params[1]);
    
    So params has the form: [Frequency, Phase, Offset]
    
    @param  xdata   The x data values.
    @type   xdata   numpy.array
    
    @param  params  A 3-vector of the form [Time Constant, Amplitude, Offset]
    @type   params  numpy.array
    
    @return Returns the ydata described by the parameters.
    @rtype  numpy.array
    """
    return params[2]*sin(multiply(xdata, params[0]) + params[1]);

def fit_g(x, y, guess=None):
    """
    Convenience method for fitting Gaussians with an offset.
    
    @param  x       The x values. Vector. Should be same size as y.
    @type   x       numpy.array
    
    @param  y       The y values. Vector. Should be same size as x.
    @type   y       numpy.array
    
    @param  guess   A best guess for the initial value of T2. Of the form [sigma, GAmp, GOff]. 
                    Pass None to estimate from data. Default: None
    @type   guess   numpy.array
    
    @return Returns an array of the form [T2, T2Amp, T2Off]. This can be passed directly to the 
            function exp_fit_o;
    @rtype  numpy.array
    """
    
    if guess == None:
        gai = max(y)
        ggi = argmin(abs(y**2-gai*exp(-1)))
        
        # If we don't see enough decay, select an index halfway between beginning and 
        # end and calculate it from that
        if ggi == len(x):
            cind = int(round(len(x)/2))
            ggi = log(gai/y[cind])/(x[cind]**2)
        else:
            ggi = 1/x[ggi]     
        
        guess = [ggi, gai]
    
    [output, e] = optimize.leastsq(lambda params: exp_sq_fit(x, params)-y, guess);
    
    if e > 4 or e < 1:
        print('Error in fit.')
    
    return output
def fit_go(x, y, guess=None):
    """
    Convenience method for fitting Gaussians with an offset.
    
    @param  x       The x values. Vector. Should be same size as y.
    @type   x       numpy.array
    
    @param  y       The y values. Vector. Should be same size as x.
    @type   y       numpy.array
    
    @param  guess   A best guess for the initial value of T2. Of the form [sigma, GAmp, GOff]. 
                    Pass None to estimate from data. Default: None
    @type   guess   numpy.array
    
    @return Returns an array of the form [T2, T2Amp, T2Off]. This can be passed directly to the 
            function exp_fit_o;
    @rtype  numpy.array
    """
    
    if guess == None:
        gai = max(y)
        ggi = 1/argmin(abs(y**2-gai*exp(-1)))
        
        # Assuming offset is zero.
        goi = 0;
        
        # If we don't see enough decay, select an index halfway between beginning and 
        # end and calculate it from that
        if ggi == len(x):
            cind = int(round(len(x)/2))
            ggi = log(gai/y[cind])/(x[cind]**2)
        else:
            ggi = 1/x[ggi]
        
        guess = [ggi, gai, goi]
    
    [output, e] = optimize.leastsq(lambda params: exp_sq_fit_o(x, params)-y, guess);
    
    if e > 4 or e < 1:
        print('Error in fit.')
    
    return output
    
def fit_t2o(x, y, guess=None):
    """
    Convenience method for fitting T2.
    
    @param  x       The x values. Vector. Should be same size as y.
    @type   x       numpy.array
    
    @param  y       The y values. Vector. Should be same size as x.
    @type   y       numpy.array
    
    @param  guess   A best guess for the initial value of T2. Of the form [T2, T2Amp, T2Off]. 
                    Pass None to estimate from data. Default: None
    @type   guess   numpy.array
    
    @return Returns an array of the form [T2, T2Amp, T2Off]. This can be passed directly to the 
            function exp_fit_o;
    @rtype  numpy.array
    """
    
    if guess == None:
        t2ai = max(y)
        t2gi = argmin(abs(y-t2ai*exp(-1)))
        
        # Assuming offset is zero.
        t2oi = 0;
        
        # If we don't see enough decay, select an index halfway between beginning and 
        # end and calculate it from that
        if t2gi == len(x):
            cind = int(round(len(x)/2))
            t2gi = x[cind]/(log(t2ai/y[cind]))
        else:
            t2gi = x[t2gi];
            
        
        guess = [t2gi, t2ai, t2oi]
    
    [output, e] = optimize.leastsq(lambda params: exp_fit_o(x, params)-y, guess);
    
    if e > 4 or e < 1:
        print('Error in fit.')
    
    return output
    
def fit_t2(x, y, guess=None):
    """
    Convenience method for fitting T2.
    
    @param  x       The x values. Vector. Should be same size as y.
    @type   x       numpy.array
    
    @param  y       The y values. Vector. Should be same size as x.
    @type   y       numpy.array
    
    @param  guess   A best guess for the initial value of T2. Of the form [T2, T2Amp, T2Off]. 
                    Pass None to estimate from data. Default: None
    @type   guess   numpy.array
    
    @return Returns an array of the form [T2, T2Amp, T2Off]. This can be passed directly to the 
            function exp_fit_o;
    @rtype  numpy.array
    """
    
    if guess == None:
        t2ai = max(y)
        t2gi = argmin(abs(y-t2ai*exp(-1)))
                
        # If we don't see enough decay, select an index halfway between beginning and 
        # end and calculate it from that
        if t2gi == len(x):
            cind = int(round(len(x)/2))
            t2gi = -x[cind]/(log(y[cind]/log(t2ai)))
        else:
            t2gi = x[t2gi]
            
        guess = [t2gi, t2ai]
    
    [output, e] = optimize.leastsq(lambda params: exp_fit(x, params)-y, guess);
    
    if e > 4 or e < 1:
        print('Error in fit.')
    
    return output
    
def fit_exp_and_g_o(x, y, guess=None):
    """
    Convenience method for fitting data to a linear combination of a gaussian and an exponential.
    
    Fits to:
    
    exp_fit(x, y, [params[0], params[1]]) + exp_sq_fit(x, y, [params[2], params[3]]) + params[4]
    
    @param  x       The x values. Vector. Should be same size as y.
    @type   x       numpy.array
    
    @param  y       The y values. Vector. Should be same size as x.
    @type   y       numpy.array
    
    @param  guess   A best guess for the initial value of the result.
    @type   guess   numpy.array
    
    @return An array with the resulting params.
    @rtype  numpy.array
    """
    
    if guess == None:
        t2ai = max(y)
        t2gi = argmin(abs(y-t2ai*exp(-1)))
        
        gai = max(y)
        ggi = argmin(abs(y**2-gai*exp(-1)))
                
        # If we don't see enough decay, select an index halfway between beginning and 
        # end and calculate it from that
        if t2gi == len(x):
            cind = int(round(len(x)/2))
            t2gi = -x[cind]/(log(y[cind]/log(t2ai)))
        else:
            t2gi = x[t2gi]
            
        if ggi == len(x):
            cind = int(round(len(x)/2))
            ggi = log(gai/y[cind])/(x[cind]**2)
        else:
            ggi = 1/x[ggi]
            
        guess = [t2gi/2, t2ai, ggi/2, gai, 0]
    
    [output, e] = optimize.leastsq(lambda params: exp_fit(x, [params[0], params[1]]) + exp_sq_fit(x, [params[2], params[3]]) + params[4]-y, guess);
    
    if e > 4 or e < 1:
        print('Error in fit.')
    
    return output
 
def get_nan_inds(ar):
    """
    Gets a list of all indices which are nan for a 1-dimensional array
    """
    
    ain = []
    for i in range(0, len(ar)):
        if isnan(ar[i]):
            ain.append(i)
    
    return ain

def get_ian_inds(ar):
    """
    Gets a list of all indices which are NOT nan for a 1-dimensional array
    """
    
    ain = []
    for i in range(0, len(ar)):
        if not isnan(ar[i]):
            ain.append(i)
    
    return ain
 
def maxf(ar, max_val=nan):
    """
    Gets the maximum value in the entire N-D array.

    @param ar  The array.
    """

    sa = shape(ar)
    np = 1
    for n in sa:
        np *= n
        
    ar2 = reshape(ar, np)
    
    ar2 = delete(ar2, get_nan_inds(ar2), 0)
    
    cinds = []
    if not isnan(max_val):
        for ii in range(0, len(ar2)):
            if ar2[ii] >= max_val:
                cinds.append(ii)
    
    if len(cinds) > 0:
        ar2 = delete(ar2, cinds, 0)
    
    if size(ar2) == 0:
        return nan
       
    return max(ar2)
        
def minf(ar, min_val=nan):
    """
    Gets the minimum value in the entire N-D array.

    @param ar  The array.
    """

    sa = shape(ar)
    np = 1
    for n in sa:
        np *= n
        
    ar2 = reshape(ar, np)
    
    ar2 = delete(ar2, get_nan_inds(ar2), 0)
    
    cinds = []
    if not isnan(min_val):
        for ii in range(0, len(ar2)):
            if ar2[ii] <= min_val:
                cinds.append(ii)
    
    if len(cinds) > 0:
        ar2 = delete(ar2, cinds, 0)
    
    if size(ar2) == 0:
        return nan
        
    return min(ar2)

 

def mina(ar, dim=0, nonzero=False):
    """
    Gets the minimum along a given dimension (0-based index)
    
    @param ar The array
    @param dim The dimension (0-based index, default: 0)
    """
    
    sa = shape(ar)
    
    sas = delete(sa, dim) 
    oar = zeros(sas)
    np = 1
    for n in sas:
        np *= n
    
    ds = sa[dim]
    
    for ii in range(0, np):
        ind = unravel_index(ii, dims=sas)
        
        ind2 = list(copy(ind))
        ind2.insert(dim, range(0, ds))        
            
        oar[ind] = min(ar[ind2])
    
    return oar


        
 