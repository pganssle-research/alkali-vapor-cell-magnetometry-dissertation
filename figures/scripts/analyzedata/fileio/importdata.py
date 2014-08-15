# Script for importing data into console.
def importdata(fpath=''):
    import scipy.io
    import os
    import readhistory;
    from Tkinter import Tk
    from tkFileDialog import askopenfilename
    
    histdir = os.path.dirname(__file__)
    histname = 'load_hist.pkl'
    
    if not os.path.exists(histdir):
        histdir = os.getcwd()
    
    histpath = os.path.join(os.path.dirname(histdir), histname);
    
    # If it's not there, use CWD.
    hr = readhistory.histReader(histpath, 20)
    prevdir = os.getcwd()
    
    for p in hr.histpaths:
        if os.path.exists(p):
            prevdir = p
            break
    
    if fpath != '' and os.path.exists(fpath):
        filename = fpath
    else:
        Tk().withdraw();   # No root window.    
        filename = askopenfilename(defaultextension='.mat', initialdir=prevdir)
        
    fe = os.path.splitext(filename)[1]
    if fe.lower() != '.mat':
        raise Exception('Invalid file extension' + fe)
    else:
        out = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
        
        hr.addpath(os.path.dirname(os.path.abspath(filename)))
        hr.update()
        
        return out

