from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse

def scale_spec(datafile):
    '''
    Returns a set of basic statistics for an SDSS FITS spectrum
    
    datafile = full path for the FITS file to be used
    
    Output:
    
    Number of pixels is spectrum,
    Overall averaged loglam channel width,
    Average flux,
    Standard Deviation of flux.
    
    '''
    hdu_list = fits.open(datafile)
    
    # Grab axis and flux data 
    dataset1 = hdu_list[1].data["flux"]
    dataset2 = hdu_list[1].data["loglam"]
    
    # Calculate stats
    numpix = len(list(dataset2))
    chanwidth = (dataset2[-1]-dataset2[0])/numpix
    avflux = np.average(dataset1)
    std_dev = np.std(dataset1)
    
    # Return calculated stats
    return(numpix,chanwidth,avflux,std_dev)