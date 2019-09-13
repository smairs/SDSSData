from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class Spectrum():
    def __init__(self,datafile):
        self._datafile = datafile

        # get flux information and harvest the header
        hdu_list      = fits.open(self._datafile)
        data          = hdu_list[1].data
        hd            = hdu_list[0].header
        self.loglam   = data['loglam']
        self.flux     = data['flux']        

        # extract header info
        self.hd       = hd
        self.ra       = hd['RA']
        self.dec      = hd['DEC']
        self.telesope = hd['TELESCOP']
        self.mjd      = hd['MJD']

    def get_spectrum(self):
        '''
        Get the wavelength and flux data of an SDSS spectrum
    
        loglam, flux = get_spectrum()
        '''
        return(self.loglam,self.flux)

    def find_lines(self,threshold,noise_range=[3.68,3.73],loglam_range=[3.60,3.96],display_lines=True):
        '''
        Find lines with peaks above a specified threshold (in emission only!)

        threshold    = Given in terms of multiples of the noise, how high the peak must be to be considered a line

        KWARGS:

        noise_range  = The range within which to measure the noise
        loglam_range = The range within which to look for interesting lines
        display_lines= Show a graphical representation of where the lines are


        Out:

        Locations of >=5*sigma detections

        Use:

        line_locations, peak_flux = find_lines(datafile)
        '''

        # First, get the data:
        loglam, flux = self.loglam, self.flux

        # Next, trim the edges:
        good_ind     =  np.where(np.logical_and(loglam>=loglam_range[0],loglam<=loglam_range[1]))
        loglam_good  = loglam[good_ind]
        flux_good    = flux[good_ind]

        # Now measure the noise:
        noise        = np.std(flux_good[np.where(np.logical_and(loglam_good>=noise_range[0],loglam_good<=noise_range[1]))],ddof=1)
        mean         = np.mean(flux_good[np.where(np.logical_and(loglam_good>=noise_range[0],loglam_good<=noise_range[1]))])

        peak_cutoff  = mean+(noise*threshold)

        # Now search for lines!
        line_locations     = []
        peak_fluxes        = []
        line_counter       = []
        line_counter_dummy = 0
        flux_dummy         = 0
        for eachflux,eachloglam in zip(flux_good,loglam_good):
            if eachflux>=peak_cutoff:
                if flux_dummy > 0 and flux_dummy < len(flux_good)-2:
                    if eachflux > flux_good[flux_dummy+1] and eachflux > flux_good[flux_dummy-1]:
                        line_counter.append(line_counter_dummy+1)
                        line_locations.append(eachloglam)
                        peak_fluxes.append(eachflux)
                        line_counter_dummy += 1
            flux_dummy += 1

        # Plot the spectrum with the line locations annotated
        if display_lines == True:
            plt.plot(loglam_good,flux_good)
            for eachline,eachloglam,eachpeakflux in zip(line_counter,line_locations,peak_fluxes):
                plt.annotate(str(eachline),(eachloglam,eachpeakflux))
            plt.xlabel('loglam')
            plt.ylabel('Flux')

        return(line_locations,peak_fluxes)        

# Statistical functions go here        
    @property
    def fluxav(self):
        return np.average(self.flux)

# We probably shouldn't allow them to set the average themselves?
#        @fluxav.setter
#        def fluxav(self,new_fluxav):
#            self.fluxav = new_fluxav 
    
    @property
    def fluxstd(self):
        return np.std(self.flux,ddof=1)
    
    @property
    def numpix(self):
        return len(self.flux)
        
    @property
    def chanwidth(self):
        chanwidth = (slef.loglam[-1]-self.loglam[0])/self.numpix
        return chanwidth


    # Gaussian fit a region!
    def gaussfit(self,loglam_range=[3.74,3.755],display=True):
        '''
        Fit a guassian to a range of data
        '''
        loglam_in_range = self.loglam[np.where(np.logical_and(self.loglam >= loglam_range[0],self.loglam <= loglam_range[1]))]
        flux_in_range   = self.flux[np.where(np.logical_and(self.loglam >= loglam_range[0],self.loglam <= loglam_range[1]))]
        n     = len(loglam_in_range)                   #the number of data
        mean  = sum(loglam_in_range*flux_in_range)/n              
        sigma = sum(flux_in_range*(loglam_in_range-mean)**2)/n
        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        popt,pcov = curve_fit(gaus,loglam_in_range,flux_in_range,p0=[1,mean,sigma])

        if display == True:
            plt.plot(loglam_in_range,flux_in_range,color='blue',linestyle='dashed',label='data')
            plt.plot(loglam_in_range,gaus(loglam_in_range,*popt),color='red',linestyle='dotted',label='fit')
            plt.legend()
            plt.title(f'Gaussian Fit for loglam range {loglam_range}')
            plt.xlabel('loglam')
            plt.ylabel('Flux')
            plt.show()

        return(popt,pcov)
    
    def get_std(self, l_min ,l_max):

        data = np.array([ x for x in zip(self.loglam, self.flux)], dtype=[('loglam', self.loglam.dtype), ('flux', self.flux.dtype)])
        #l_min=1
        #l_max=4

        return(np.std(data['flux'][(data['loglam'] > l_min) & (data['loglam'] < l_max)]))

    # 
