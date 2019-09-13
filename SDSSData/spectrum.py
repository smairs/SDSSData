from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class Spectrum():
    def __init__(self,datafile):
        
        # Keep the init method light,
        # could define a deparate "open"
        # method - but if we are loading in
        # a spectrum, we will probably
        # always want all of this information
        # immediately handy

        self._datafile = datafile

        # Extract primary data and header
        # There is no data in hdu_list[0].data
        hdu_list      = fits.open(self._datafile)
        data          = hdu_list[1].data
        hd            = hdu_list[0].header

        # Get flux and wavelength information
        self.loglam   = data['loglam']
        self.flux     = data['flux']        

        # Make full header available and extract other info
        self.hd       = hd
        self.ra       = hd['RA']
        self.dec      = hd['DEC']
        self.telesope = hd['TELESCOP']
        self.mjd      = hd['MJD']
        self.plate    = hd['NAME']
        self.survey   = hd['PLATETYP']
        self.exptime  = hd['EXPTIME']
        
        # Weather Keywords:
        self.pressure = hd['PRESSURE']
        self.windspeed= hd['WINDS']
        self.winddir  = hd['WINDD']
        self.airtemp  = hd['AIRTEMP']
        self.humidity = hd['HUMIDITY']


###########
# Quick ways to grab the data and view it
###########

    # By: smairs
    def get_spectrum(self):
        '''
        Get the wavelength and flux data of an SDSS spectrum
    
        loglam, flux = Spectrum.get_spectrum()
        '''
        return(self.loglam,self.flux)

    # By: smairs
    def showme(self):
        '''
        Quick, non-customisable view of the spectrum

        Spectrum.showme()
        '''
        plt.plot(self.loglam,self.flux)
        plt.ylabel('Flux')
        plt.xlabel('Loglam')
        plt.show()

###########
# Quick statistical methods/properties: By: mgrawlings1
###########

    @property
    def fluxav(self):
        return np.average(self.flux)
    
    @property
    def fluxstd(self):
        return np.std(self.flux,ddof=1)

    # We may want to set our own noise level
    @fluxstd.setter
    def fluxstd(self,new_std):
        self.fluxstd = new_std

    @property
    def numpix(self):
        return len(self.flux)
        
    @property
    def chanwidth(self):
        chanwidth = (slef.loglam[-1]-self.loglam[0])/self.numpix
        return chanwidth

    # By: izumizuno
    # get SD in a given loglam range
    def get_std(self, l_min ,l_max):

        data = np.array([ x for x in zip(self.loglam, self.flux)], dtype=[('loglam', self.loglam.dtype), ('flux', self.flux.dtype)])
        #l_min=1
        #l_max=4

        return(np.std(data['flux'][(data['loglam'] > l_min) & (data['loglam'] < l_max)]))

############
# Methods to Summarise Header Information
############

    # By: smairs
    def weathersum(self):
        '''
        Print a Report summarising the weather

        Spectrum.weathersum()
        '''
        print('\nSummary of the Weather:')
        print('-----------------------\n')
        print(f"Bolometric Pressure = {self.pressure}")
        print(f"Wind Speed          = {self.windspeed}")
        print(f"Wind Direction      = {self.winddir}")
        print(f"Air Temperature     = {self.airtemp}")
        print(f"Humidity            = {self.humidity}\n")

#############
# Methods to Find and analyse lines
#############

    # By: smairs
    # Find lines above a given threshold in the
    # dumbest way possible
    def find_lines(self,threshold,noise_range=[3.68,3.73],loglam_range=[3.60,3.96],display_lines=False):
        '''
        Find lines with peaks above a specified threshold (in emission only!)

        threshold    = Given in terms of multiples of the noise (integer), 
                       how bright the peak must be to be considered a line

        KWARGS:

        noise_range  = The loglam range within which to measure the noise (list, int or float)
        loglam_range = The loglam range within which to look for interesting lines (list, int or float)
        display_lines= Show a graphical representation of where the lines are (Bool)

        Out:

        Locations of >=threshold detections

        Use:

        line_locations, peak_flux = Spectrum.find_lines(threshold)
        '''

        # First, get the data:
        loglam, flux = self.loglam, self.flux

        # Next, trim the edges:
        good_ind     = (loglam >= loglam_range[0]) & \
                       (loglam <= loglam_range[1])
        loglam_good  = loglam[good_ind]
        flux_good    = flux[good_ind]

        # Now measure the noise and get the mean flux in that range:
        noise_range_index = (loglam_good >= noise_range[0]) & \
                            (loglam_good <= noise_range[1])

        noise        = np.std(flux_good[noise_range_index],ddof=1)
        mean         = np.mean(flux_good[noise_range_index])

        # Define the peak_cutoff based on the noise
        # and user supplied threshold
        peak_cutoff  = mean+(noise*threshold)

        # Now search for lines! 
        # Just do the dumbest thing possible
        line_locations     = []
        peak_fluxes        = []
        line_counter       = []
        line_counter_dummy = 0
        flux_dummy         = 0
        for eachflux,eachloglam in zip(flux_good,loglam_good):
            if eachflux >= peak_cutoff:
                if flux_dummy > 0 and flux_dummy < len(flux_good)-2:
                    if eachflux > flux_good[flux_dummy+1] and eachflux > flux_good[flux_dummy-1]:
                        line_counter_dummy+=1
                        line_counter.append(line_counter_dummy)
                        line_locations.append(eachloglam)
                        peak_fluxes.append(eachflux)
            flux_dummy += 1

        # Plot the spectrum with the line locations annotated
        if display_lines == True:
            plt.plot(loglam_good,flux_good)
            for eachline,eachloglam,eachpeakflux in zip(line_counter,line_locations,peak_fluxes):
                plt.annotate(str(eachline),(eachloglam,eachpeakflux))
            plt.xlabel('loglam')
            plt.ylabel('Flux')
            plt.show()

        return(line_locations,peak_fluxes)        

    # By: smairs
    # Show lines detected above a gien threshold
    def show_lines(self,thresh):
        '''
        Show the lines above some SNR threshold
        This does not include noisy edges
        '''
        line_locations,peak_fluxes = self.find_lines(thresh)
        plt.plot(self.loglam,self.flux)
        line_number = 0
        for eachloglam,eachflux in zip(line_locations,peak_fluxes):
            line_number += 1
            plt.annotate(str(line_number),(eachloglam,eachflux))
        plt.xlabel('Loglam')
        plt.ylabel('Flux')
        plt.show()

    # By: smairs
    # Gaussian fit a region of the data!
    def gaussfit(self,loglam_range=[3.74,3.755],display=True):
        '''
        Fit a gaussian to a range of data

        amp,cen_lam,sigma,pcov = gaussfit()
        '''
        # Ignore everything outside the specified range
        range_ind = (self.loglam >= loglam_range[0]) & \
                    (self.loglam <= loglam_range[1])


        loglam_in_range = self.loglam[range_ind]
        flux_in_range   = self.flux[range_ind]

        # Perform the fit
        n     = len(loglam_in_range)                   #the number of data
        mean  = sum(loglam_in_range*flux_in_range)/n              
        sigma = sum(flux_in_range*(loglam_in_range-mean)**2)/n

        # Define a quick Guassian Function
        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        popt,pcov = curve_fit(gaus,loglam_in_range,flux_in_range,p0=[1,mean,sigma])

        # Turn the output into more useful variable names
        # For return
        amp,cen_lam,sigma = popt[0],popt[1],popt[2]

        if display == True:
            plt.plot(loglam_in_range,flux_in_range,color='blue',linestyle='dashed',label='data')
            plt.plot(loglam_in_range,gaus(loglam_in_range,*popt),color='red',linestyle='dotted',label='fit')
            plt.legend()
            plt.title(f'Gaussian Fit for loglam range {loglam_range}')
            plt.xlabel('loglam')
            plt.ylabel('Flux')
            plt.show()

        return(amp,cen_lam,sigma,pcov)