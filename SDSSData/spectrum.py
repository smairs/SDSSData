from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

class Spectrum():
    def __init__(self):
        self.datafile = None

    def get_spectrum(self):
        '''
        Get the wavelength and flux data of an SDSS spectrum
    
        datafile = full path to the FITS file you want to get the spectrum for
    
        Output:
    
        The loglam and flux columns of the fits file.
    
        Use:
    
        loglam, flux = get_spectrum(datafile)
        '''
        hdu_list = fits.open(self.datafile)
        data     = hdu_list[1].data
        hd       = hdu_list[1].header
        loglam   = data['loglam']
        flux     = data['flux']
        return(loglam,flux)

    def find_lines(self,threshold,noise_range=[3.68,3.73],loglam_range=[3.60,3.96],display_lines=True):
        '''
        Find lines with peaks above a specified threshold (in emission only!)

        datafile = Full path to SDSS FITS spectrum

        kwargs:

        threshold    = Given in terms of the noise, how high the peak must be to be considered a line
        noise_range  = The range within which to measure the noise
        loglam_range = The range within which to look for interesting lines
        display_lines= Show a graphical representation of where the lines are


        Out:

        Locations of >=5*sigma detections

        Use:

        line_locations, peak_flux = find_lines(datafile)
        '''

        # First, get the data:
        loglam, flux = self.get_spectrum()

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
            wave,flux = self.get_spectrum()
            return np.average(flux)
    	
    	@property
        def fluxstd(self):
            wave,flux = self.get_spectrum()
            return np.std(flux,ddof=1)
    	
    	@property
        def numpix(self):
            wave,flux = self.get_spectrum()
            return len(wave)
        
        @property
        def chanwidth(self):
            wave,flux = self.get_spectrum()
            chanwidth = (wave[-1]-wave[0])/self.numpix
            return chanwidth

    def fit_lines(self):
        pass
