# SDSSData
SciCoder-2019 Exercise to Analyse SDSS Spectral Data:

## Usage

from spectrum import Spectrum

my\_spectrum = Spectrum('Path/To/Data/File/')

### returns the loglam and flux columns of the fits data

loglam,flux = my\_spectrum.get\_spectrum()

### Finds lines above a 5 sigma threshold:

line\_locations,peak\_fluxes = my\_spectrum.find\_lines(5)

### Get statistical properties:

my\_spectrum.fluxav

my\_spectrum.fluxstd


### Get number of channels in the spectrum

my\_spectrum.numpix
