# SDSSData
SciCoder-2019 Exercise to build a module for analysing SDSS Spectral Data

## Demo

Head into the demo directory and run "python demo.py" to see a 
demonstration of the Spectrum class and its various methods!

## Usage

from spectrum import Spectrum

my\_spectrum = Spectrum('Path/To/Data/File/')

### Quickly show spectrum

my\_spectrum.showme()

### returns the loglam and flux columns of the fits data

loglam,flux = my\_spectrum.get\_spectrum()

### Finds lines above a 5 sigma threshold:

line\_locations,peak\_fluxes = my\_spectrum.find\_lines(5)

#### Show lines above a 5 sigma threshold:

my\_spectrum.show\_lines(5)

#### Fit a line with a Gaussian model:

my\_spectrum.gaussfit([lower\_loglam,upper\_loglam])

### Get statistical properties:

my\_spectrum.fluxav

my\_spectrum.fluxstd

my\_spectrum.get\_std(lower\_loglam,upper\_loglam)

### Get number of channels and channel width in the spectrum

my\_spectrum.numpix

my\_spectrum.chanwidth

### Print a summary of the weather

my\_spectrum.weathersum()

