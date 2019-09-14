#!/usr/bin/python env

import sys
sys.path.append('../')
from spectrum import Spectrum

example_spectrum = '../test/spectrum.fits'

a = Spectrum(example_spectrum)


###
# Begin

print('\n\nThis is a demo for the Spectrum class.')

print('\n\n###################')
print('# Initialisation  #')
print('###################')

print(f'\nThe SDSS spectrum fits file "{example_spectrum}" \n'+ \
          'has been loaded into the Spectrum class and \n'+ \
          'initiated using the variable "a":\n\n'+ \
          f'>>>a = Spectrum("{example_spectrum}")')

print('\nSeveral variables from the header are already \n'+ \
         'initialised and can be called as follows:\n'+ \
         '>>>a.ra\n'+ \
         '>>>a.dec\n'+ \
         '>>>a.mjd')

print('')
print(f'a.ra  = {a.ra}')
print(f'a.dec = {a.dec}')
print(f'a.mjd = {a.mjd}')
print('')

print('\n\n--------------------\n\n')

input('Press Enter to continue...\n\n')

###
# Weather Summary

print('###################')
print('# Weather Summary #')
print('###################')

print('\nYou can also get a summary of the weather! \n\n' + \
         '>>>a.weathersum()')

print(a.weathersum())
print('^Prints "None" as it is a function that returns nothing')

print('')

print('\n\n--------------------\n\n')

input('Press Enter to continue....')

###
# Properties that are really methods

print('\n\n###################')
print('#   Statistics    #')
print('###################')

print('\nDecorators are used to turn methods into properties: \n' + \
         'Average Flux, Standard Dev, Numpix (Number of data points), and Chanwidth (Channel Width): \n')

print('>>>a.fluxav\n'+ \
      '>>>a.fluxstd\n' + \
      '>>>a.numpix\n' + \
      '>>>a.chanwidth\n')

print(f'a.fluxav    = {a.fluxav}')
print(f'a.fluxstd   = {a.fluxstd}')
print(f'a.numpix    = {a.numpix}')
print(f'a.chanwidth = {a.chanwidth}')

print('')

print('You can also specify a loglam range over which to measure the standard deviation:\n')

print('>>>a.get_std(3.6,3.9)')

print(a.get_std(3.6,3.9))

print('')

print('\n\n--------------------\n\n')

input('Press Enter to continue...')

###
# showme!

print('\n\n###################')
print('#   Quick Show    #')
print('###################')

print('')

print('What if you just want to see the spectrum?\n')

print('>>>a.showme()')

input('\nPress Enter to run a.showme()....')

a.showme()

print('\n\n--------------------\n\n')

input('Press Enter to continue...')

###
# get_spectrum

print('\n\n###################')
print('#Get Loglam & Flux#')
print('###################')


print('\nSave the loglam and flux arrays into variables:\n')

print('>>>loglams,fluxes = a.get_spectrum()')

print('\n\n--------------------\n\n')

input('Press Enter to continue...')

###
# find_lines and show_lines

print('\n\n###################')
print('#Find & Show Lines#')
print('###################')

print('\nUse find_lines and show_lines to identify peaks above a supplied SNR threshold:')

print('\n>>>line_locations,peak_fluxes = a.find_lines(SNR)')
print('\n>>>a.show_lines(SNR)')

input('\nPress Enter to run a.show_lines(5)...')

a.show_lines(5)

print('')

print('\n\n--------------------\n\n')

input('Press Enter to continue...')

###
# Gaussian fits

print('\n\n###################')
print('#Gaussian Fit Area#')
print('###################')

print('\nYou can also perform a Gaussian fit over a specified range:')

print('\n>>>a.gaussfit([loglam_min,loglam_max])')

input('\nPress Enter to run a.gaussfit(3.74,3.755)...')

a.gaussfit()

print("\n\nThat's all for now! Have a great day!\n\n")
