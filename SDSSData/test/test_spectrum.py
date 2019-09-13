#!/usr/bin/env python

# This file is for the purposes of testing the
# spectrum class file - ensuring each function works
# as expected

from ..spectrum import Spectrum

example_fits = 'spectrum.fits'

# Check to ensure header keys exist:

def test_pressure():
    a = Spectrum(example_fits)
    assert a.hd.get('PRESSURE') != None, "PRESSURE keyword doesn't exist!"




