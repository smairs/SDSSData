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

def test_airtemp():
    a = Spectrum(example_fits)
    assert a.hd.get('AIRTEMP') != None, "AIRTEMP keyword doesn't exist!"

def test_humidity():
    a = Spectrum(example_fits)
    assert a.hd.get('HUMIDITY') != None, "HUMIDITY keyword doesn't exist!"

def test_winds():
    a = Spectrum(example_fits)
    assert a.hd.get('WINDS') != None, "WINDS keyword doesn't exist!"

def test_windd():
    a = Spectrum(example_fits)
    assert a.hd.get('WINDD') != None, "WINDD keyword doesn't exist!"

