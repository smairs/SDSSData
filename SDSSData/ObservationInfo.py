# By: li-jcmt
from astropy.io import fits

class Observation():
    '''
    Create an observation class that can be used to retrieve basic observation info.
    3 info added for now:
        1) telescope name
        2) RA
        3) DEC
    '''
    def __init__(self, filename=None):
        if filename is None:
            raise FileNotSpecified("A spectrum file must be specified")
        self.filename = filename
        hdu_list = fits.open(self.filename)
        self.hdr = hdu_list[0].header
        hdu_list.close()

    @property
    def telescope(self):
        return self.hdr['TELESCOP']

    @property
    def ra(self):
        return (self.hdr['RA'])

    @property
    def dec(self):
        return (self.hdr['DEC'])



