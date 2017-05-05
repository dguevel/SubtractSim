from astropy.io import fits
from astropy import wcs
from astropy import coordinates
from astropy import units
import numpy as np
from PyZOGY.subtract import ImageClass
from PyZOGY.util import solve_iteratively
import os


def degs2coords(ra, dec, header):
    """Converts RA and DEC into coordinates on image; taken from ryan's code"""

    w = wcs.WCS(header)
    return w.all_world2pix(ra, dec, 0)


def find_flux_ratio(image, image_psf, reference, reference_psf, image_mask='', reference_mask='', n_stamps=4):
    """Calculate flux required to generate a transient"""

    science = ImageClass(image, image_psf, mask_filename=image_mask, n_stamps=n_stamps)
    reference = ImageClass(reference, reference_psf, mask_filename=reference_mask, n_stamps=n_stamps)
    flux_ratio = solve_iteratively(science, reference)
    return flux_ratio


def implant_transient(flux, image_filename, psf_filename):
    """Make an image containing a fake transient; contains code from ryan"""

    image = fits.getdata(image_filename)
    image_header = fits.getheader(image_filename)
    sexa_ra, sexa_dec = image_header['RA'], image_header['DEC']
    ra, dec = sexa2deg(sexa_ra, sexa_dec)
    x, y = degs2coords(ra, dec, image_header)
    image_coords = [round(float(coord)) for coord in (x, y)]
    psf = fits.getdata(psf_filename)
    psf /= np.sum(psf)
    fake_transient = make_transient(flux, psf, image.shape, image_coords)
    return image + fake_transient


def make_fake_image(image, image_psf, reference, reference_psf, flux):
    """Make an image with a fake transient of given flux as seen in the reference image"""

    output_list = ''

    hdu = fits.open(image)
    ra, dec = hdu[0].header['ra'], hdu[0].header['dec']
    print('Inserting transient with {0} counts at {1}, {2}.'.format(flux, ra, dec))
    fake_image = implant_transient(flux, image, image_psf)

    output_file = image.replace('.fits', '.{0}.fake.fits'.format(int(flux)))
    image_list = [output_file, image_psf, reference, reference_psf]

    hdu = fits.open(image)

    hdu[0].data = fake_image
    hdu.writeto(output_file, overwrite=True, output_verify='warn')

    output_list += '{0} {1} {2} {3}\n'.format(*image_list)


    file = open('images.txt', 'a')
    file.write(output_list)
    file.close()


def make_transient(flux, psf, shape, coords):
    """Plant a make a fake transient"""

    transient = flux * psf
    fake_image = np.zeros(shape)
    s = psf.shape
    ystart = int(coords[0] - s[0] / 2)
    xstart = int(coords[1] - s[1] / 2)
    fake_image[ystart: ystart + s[0], xstart: xstart + s[1]] = transient
    return fake_image


def sexa2deg(ra, dec):
    """Convert sexagesimal to degree; taken from ryan's code"""

    ra = coordinates.Angle(ra, units.hour).degree
    dec = coordinates.Angle(dec, units.degree).degree
    return ra, dec
