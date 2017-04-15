from astropy.io import fits
from astropy import wcs
from astropy import coordinates
from astropy import units
import numpy as np
from PyZOGY.subtract import ImageClass
from PyZOGY.util import solve_iteratively
import os


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


def make_transient(flux, psf, shape, coords):
    """Plant a make a fake transient"""

    transient = flux * psf
    fake_image = np.zeros(shape)
    s = psf.shape
    ystart = int(coords[0] - s[0] / 2)
    xstart = int(coords[1] - s[1] / 2)
    fake_image[ystart:ystart+s[0],xstart:xstart+s[1]] = transient
    return fake_image


def find_flux_ratio(image, image_psf, reference, reference_psf, image_mask='', reference_mask='', n_stamps=4):
    """Calculate flux required to generate a transient"""

    science = ImageClass(image, image_psf, mask_filename=image_mask, n_stamps=n_stamps)
    reference = ImageClass(reference, reference_psf, mask_filename=reference_mask, n_stamps=n_stamps)
    flux_ratio = solve_iteratively(science, reference)
    return flux_ratio


def make_fake_image(image, image_psf, reference, reference_psf, magnitude,
                    flux_ratio=np.inf, image_mask='None', reference_mask='None'):
    """Make an image with a fake transient of given magnitude as seen in the reference image"""

    if flux_ratio == np.inf:
        flux_ratio = find_flux_ratio(image, image_psf, reference, reference_psf,
                                     image_mask=image_mask, reference_mask=reference_mask)
    output_list = ''
    path = os.getcwd() + '/'

    for mag in magnitude:
        flux = 10 ** (-mag / 2.5)
        print(flux)
        flux *= flux_ratio
        print(flux, flux_ratio)
        fake_image = implant_transient(flux, image, image_psf)

        output_file = 'mag{}.fits'.format(str(mag))
        image_list = [output_file, image_psf, image_mask, reference, reference_psf, reference_mask]
        image_list_path = []
        for image in image_list:
            if image != 'None':
                image_list_path.append(path + image)
            else:
                image_list_path.append(image)
        fits.writeto(output_file, fake_image, overwrite=True, output_verify='warn')
        output_list += '{0},{1},{2},{3},{4},{5}\n'.format(*image_list_path)

    file = open('images.txt', 'w')
    file.write(output_list)
    file.close()


def degs2coords(ra, dec, header):
    """Converts RA and DEC into coordinates on image; taken from ryan's code"""

    w = wcs.WCS(header)
    return w.all_world2pix(ra,dec,0)


def sexa2deg(ra,dec):
    """Convert sexagesimal to degree; taken from ryan's code"""

    ra = coordinates.Angle(ra, units.hour).degree
    dec = coordinates.Angle(dec, units.degree).degree
    return ra, dec