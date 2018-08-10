#!/usr/bin/env python
#
# Script `make_sca_file` reads detector illumination data either from an ASCII
# file or from a standard image format file (JPEG, TIFF, GIF, PNG, PPM, etc...)
# and writes a FITS file in the format readable by the MIRI SCA simulator.
#
# :History:
# 
# 05 Jul 2010: Created
# 16 Jul 2010: Try importing the IlluminationMap module in two different ways.
#             Corrected a bool/int incompatibility problem with the
#             verbose parameter. Added --silent option.
# 04 Aug 2010: Optional arguments converted to parser options.
#              Wavelength parameters added.
# 06 Aug 2010: Documentation formatting problems corrected.
# 03 Sep 2010: illumination_maps renamed to data_maps
# 27 Sep 2010: Python environment for windows verified.
# 12 Oct 2010: Unwanted option strings removed.
# 10 Jan 2011: Added optional scaid and subarray parameters.
# 13 Nov 2012: Major restructuring of the package folder.
#              Import statements updated.
# 18 Jun 2013: Modified to match the new parameter names used
#              by data_maps. Fixed header->metadata typo.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#              Corrected detector names.
# 14 Aug 2015: Tidied up line breaks.
# 08 Sep 2015: Removed dependency on legacy data maps.
#              check_wavelength and load_illumination_map functions
#              moved here. Made compatible with Python 3.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
# 22 Jan 2018: Protect against incompatibility problems from legacy Image
#              package.
#
# @author: Steven Beard (UKATC)

"""

The make_sca_file script reads detector illumination data from a file in
the specified format (IMAGE or ASCII) and writes a FITS file in the format
readable by the MIRI SCA simulator. This utility is useful for creating
test data for the SCA simulator. The compulsory command arguments are:

    inputfile
        The path+name of the file to be read.
    outputfile
        The path+name of the SCA format FITS file to be created.

The following optional parameters may be provided by keyword:

    --filetype
        The type of input file provided.

        * If 'ASCII', intensity data are read from an ASCII file and
          the wavelength and direction data are ignored.
        * If 'IMAGE', intensity data are read from a standard image
          format file (such as JPEG, GIF or PPM) whose contents are
          converted into a composite image. The wavelength data is
          generated from the minimum and maximum wavelengths provided
          and the direction data are ignored.
        * A file type of 'FITS' is also accepted, but converting FITS
          format data into FITS format data is not useful.
          
        If a file type is not specified, the default is 'IMAGE'.
    --scale
        An optional scale factor to apply to the intensity data
        read from the file. This is more useful for 'IMAGE' data in
        cases where the contents of the file do not give a photon
        flux.
    --wavmin
        For data obtained from an IMAGE file, the minimum wavelength
        in microns associated with the data.
        A wavelength extension will only be created if wavmin and
        wavmax are both specified.
    --wavmax
        For data obtained from an IMAGE file, the maximum wavelength
        in microns associated with the data.
        A wavelength extension will only be created if wavmin and
        wavmax are both specified.
    --scaid
        The identify of the sensor chip assembly module with which
        the data are associated (e.g. 'MIRIMAGE', 'MIRIFULONG' or 'MIRIFUSHORT').
        Any digits provided will be written to an integer SCA_ID
        keyword in the FITS header of the output file. If not
        provided, no keyword will be written. 
    --subarray
        The subarray mode for the data. Any string provided will be
        written to a SUBMODE keyword in the FITS header of the
        output file. If not provided, no keyword will be written.

The command also takes the following options:

    --silent or -s:
        Generate no output.
    --verbose or -v:
        Generate more output and some data plots.
    --debug or -d:
        Generate maximum output.
    --overwrite or -o:
        Overwrite any existing SCA format FITS file.

"""

import optparse
import os, sys, time
import re
import numpy as np
import astropy.io.fits as pyfits

try:
    import Image
    _PIL_AVAILABLE = True
except Exception:
    # NOTE: Incompatibility problems may cause a number of different exceptions
    # while attempting to import the Image package. Catch all of them.
    Image = None
    _PIL_AVAILABLE = False

from miri.datamodels.sim import MiriIlluminationModel

def _np_from_jpeg( filename ):
    """
    
    Open a JPEG file and extract a numpy array from its contents.
    
    Requires the Python Imaging Library, PIL.
    
    """
    if not _PIL_AVAILABLE:
       raise NotImplementedError("Cannot open JPEG file without Python Imaging Library")
    # Open the image file, convert the contents to a composite image
    # and transpose it to numpy orientation.
    jpim1 = Image.open(filename)
    jpim2 = jpim1.convert('L')
    jpim3 = jpim2.transpose(Image.ROTATE_180)
    jpim4 = jpim3.transpose(Image.FLIP_LEFT_RIGHT)
    # Convert to a numpy array of the correct shape, reducing the
    # amplitude by the given scale factor.
    data = np.array(list(jpim4.getdata()))
    data.shape = (jpim4.size[1], jpim4.size[0])

    # Explicitly delete the image copies rather than leaving it
    # to the Python garbage collector, as they could be large.
    del jpim1, jpim2, jpim3, jpim4
    return data

def check_wavelength(data, shape, slices):
    """
    
    Verifies that a wavelength array is compatible with the shape
    of the illumination array and number of slices.
    
    :Parameters:
    
    data: array_like
        The wavelength data array to be checked.
    shape: tuple of 2 ints
        The illumination data shape.
    slices: int
        The number of slices of illumination data included
        in the intensity array.
    
    :Returns:
    
    valid: bool
        True if the array is compatible and False if not compatible.
        
    """
    # The wavelength array can be null.
    if data is None: return True
    
    npdata = np.asarray(data)
    
    # If the intensity array is 2 dimensional
    if slices == 0:
        # The wavelength array is compatible if it is the same shape
        if npdata.shape == shape:
            return True
    else:
        # For 3 dimensional data the wavelength array must match the
        # first 2 dimensions and its 3rd dimension must either be absent,
        # be 1 or be the same as the number of slices.
        if npdata.shape[-2:] == shape:
            if len(npdata.shape) < 3:
                return True
            if (npdata.shape[0] == 1) or (npdata.shape[0] == slices):
                return True

    # The wavelength array is also valid if 3 dimensional with the rows and
    # columns exactly 1 and the number of slices matching the intensity array.
    if len(npdata.shape) == 3:
        if (npdata.shape[1] == 1) and (npdata.shape[2] == 1) and \
            (npdata.shape[0] == slices):
            return True
    # Otherwise the wavelength array is not compatible.
    return False

def load_illumination_map(filename, ftype='FITS', metadata=None, scale=1.0,
                          add_wavelength=True, wavmin=1.0, wavmax=30.0):
    """
    
    Create an IlluminationMap object from a file of detector illumination data.
                 
    :Parameters:
    
    filename: string or tuple of strings
        The name(s) of the file(s) to be opened.
    ftype: string, optional, default='FITS'
        The type of file from which to read the data.
        
        * If 'STSCI', the illumination map is read from a standard MIRI
          illumination data file, based on the STScI data model.
        * If 'FITS', intensity and wavelength data are obtained
          from a single FITS file containing INTENSITY, WAVELENGTH and
          DIRECTION extensions (with the latter two extensions being
          optional).
        * If 'ASCII', intensity and wavelength data are read from
          a list of 1, 2 or 3 ASCII files (with the latter two files
          being optional). Intensity values are multiplied by the given
          scale factor.
          ASCII files describing an illumination map define a value for every
          pixel. They should contain nrows lines, each of which contains
          ncolumns of blank-separated numbers; i.e. the default format
          expected by the numpy.loadtxt() function.
        * If 'IMAGE' (or 'JPEG'), intensity data are read from a standard
          image format file (such as JPEG, GIF or PPM) whose contents are
          converted into a composite image. A simple wavelength image may be
          added varying between a specified minimum and maximum wavelength,
          or wavelength data may be ignored. Intensity values are multiplied
          by the given scale factor.
          
    metadata: dictionary-like object, optional
        An object containing metadata to be associated with the data.
        It could be a plain Python dictionary, a Metadata object or
        a pyFits Header object, as long as it supports keyword operations.
        (Only used if the input file is ASCII or IMAGE.)
    scale: float, optional, default=1.0
        An optional scale factor to apply to the intensity data. This is
        more useful for ftype='IMAGE' data in cases where the data do not
        contain a photon flux.
    add_wavelength: boolean, optional, default=True
        Add a test wavelength map to the data. This parameter is valid
        only for ftype='IMAGE' data.
    wavmin: float, optional, default=1.0
        Minimum wavelength associated with data in microns. This parameter
        is valid only for ftype='IMAGE' data.
    wavmax: float, optional, default=30.0
        Maximum wavelength associated with data in microns. This parameter
        is valid only for ftype='IMAGE' data.
       
    :Raises:
    
    ValueError
        Raised if any of the parameters are out of range.
    IOError
        Raised if there is an error opening, reading or interpreting
        the data from the input file.
    ImportError
        Raised if the Python Imaging Library (PIL) is not available.
        
    :Returns:
    
    A new MiriIlluminationModel object.
        
    """
    # Attempt to open and then read the given file
    if ftype == 'STSCI' or ftype == 'STScI':
    
        # Read a MIRI data model in standard STScI format.
        illumination_map = MiriIlluminationModel( filename )
        illumination_map.apply_scale(scale)
        return illumination_map

    elif ftype == 'FITS':
    
        # Read a FITS file in the agreed detector illumination format.
        # The file must contain a primary FITS header and an INTENSITY
        # extension. It may also contain WAVELENGTH and DIRECTION
        # extensions.
        try:
            hdulist = pyfits.open(filename)
 
            # Read the contents of the file       
            header = hdulist[0].header
            description = 'Metadata extracted from %s' % filename
            metadata = Metadata(description)
            metadata.from_fits_header(header)
            
            intensity = scale * hdulist['INTENSITY'].data
            intensity_header = hdulist['INTENSITY'].header
            intensity_metadata = Metadata('Intensity metadata')
            intensity_metadata.from_fits_header(intensity_header)
        except Exception as e:
            # If the file could not be read re-raise the exception
            # with a more meaningful error message.
            strg = "%s: Could not read illumination data file %s.\n   %s" % \
                (e.__class__.__name__, filename, e)
            raise IOError(strg)

        # The wavelength data is optional
        try:
            wavelength = hdulist['WAVELENGTH'].data
            wavelength_header = hdulist['WAVELENGTH'].header
            wavelength_metadata = Metadata('Wavelength metadata')
            wavelength_metadata.from_fits_header(wavelength_header)
        except (KeyError, AttributeError):
            wavelength = None
            wavelength_metadata = None

        hdulist.close()

    elif ftype == 'ASCII':
        if isinstance(filename, str):
            # A single string is provided - there is just an intensity file
            intensity = scale * np.loadtxt(filename)
            intensity_metadata = None
            wavelength = None
            wavelength_metadata = None
        else:
            # A list of strings has been provided. Attempt to open each
            # file, ignoring the wavelength file if not specified or it
            # doesn't exist.
            intensity = scale * np.loadtxt(filename[0])
            intensity_metadata = None
            try:
                wavelength = np.loadtxt(filename[1])
            except Exception:
                wavelength = None
            wavelength_metadata = None

    elif ftype == 'IMAGE' or ftype == 'JPEG':
        if _PIL_AVAILABLE:
            # Convert to a numpy array of the correct shape, reducing the
            # amplitude by the given scale factor.
            data = _np_from_jpeg(filename)
            intensity = scale * data.astype(np.float32)
            intensity_metadata = None

            # Add some test wavelength data if required. The wavelength will
            # increase linearly from bottom to top over the range specified.
            if add_wavelength:
                wavelength = np.empty_like(intensity)
                wav = np.linspace(wavmin, wavmax, wavelength.shape[1])
                wavelength[:,:] = wav
                wavelength = np.transpose(wavelength)
            else:
                wavelength = None
            wavelength_metadata = None

        else:
            strg = "Sorry, file format %s can't be processed " \
                "because the Python Imaging Library is not available." % ftype
            raise ImportError(strg)
           
    else:
        # Other file formats are not yet supported
        strg = "Sorry, file format %s is not supported." % ftype
        raise ValueError(strg)

    # Bail out if the intensity data has not been read successfully.
    if intensity is None:
        strg = "No intensity data could be found within file.\n   %s" % \
            filename
        raise ValueError(strg)

    # Ensure that the intensity and wavelength arrays are
    # compatible.
    illumination_shape = intensity.shape[-2:]
    if len(intensity.shape) > 2:
        slices = intensity.shape[0]
    else:
        slices = 0

    # If the wavelength array is 1-D and the intensity array is 3-D
    # the wavelength array needs to be reshaped into a 3-D array.
    if wavelength is not None:
        if len(wavelength.shape) == 1 and len(intensity.shape) > 2:
            sz1 = wavelength.shape[0]
            wavelength.shape = [sz1, 1, 1]

    if check_wavelength(wavelength, illumination_shape, slices) == False:
        print( "WARNING: Wavelength array has unexpected size - ignoring it." )    
        wavelength = None
 
    # Finally, create and return an IlluminationMap object
    illum = MiriIlluminationModel(intensity=intensity, wavelength=wavelength)  
    return illum

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile"
    usage += "\n\t[--filetype] [--scale] [--wavmin] [--wavmax]"
    parser = optparse.OptionParser(usage)
    
    # Optional arguments (long option strings only).
    parser.add_option("", "--filetype", dest="filetype", type="string",
                     default='IMAGE', help="File type"
                     )
    parser.add_option("", "--scale", dest="scale", type="float",
                     default=1.0, help="Intensity scale factor"
                     )
    parser.add_option("", "--wavmin", dest="wavmin", type="float",
                     default=None, help="Minimum wavelength (microns)"
                     )
    parser.add_option("", "--wavmax", dest="wavmax", type="float",
                     default=None, help="Maximum wavelength (microns)"
                     )
    parser.add_option("", "--scaid", dest="scaid", type="string",
                     default='', help="Sensor Chip Assembly ID"
                     )
    parser.add_option("", "--subarray", dest="subarray", type="string",
                     default='', help="Detector subarray mode for output"
                     )

    # Boolean options (short and long option strings).
    parser.add_option("-d", "--debug", dest="debug", action="store_true",
                      help="Debugging mode"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-s", "--silent", dest="silent", action="store_true",
                      help="Silent mode"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite existing file"
                     )

    (options, args) = parser.parse_args()
    
    # Compulsory arguments
    try:
        inputfile = args[0]
        outputfile = args[1]
    except IndexError:
        print( help_text )
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)
        
    # Optional arguments.
    filetype = options.filetype
    scale = float(options.scale)
    wavmin = options.wavmin
    wavmax = options.wavmax
    scaid = options.scaid
    subarray_str = options.subarray
    
    if wavmin is not None and wavmax is not None:
        add_wavelength = True
        wavmin = float(wavmin)
        wavmax = float(wavmax)
    else:
        add_wavelength = False

    # Boolean flags.
    overwrite = options.overwrite
    verb = options.verb
    debug = options.debug
    silent = options.silent
    
    # Set the verbosity level according to the --verbose and --silent
    # options. (Note that --debug wins over --verbose and --silent wins
    # over all the other options if they are provided together.)
    verbose = 1
    if verb: verbose = 2
    if debug: verbose = 4
    if silent: verbose = 0

    # Describe the file using HISTORY keywords.
    # Note: This header info is only added when data are read from an
    # ASCII or JPEG file. If read from a FITS file, the header contained
    # in the FITS file is used.
    keyw01 = pyfits.Card('HISTORY', 'Written by the make_sca_file script.',
                         '')
    str02 = "Test data %s file %s" % (filetype, inputfile)
    kws02 = "%.64s" % str02
    keyw02 = pyfits.Card('HISTORY', kws02, '')
    str03 = "Original intensity data multiplied by %f" % scale
    kws03 = "%.64s" % str03
    keyw03 = pyfits.Card('HISTORY', kws03, '')
    pheader = pyfits.Header(cards=[keyw01, keyw02, keyw03])
    
    # If specified, add the SCA_ID and/or SUBMODE keywords.
    if scaid:
        # Extract only digits from the scaid to make an integer.
        digits = re.compile(r'([0-9]+)')
        scaid_int = digits.findall(scaid)
        if len(scaid_int) > 0 and scaid_int[0]:
            pheader.update("SCA_ID", int(scaid_int[0]),
                           "Sensor Chip Assembly ID")
    if subarray_str:
        pheader.update("SUBMODE", subarray_str, "Detector subarray mode")

    if verbose > 0:
        print( "Converting the %s file %s into FITS file %s." % \
            (filetype, inputfile, outputfile) )
        if abs(scale-1.0) > 0.0001:
            print( "Data will be scaled by %f." % scale )
        if (filetype == 'IMAGE' or filetype == 'JPEG') and add_wavelength:
            print( "A wavelength extension will be added covering the range " \
                "%.2f to %.2f microns." % (wavmin,wavmax) )

    # Create a new illumination map object from the data,
    map_object = load_illumination_map(inputfile, ftype=filetype,
                                       metadata=pheader,
                                       scale=scale,
                                       add_wavelength=add_wavelength,
                                       wavmin=wavmin, wavmax=wavmax)
    if verbose > 1:
        print( map_object )

    # Save the illumination map to a FITS file.
    map_object.save(outputfile, overwrite=overwrite)
    
    if verbose > 0:
        print( "New illumination map saved." )
