#!/usr/bin/env python
#
# :History:
# 
# 10 Mar 2015: Created.
# 08 Sep 2015: MiriImageModel changed to MiriMeasuredModel
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
# 15 Jun 2017: Added option to average multiple frames/groups
#              together. Use the MiriDarkReferenceModel methods.
# 26 Jun 2017: Added ability to extract the groups to a smaller 3-D dark.
#
# @author: Steven Beard (UKATC)
#
"""

Script `cdp_reduce_dark` combines one or more groups from a DARK data product.
DARK data products contain dark measurements as a function of both group
and integration. The files can be very large and reading them very memory
intensive. If the variation between groups is not important, this script
can produce a much smaller summary DARK.

The "averaging" option averages the groups and generates output data
in DN units.

The "slope" option calculates a mean slope and generates output data
in DN/group units.

The following command arguments are defined by position:

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with
        "_<start group>_<end group>" or "_allgroups" appended.

The command also takes the following options:

    --groupstart
        The first group/frame number to be extracted/averaged from the dark.
        Specifying only groupstart will extract only that group.
        Specifying neither groupstart nor groupend will average all groups.
    --groupend
        The last group/frame number to be extracted/averaged from the dark.
    --average or -a
        Average the data. This option overrides --slope if both are specified.
    --slope or -l
        Determine a slope of the data. If neither --average nor --slope are
        specified, the data will be extracted unchanged to a smaller 3-D
        DARK product.
    --simplefits or -s
        Use simple FITS I/O (in an attempt to save memory).
    --normalize
        normalize the resultant data so the average is 1.0.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""

import optparse
import sys, time

import numpy as np

# Import all the data models that might be contained in the file.
from miri.datamodels.cdp import MiriDarkReferenceModel
#from miri.datamodels.miri_measured_model import MiriMeasuredModel

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Extracts a frame from a MIRI dark data product."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--groupstart", dest="groupstart", type="int",
                     default=None, help="First frame number to extract"
                     )
    parser.add_option("", "--groupend", dest="groupend", type="int",
                     default=None, help="Last frame number to extract"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the copy of the file if it already exists"
                     )
    parser.add_option("-a", "--average", dest="average", action="store_true",
                      help="Average the data?"
                     )
    parser.add_option("-l", "--slope", dest="slope", action="store_true",
                      help="Determine a slope from the data?"
                     )
    parser.add_option("-s", "--simplefits", dest="simplefits", action="store_true",
                      help="Run simplified FITS code to save memory?"
                     )
    parser.add_option("-n", "--normalize", dest="normalize", action="store_true",
                      help="normalize the data so the average is 1.0?"
                     )

    (options, args) = parser.parse_args()

    groupstart = options.groupstart
    groupend = options.groupend
    overwrite = options.overwrite
    average = options.average
    slope = options.slope
    simplefits = options.simplefits
    normalize = options.normalize
    try:
        inputfile = args[0]
        if len(args) > 1:
            outputfile = args[1]
        else:
            if groupstart is not None and groupend is not None:
                if groupstart == groupend:
                    outputfile = inputfile + "_%d.fits" % groupstart
                else:
                    outputfile = inputfile + "_%d_%d.fits" % (groupstart, groupend)
            elif groupstart is not None:
                outputfile = inputfile + "_%d.fits" % groupstart
            else:
                outputfile = inputfile + "_allgroups.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    strg = ""
    if groupstart is not None and groupend is not None:
        if groupstart == groupend:
            strg += "Extracting group %d" % groupstart
        else:
            if average:
                strg += "Averaging groups %d-%d" % (groupstart, groupend)
            elif slope:
                strg += "Determining slope for groups %d-%d" % (groupstart, groupend)
            else:
                strg += "Extracting groups %d-%d" % (groupstart, groupend)
    elif groupstart is not None:
        strg += "Extracting group %d" % groupstart
        groupend = groupstart
    else:
        if average:
            strg += "Averaging all groups"
        elif slope:
            strg += "Determining slope for all groups"
        else:
            strg += "No group limits. Copying the whole data model."

    # Open the MIRI dark data model, extract the required frame
    # and save it to a MIRI image model.
    # THIS CODE IS SENSITIVE TO MEMORY ERROR. FULL FRAME DARKS MAY FAIL!
    if simplefits:
        print("%s to %s (using pyfits)" % (strg, outputfile))
        import pyfits
        # Simplified code which reads the FITS file directly
        hdulist = pyfits.open(inputfile)
  
        # Read the primary header
        fitsheader = hdulist[0].header
        # Extract the required frame from the SCI and ERR data arrays
        if groupstart is None:
            groupstart = 0
        if groupend is None:
            groupend = hdulist[1].data.shape[1]
        ngroups = groupend - groupstart + 1
        dataframe = hdulist[1].data[:,groupstart:groupend,:,:]
        dataframe = np.sum(dataframe, 1)
        if normalize:
            posdata = np.where(dataframe > 0.0)
            factor = np.mean(dataframe[posdata])
            print("Normalizing by a factor of %f" % factor)
            if factor > 0.0:
                dataframe = dataframe / factor  
        errframe = hdulist[2].data[:,groupstart:groupend,:,:]
        errsq = errframe * errframe
        errsq = np.sum(errsq, 1)
        errframe = np.sqrt( errsq ) / float(ngroups)
        if normalize:
            if factor > 0.0:
                errframe = errframe / factor
        dqframe = hdulist[3].data
        hdulist.close()
        
        newmodel = MiriDarkReferenceModel(data=dataframe, err=errframe, dq=dqframe)
        newmodel.save( outputfile, overwrite=overwrite )
        del newmodel
        del hdulist
       
    else:
        print("%s to %s (using data model)" % (strg, outputfile))
        with MiriDarkReferenceModel( init=inputfile ) as datamodel:
            if groupstart is None and groupend is None:
                if average:
                    newmodel = datamodel.average_groups(startgroup=None,
                                                        endgroup=None,
                                                        normalize=normalize)
                elif slope:
                    newmodel = datamodel.slope_data(startgroup=None,
                                                    endgroup=None)
                else:
                    newmodel = datamodel
                    
            elif groupstart == groupend:
                newmodel = datamodel.extract_group(group=groupstart,
                                                   normalize=normalize)
            else:
                if average:
                    newmodel = datamodel.average_groups(startgroup=groupstart,
                                                        endgroup=groupend,
                                                        normalize=normalize)
                elif slope:
                    newmodel = datamodel.slope_data(startgroup=groupstart,
                                                    endgroup=groupend)
                else:
                    newmodel = MiriDarkReferenceModel( \
                                    data=datamodel.data[:,groupstart:groupend,:,:],
                                    err=datamodel.err[:,groupstart:groupend,:,:],
                                    dq=datamodel.dq )
                    
            newmodel.save( outputfile, overwrite=overwrite )
            del newmodel
            del datamodel
