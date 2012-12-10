#-----------------------------------------------------------------------------
# Copyright (c) 2008  Raymond L. Buvel
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#-----------------------------------------------------------------------------

'''Utility functions for reading ephemeris files.

The following functions are available:

readHeader - open the ephemeris file, read the header, and return a tuple
    containing the header information.

julian - convert a calendar date to Julian date.
date - convert a Julian date to calendar date.
'''
__all__ = () # Don't export anything for "from ephemUtils import *"
import os,struct,numpy

bytesPerDouble = struct.calcsize('d')

#-----------------------------------------------------------------------------
# Assume the ephemeris data files are stored in the same directory with this
# module.  If you need a different directory structure, modify the way that the
# _ephemerisInfo dictionary is initialized.

_pathToData = os.path.dirname(os.path.abspath(__file__))

_ephemerisInfo = {}
for num, size in [(200,826), (403,1018), (405,1018), (406,728)]:
    _ephemerisInfo[num] = os.path.join(_pathToData, 'DE%03d.bin' % num), size

del _pathToData, num, size

#-----------------------------------------------------------------------------
def readHeader(ephemerisNumber):
    '''Open the ephemeris and read the header records.

    ephemerisNumber - numeric designation of the ephemeris to open.

    Returns a tuple containing the following:

    title - string of lines containing the ephemeris title information.
    timeData - tuple of (startTime, endTime, timeStep) units of Julian days
    dataStruct - structure parameters for the data records
    const - class instance containing the constants from the header.
    numRecords - number of data records in the ephemeris file.
    efile - open file object for the ephemeris.
    arrayBytes - number of bytes in each data record.
    '''
    fileName, arraySize = _ephemerisInfo[int(ephemerisNumber)]
    efile = file(fileName, 'rb')
    arrayBytes, numRecords = _getSizes(efile, arraySize)

    # Define the format for the structures in the header.  Specify Native
    # format so that no alignment is performed.
    fmt = ('=3d' # startTime, endTime, timeStep; units are Julian days
           'i'   # Number of constants
           '2d'  # AU, EMRAT - these values are not used (contained in const)
           '36i' # First 12 sets of data record structure
           'i'   # DENUM - this value is not used (contained in const)
           '3i') # Set 13 of data record structure
    size = struct.calcsize(fmt)

    partsFmt = ''.join([
        '%ds' % (84*3 + 6*400),         # Title and constant names
        '%ds' % size,                   # Data structure defined above
        '%dx' % (arrayBytes - (84*3 + 6*400) - size), # Skip pad
        '%ds' % (bytesPerDouble*400),  # Constant values
        '%dx' % (arrayBytes - bytesPerDouble*400), # Skip pad
    ])

    # Read the header records and unpack the major parts
    efile.seek(0)
    parts = struct.unpack(partsFmt, efile.read(2*arrayBytes))

    # Unpack the structures in the first record
    data = struct.unpack(fmt, parts[1])
    timeData = data[:3]
    numConst = data[3]

    # Validate the time data against the number of records in the file.
    numRec = (timeData[1] - timeData[0])/timeData[2]
    if abs(numRecords - numRec) > 1e-6:
        print 'numRecords = %s, numRec = %s' % (numRecords, numRec)
        raise TypeError('Ephemeris is invalid')

    # Pack the record structure values into a two dimensional array.
    dataStruct = numpy.array(data[6:6+36] + data[-3:])
    dataStruct = numpy.reshape(dataStruct, (13,3))

    # Compensate for the Fortran indexing
    for i in xrange(len(dataStruct)):
        dataStruct[i,0] -= 1

    # Unpack the title into a string of lines terminated by newlines.
    lst = struct.unpack('84s'*3, parts[0][:84*3])
    title = ''.join([x.strip()+'\n' for x in lst])

    # Unpack the constant names and values.  Make the information available as
    # instance attributes for user convenience.
    const = _ConstData()
    lst = struct.unpack('6s'*numConst, parts[0][84*3:84*3+6*numConst])
    names = tuple([x.strip() for x in lst])
    values = numpy.fromstring(parts[2], dtype=numpy.float64)
    for i in xrange(numConst):
        # Bypass the protection mechanism to set the values
        const.__dict__[names[i]] = values[i]
    # Save the names so they can be used to print the const values.
    const.__dict__['_constNames'] = names

    # Verify the unused header values since the JPL Fortran program uses them.
    # This ensures that the translator programs have operated correctly.
    if data[4] != const.AU or data[5] != const.EMRAT or data[-4] != const.DENUM:
        raise ValueError('Header is invalid')

    return title, timeData, dataStruct, const, numRecords, efile, arrayBytes

#-----------------------------------------------------------------------------
def julian(*args):
    '''Convert a calendar date to Julian date.

    This function uses the calendar date format used in the JPL test data.

    The date can be input in one of two formats:

    string - a string having the format year.month.day.  Note that abbreviated
        years are not allowed.  Years continue through zero to negative.

    year, month, day - three integers
    '''
    if len(args) == 3:
        year,month,day = args
    elif len(args) == 1:
        year,month,day = [int(x) for x in args[0].split('.')]
    else:
        raise ValueError('Invalid format')

    if year > 1582:
        greg = True
    elif year == 1582 and month > 10:
        greg = True
    elif year == 1582 and month == 10 and day >= 15:
        greg = True
    else:
        greg = False

    a = (14-month)//12
    y = year + 4800 - a
    m = month + 12*a - 3

    if greg:
        JD = day + (153*m+2)//5 + 365*y + y//4 - y//100 + y//400 - 32045.5
    else:
        JD = day + (153*m+2)//5 + 365*y + y//4 - 32083.5

    return JD

#-----------------------------------------------------------------------------
def date(JD):
    '''Convert a Julian date to calendar date

    returns: (year, month, day)
    '''
    if JD <= 2299159.5:
        # Julian calendar
        b = 0
        c = int(JD + 32082.5)
    else:
        # Gregorian calendar
        a = int(JD + 32044.5)
        b = (4*a+3)//146097
        c = a - (146097*b)//4

    d = (4*c+3)//1461
    e = c - (1461*d)//4
    m = (5*e+2)//153

    m10 = m//10

    day = e - (153*m+2)//5 + 1
    month = m + 3 - 12*m10
    year = 100*b + d - 4800 + m10

    return year, month, day

#-----------------------------------------------------------------------------
class _ConstData:
    '''Container class used to hold the constant data read from the header.
    The values are protected from writing using normal asignment.
    '''
    def __setattr__(self, name, value):
        raise AttributeError('Attributes are read-only')

    def __delattr__(self, name, value):
        raise AttributeError('Attributes are read-only')

    def __str__(self):
        names = self._constNames
        lst = ['%s = %s\n' % (name, getattr(self,name)) for name in names]
        return ''.join(lst)

#-----------------------------------------------------------------------------
def _getSizes(efile, arraySize):
    arrayBytes = arraySize*bytesPerDouble

    # Find the number of bytes in the file and verify that the file contains
    # an integral number of records.
    efile.seek(0,2)
    bytes = efile.tell()
    if (bytes % arrayBytes) != 0:
        raise TypeError('Ephemeris size is invalid')

    # Compute the number of data records.
    numRecords = bytes/arrayBytes - 2  # Don't include the header records    
    return arrayBytes, numRecords

