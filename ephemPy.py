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

'''Module for interpolating JPL ephemeris tables.

To use this module, create an instance of the class Ephemeris and call the
appropriate methods to extract the required information.  You may also want to
sub-class Ephemeris if the results are required in a different form.  See the
test program testEphem.py for an example.
'''
__all__ = () # Don't export anything for "from ephemPy import *"

import ephemUtils
import numpy

#-----------------------------------------------------------------------------
class Ephemeris:
    '''Interpolate JPL ephemeris tables.

    An ephemeris object is created with the following call.

    Ephemeris(ephemerisNumber)
        ephemerisNumber - numeric designation of the ephemeris to open.

    The following attributes are extracted from the header records of the
    selected ephemeris.

    title - a string containing the three lines in the title section.
    eStartTime - starting Julian day of the ephemeris
    eEndTime - ending Julian day of the ephemeris
    eTimeStep - time interval covered by each record
    constants - class instance containing the constants found in the header
        records.  For example, to get the number of kilometers in an
        Astronomical Unit, use constants.AU.
    dataStruct - array containing the structure parameters for a data
        record.  See the JPL documentation for details.
    numRecords - number of data records in the ephemeris file.

    Additional attributes:

    record - array containing the current data record
    rStartTime - starting Julian day of the current record
    rEndTime - ending Julian day of the current record
    refTime - reference time set by the user with the setRefTime method.  The
        time parameter to the other methods is relative to this value.
    arrayBytes - number of bytes in a data record.

    The following attributes are used to select the target.

    Note: the EARTH target is the Earth-Moon barycenter and the MOON is
    relative to the the geocenter.

    MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO,
    MOON, SUN

    The following are required by the test program but must not be used as
    targets for any of the methods.

    SS_BARY, EM_BARY, NUTATIONS, LIBRATIONS

    The following methods are available.

    position(t, target)
        Interpolate the position vector of the target.

    state(t, target)
        Interpolate the state vector of the target.

    nutations(t)
        Interpolate nutations

    librations(t)
        Interpolate librations

    setRefTime(t)
        Set the reference time to the start of the record containing t.

    getRecord(t)
        Get the record corresponding to the specified time.
    '''

    # The following constants can be used for the target parameter of the
    # methods requiring a target.
    MERCURY = 0
    VENUS = 1
    EARTH = 2        # Earth-Moon Barycenter
    MARS = 3
    JUPITER = 4
    SATURN = 5
    URANUS = 6
    NEPTUNE = 7
    PLUTO = 8
    MOON = 9         # Relative to geocenter
    SUN = 10

    # The following are required by the test program but must not be used as
    # targets for any of the methods.
    SS_BARY = 11
    EM_BARY = 12
    NUTATIONS = 13
    LIBRATIONS = 14

    def __init__(self, ephemerisNumber):
        # Read the header and assign the result to instance attributes.
        hdr = ephemUtils.readHeader(ephemerisNumber)
        self.title = hdr[0]
        self.eStartTime, self.eEndTime, self.eTimeStep = hdr[1]
        self.dataStruct = hdr[2]
        self.constants = hdr[3]
        self.numRecords = hdr[4]
        self._efile = hdr[5]
        self.arrayBytes = hdr[6]

        # Read the first record to set the record times.
        self._readRecord(0)

        # Initialize the reference time so that the time parameter represents
        # the Julian day values in the ephemeris file.
        self.refTime = 0.0


    def position(self, t, target):
        '''Interpolate the position vector of the target

        Returns an array containing the position measured in kilometers.  For
        all targets except the Moon, the position is relative to the Solar
        System barycenter.  For the Moon, the position is relative to the
        geocenter.

        t - time in Julian days at which the position is desired.
        target - object for which the position is desired [0,...,10].
        '''
        if not (0 <= target <= 10):
            raise ValueError('target out of range')
        Tc, dt, A = self._getParms(t, target)
        pos = numpy.zeros((3,), numpy.float64)
        for i in xrange(3):
            pos[i] = chebeval(A[i], Tc)
        return pos


    def state(self, t, target):
        '''Interpolate the state vector of the target

        Returns an array containing the state vector of the target.  The
        position is in the first three elements and is measured in kilometers.
        The velocity is in the last three elements and is measured in
        kilometers per Julian day.  For all targets except the Moon, the
        position is relative to the Solar System barycenter.  For the Moon, the
        position is relative to the geocenter.

        t - time in Julian days at which the state vector is desired.
        target - object for which the state vector is desired [0,...,10].
        '''
        if not (0 <= target <= 10):
            raise ValueError('target out of range')
        Tc, dt, A = self._getParms(t, target)
        dt2 = 2.0/dt
        PV = numpy.zeros((6,), numpy.float64)
        for i in xrange(3):
            PV[i] = chebeval(A[i], Tc)
            PV[i+3] = chebder(A[i], Tc)*dt2
        return PV


    def nutations(self, t):
        '''Interpolate nutations

        t - time in Julian days at which the nutations are desired.
        '''
        if self.dataStruct[11,1] < 2:
            raise TypeError('Ephemeris does not contain nutations')
        Tc, dt, A = self._getParms(t, 11)
        dt2 = 2.0/dt
        NU = numpy.zeros((4,), numpy.float64)
        for i in xrange(2):
            NU[i] = chebeval(A[i], Tc)
            NU[i+2] = chebder(A[i], Tc)*dt2
        return NU


    def librations(self, t):
        '''Interpolate librations

        t - time in Julian days at which the Librations are desired.
        '''
        if self.dataStruct[12,1] < 2:
            raise TypeError('Ephemeris does not contain librations')
        Tc, dt, A = self._getParms(t, 12)
        dt2 = 2.0/dt
        LI = numpy.zeros((6,), numpy.float64)
        for i in xrange(3):
            LI[i] = chebeval(A[i], Tc)
            LI[i+3] = chebder(A[i], Tc)*dt2
        return LI


    def setRefTime(self, t):
        '''Set the reference time to the start of the record containing t.

        If t == 0, the reference time is removed.

        Returns the difference from the value set as the reference time.

        t - time in Julian days
        '''
        if t == 0:
            self.refTime = 0.0
            return 0.0

        if t < self.eStartTime or t > self.eEndTime:
            raise ValueError('Time out of range')

        self._readRecord((t-self.eStartTime)/self.eTimeStep)
        self.refTime = self.rStartTime
        return t - self.refTime


    def getRecord(self, t):
        '''Get the record corresponding to the specified time.

        t - time in Julian days
        '''
        t = t + self.refTime
        if t >= self.rStartTime and t <= self.rEndTime:
            return self.record

        if t < self.eStartTime or t > self.eEndTime:
            raise ValueError('Time out of range')

        self._readRecord((t-self.eStartTime)/self.eTimeStep)

        if t < self.rStartTime or t > self.rEndTime:
            raise ValueError('Invalid record')
        return self.record


    def _readRecord(self, num):
        num = int(num)
        if num < 0 or num >= self.numRecords:
            raise IndexError('Record number out of range')

        self._efile.seek((num+2)*self.arrayBytes)
        data = self._efile.read(self.arrayBytes)
        self.record = numpy.fromstring(data, numpy.float64)

        self.rStartTime = self.record[0]
        self.rEndTime = self.record[1]


    def _getParms(self, t, target):
        rec = self.getRecord(t)
        timeStep = self.eTimeStep
        t0 = rec[0] - self.refTime

        # Get structure parameters for the specified target
        C,N,G = self.dataStruct[target]

        if G == 1:
            dt = timeStep
            Tc = 2.0*(t-t0)/dt - 1.0
        else:
            dt = timeStep/G  # Time step per granule
            i = int((t-t0)/dt)
            if i == G: i = i-1 # This can happen if the time is the endpoint

            if target == 11:
                # Nutations only have two entries
                C = C + i*2*N
            else:
                C = C + i*3*N
            Tc = 2.0*((t-t0)-i*dt)/dt - 1.0

        if target == 11:
            # Nutations only have two entries
            A = numpy.reshape(rec[C:C+2*N],(2,N))
        else:
            A = numpy.reshape(rec[C:C+3*N],(3,N))
        return Tc, dt, A

#-----------------------------------------------------------------------------
def chebeval(coef, x):
    '''Evaluate a Chebyshev polynomial

    coef - sequence of coefficients
    x - point where the polynomial is evaluated -1 <= x <= 1
    '''
    if x < -1.0 or x > 1.0:
        raise ValueError('parameter out of range')

    # Use Clenshaw's recurrence to evaluate the polynomial.  See Numerical
    # Recipes for a discussion.
    x2 = 2.0*x
    d = dd = 0.0
    k = len(coef)-1
    while k >= 1:
        dd, d = d, x2*d - dd + coef[k]
        k -= 1
    return x*d - dd + coef[0]

#-----------------------------------------------------------------------------
def chebder(coef, x):
    '''Evaluate the derivative of a Chebyshev polynomial

    coef - sequence of coefficients
    x - point where the derivative is evaluated -1 <= x <= 1
    '''
    if x < -1.0 or x > 1.0:
        raise ValueError('parameter out of range')

    # Use Clenshaw's recurrence to evaluate the polynomial.  See Numerical
    # Recipes for a discussion.
    x2 = 2.0*x
    d = dd = 0.0
    k = len(coef)-1
    while k >= 2:
        dd, d = d, x2*d - dd + k*coef[k]
        k -= 1
    return x2*d - dd + coef[1]

