'''
Pratheek Nagaraj
July 11, 2011
Astronomy Physics Homework 4 Question 4

Orbital Elements
'''

from visual import *
from math import *
from ephemPy import Ephemeris as Ephemeris_BC
from Function_Library import *

k = 0.01720209895
G = 6.673e-11
epsilon = 23.43799 #degrees

class Ephemeris(Ephemeris_BC):
    def __init__(self, *args, **kwargs):
        Ephemeris_BC.__init__(self, *args, **kwargs)
        self.AUFAC = 1.0/self.constants.AU
        self.EMFAC = 1.0/(1.0+self.constants.EMRAT)
    def position(self, t, target, center):
        pos = self._position(t, target)
        if center != self.SS_BARY:
            pos = pos - self._position(t, center)
        return pos    
    def _position(self, t, target):
        if target == self.SS_BARY:
            return numpy.zeros((3,), numpy.float64)
        if target == self.EM_BARY:
            return Ephemeris_BC.position(self, t, self.EARTH)*self.AUFAC
        pos = Ephemeris_BC.position(self, t, target)*self.AUFAC
        if target == self.EARTH:
            mpos = Ephemeris_BC.position(self, t, self.MOON)*self.AUFAC
            pos = pos - mpos*self.EMFAC
        elif target == self.MOON:
            epos = Ephemeris_BC.position(self, t, self.EARTH)*self.AUFAC
            pos = pos + epos - pos*self.EMFAC
        return pos

ephem = Ephemeris('405')

def isLeap( year ):

    if year % 400 == 0:
        return True
    elif year % 100 == 0:
        return False
    elif year % 4 == 0:
        return True
    else:
        return False

def GDInput():
    year = input("Please enter year (2000,2001,...): ")
    month = input("Please enter month (1,2,3,...,12): ")
    day = input("Please enter day (1,2,3,...): ")
    hour = input("Please enter hours (0,1,2,...23): ")
    minutes = input("Please enter minutes (0,1,...,59): ")
    seconds = input("Please enter seconds (0,1,...,59): ")

    GD( year, month, day, hours, minutes, seconds )

def GD( year, month, day, hours, minutes, seconds ):

    monthArray = [31,28,31,30,31,30,31,31,30,31,30,31]
    monthLeapArray = [31,29,31,30,31,30,31,31,30,31,30,31]

    leapCount = floor( (year - 2001)/4 + 1 )
    nonLeapCount = year - 2000 - leapCount

    days = 0
    days = days + leapCount*366 + nonLeapCount*365

    if isLeap(year):
        days = days + sum( [monthLeapArray[i] for i in range(0,month-1)] )
    else:
        days = days + sum( [monthArray[i] for i in range(0,month-1)] )

    days = days + day - 1
    hours = hours + days * 24
    JDN = 2451545.0

    JD = JDN + (hours-12)/24. + minutes/1440. + seconds/86400.
    return JD

def toDegrees( rad ):
    return rad*180/pi

def toRadians( deg ):
    return deg*pi/180

def NewtonMethod( e, M, tolerance = 1e-10 ):
    xo = M
    zero = 1
    
    while fabs(xo - e*sin(xo) - M) > tolerance:
        zero = ( ( xo - e*sin(xo) - M)/(1 - e*cos(xo) ) )
        x1 = xo - zero
        xo = x1

    return xo

def AUtoMeters( AU ):
    return 149598000000*AU

def MeterstoAU( meters ):
    return meters/149598000000

fileIn = open( "data.txt", "r" )

allList = []

for line in fileIn:
    allList.append( line.split() )

mu = 1
e = float(allList[0][0])
i = toRadians(float(allList[1][0]))
omega = toRadians(float(allList[2][0]))
w = toRadians(float(allList[3][0]))
a = float(allList[4][0])
Tp = float(allList[5][0])
TobsArray = []
for element in allList[6]:
    TobsArray.append(int(element))
Tobs = GD( TobsArray[0], TobsArray[1], TobsArray[2], TobsArray[3], TobsArray[4], TobsArray[5] )
#Tobs = float(allList[6][0])

M = k * sqrt( (mu)/(a**3) ) * ( Tobs - Tp )
#M = M%(2*pi)
E = NewtonMethod( e, M )
#E = E%(2%pi)
xbar = a*(cos(E)-e)
ybar = a*(1-e**2)**0.5*sin(E)
zbar = 0
r = a*(1-e*cos(E))
nu = acos( (a*(cos(E)-e))/r )
if ybar < 0:
    nu = 2*pi - nu

rvector = vector(xbar, ybar, zbar)

r1 = rotate(rvector, angle = w, axis = (0,0,1) )
r1 = rotate(r1, angle = i, axis = (1,0,0) )
recliptic = rotate(r1, angle = omega, axis = (0,0,1) )

times = [Tobs, Tobs]
Rvector = getBigRVector( times, ephem )
Rvector = newR( Rvector, times )
Rvector1 = -Rvector[0]
req = rotate(recliptic, angle = toRadians(epsilon), axis = (1,0,0) )
pvector = req - Rvector1


punitvector = norm(pvector)

dec = asin(punitvector.z)
#ra1 = asin((punitvector.y)/(cos(dec)))
ra1 = acos((punitvector.x)/(cos(dec)))

ra = ra1
if punitvector.y < 0:
    ra = 2*pi-ra1

ra = ra%(2*pi)

raDec = toDegrees(ra)/15
decDec = toDegrees(dec)
if decDec < 0:
    mult = -1
    decDec = -1*decDec
else:
    mult = 1

raArray = [ int(raDec), int((raDec-int(raDec))*60), round((((raDec-int(raDec))*60) - int((raDec-int(raDec))*60))*60,4)]
decArray = [ mult*int(decDec), int((decDec-int(decDec))*60), round((((decDec-int(decDec))*60)-int((decDec-int(decDec))*60))*60,4)]




    
