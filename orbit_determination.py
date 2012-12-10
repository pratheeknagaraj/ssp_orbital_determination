##################################################################################################################
##### ------------------------------------------------------------------------------------------------------ #####
'''
                                                    Pratheek Nagaraj
                                                    July 26, 2011
                                                    Orbit Determination

                                    Find the orbital elements using a series of observations
'''
##### ------------------------------------------------------------------------------------------------------ #####
##################################################################################################################

##### Import #####

from visual import *
from math import *
from numpy import *
from ephemPy import Ephemeris as Ephemeris_BC
from Function_Library import *

##### Class Initializations #####

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

############################################## --------------------- #############################################
############################################## - Program Mainframe - #############################################
############################################## --------------------- #############################################

#- Input Data -#

#choice = raw_input("Would you like to import text data (Y/N)?: ")
choice = "Y"
ra = []
dec = []
time = []
timeArray = []
obs = 0
    
if choice == "Y" or choice == "y" or choice == "Yes" or choice == "yes":    
     ra, dec, time, obs = readInput()
else:
     print "Goodbye!"

#- Get p hat, taus, eath sun vector -#

phat = getPHats( ra, dec )
taus = getTaus( time )
RVector = getBigRVector( time, ephem )

#- Parallax Correction (Center of Earth to Observation Location) -#

RVector = newR( RVector, time )

#- Triple Products -#

tP = []
tP = tripleScalarProductFunction( RVector, phat )

#- Initial Guesses for p -#

a1 = ( time[2] - time[1] )/( time[2] - time[0] )
a3 = -( time[0] - time[1] )/( time[2] - time[0] )

#- Calculate p1, p2, p3 -#

p1 = ( a1*tP[0] - tP[1] + a3*tP[2] )/( a1*tP[9] )
p2 = ( a1*tP[3] - tP[4] + a3*tP[5] )/( -1*tP[9] )
p3 = ( a1*tP[6] - tP[7] + a3*tP[8] )/( a3*tP[10] )
p = [p1,p2,p3]

#- Light Time Correction -#
taus, timeArray = correctTime( p, taus, time )

#- Calculate r -#

rvector = calculateRVectors( RVector, p, phat )

#- F and G series -#

f, g = calculateFandGSmall( rvector, taus )

#- Recalculate a1 and a3 -#

a1, a3 = recalculateA( f, g )

#- Full F and G -#

p1 = ( a1*tP[0] - tP[1] + a3*tP[2] )/( a1*tP[9] )
p2 = ( a1*tP[3] - tP[4] + a3*tP[5] )/( -1*tP[9] )
p3 = ( a1*tP[6] - tP[7] + a3*tP[8] )/( a3*tP[10] )
p = [p1,p2,p3]

#- Light Time Correction -#
taus, timeArray = correctTime( p, taus, time )

#- Initialize F and G Vectors -#

rvector = calculateRVectors( RVector, p, phat )
rdotvector = calculateRDotVector( f, g, rvector )

#####- Main Iteration Loop -#####
cont = True
tolerance = 1e-8
while cont == True:
     pastR = rvector[len(rvector)/2]  #Get Previous rvector
     f, g = calculateFandGLarge( rvector, rdotvector, taus )  #Perform f and g series
     a1, a3 = recalculateA(f, g)  #Get new a1 and a3
     p1 = ( a1*tP[0] - tP[1] + a3*tP[2] )/( a1*tP[9] )
     p2 = ( a1*tP[3] - tP[4] + a3*tP[5] )/( -1*tP[9] )
     p3 = ( a1*tP[6] - tP[7] + a3*tP[8] )/( a3*tP[10] )
     p = [p1,p2,p3]  #Get new rho magnitudes
     taus, timeArray = correctTime( p, taus, time )  #Light time corrected taus
     rvector = calculateRVectors( RVector, p, phat )  #New r vector
     rdotvector = calculateRDotVector( f, g, rvector )  #New r dot vector
     difference = [fabs(a - b) for a, b in zip(pastR, rvector[1])]
     if difference[0] < tolerance and difference[1] < tolerance and difference[2] < tolerance:  #Tolerance check
          cont = False
          
#- Convert to orbital elements -#

a, e, inc, omega, littleOmega, E, T, mu = orbitalElements( rvector, rdotvector, timeArray )

################################################## ------------------ ############################################
################################################## - Visualization  - ############################################
################################################## ------------------ ############################################

scene = display( title = "Orbital Determination", width = 800, height = 600, up = (0,0,1) )

sunRadius = metersToAU(6.955e8)
earthRadius = metersToAU(6371000)

scene.range = (mag(RVector[0])*1.5,mag(RVector[0])*1.5,mag(RVector[0])*1.5)

sun = sphere( pos = (0,0,0), color = color.yellow, radius = sunRadius*2.5e1, material = materials.emissive )

xAxis = arrow( pos = (0,0,0), axis=(.25,0,0), shaftwidth = 0.005, color = color.white )
yAxis = arrow( pos = (0,0,0), axis=(0,.25,0), shaftwidth = 0.005, color = color.white )
zAxis = arrow( pos = (0,0,0), axis=(0,0,.25), shaftwidth = 0.005, color = color.white )

ecliptic1 = box( pos = (0,0,0), width = 0.001, length = 5, height = 5, color = color.green, opacity = 0.1 )


EValue = toRadians(E)
distMin = 1.0
t = []
t.append(T)

while True:
    M = (EValue - e*sin(EValue))
    t[0] = T + M/(k*sqrt((mu/(a**3))))
    
    xbar = a*(cos(EValue)-e)
    ybar = a*(sin(EValue))*(1-e**2)**0.5
    zbar = 0
    rvec = vector( xbar, ybar, zbar )

    rvec = rotate( rvec, angle = toRadians(littleOmega), axis = (0,0,1) )
    rvec = rotate( rvec, angle = toRadians(inc), axis = (1,0,0) )
    rvec = rotate( rvec, angle = toRadians(omega), axis = (0,0,1) )

    if EValue == toRadians(E):
        asteroid = sphere( pos = rvec, color = color.white, radius = earthRadius*5e2, material = materials.rough )
        asteroid.orbit = curve( pos = asteroid.pos, color = asteroid.color,  radius = earthRadius*1e2 )
        earth = sphere( pos = -(getBigRVector( t, ephem ))[0], radius = earthRadius*1e3, material = materials.earth )
        earth.orbit = curve( pos = earth.pos, color = color.cyan,  radius = earthRadius*1e2 )
    
    asteroid.pos = rvec
    asteroid.orbit.append( pos = ( asteroid.pos ) )
    earth.pos = -(getBigRVector( t, ephem ))[0]
    earth.orbit.append( pos = ( earth.pos ) )
    
    dist = mag(asteroid.pos - earth.pos)
    if distMin > dist:
        distMin = dist
        print "Closest Approach: ", distMin, " | Time: ", round((t[0] - T)/365.24,4)

    EValue = EValue + 0.001


    
