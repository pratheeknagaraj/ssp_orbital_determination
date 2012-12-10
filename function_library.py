##################################################################################################################
##### ------------------------------------------------------------------------------------------------------ #####
'''
                                                    Pratheek Nagaraj
                                                    July 26, 2011
                                                    Orbit Determination

                                    Helper library for Orbital Determination Project
'''
##### ------------------------------------------------------------------------------------------------------ #####
##################################################################################################################

#### Imports #####

from visual import *
from math import *
from numpy import *
from ephemPy import Ephemeris as Ephemeris_BC

##### Constants #####

k = 0.01720209895
epsilon = 23.43799 #degrees
c = 0.00200398707

##### Variables #####

latitude = 34.44556
longitude = -119.66167
rEarth = 4.26349283e-5
LST2011July = 17.620083

################################################## ------------- #################################################
################################################## - Functions - #################################################
################################################## ------------- #################################################

#- Conversion Functions -#

def convertHMStoDegDec( hours, minutes, seconds ): #Convert HH:MM:SS to Decimal Degrees
     if hours < 0:
          return (hours - minutes/60.0 - seconds/3600.0) * 15
     else:
          return (hours + minutes/60.0 + seconds/3600.0) * 15

def convertDMStoDegDec( degrees, minutes, seconds ): #Convert DD:MM:SS to Decimal Degrees
     if degrees < 0:
          return degrees - minutes/60.0 - seconds/3600.0
     else:
          return degrees + minutes/60.0 + seconds/3600.0

def convertDegDectoHMS( degrees, default = "HMS" ):  #Convert Decimal Degrees to HH:MM:SS    
     if degrees > 0:
          hours = int(degrees)
     else:
          hours = int(degrees)
          degrees = abs(degrees)
     minutes = int(degrees*60%60)
     seconds = round(degrees*3600%3600%60,4)

     #Custom Output Settings
     if default == "HMS":
          return hours, minutes, seconds
     elif default == "MS":
          return minutes + 60*hours, seconds
     elif default == "S":
          return seconds + minutes*60 + hours*3600

def convertDegDectoDMS( degreesDec, default = "DMS" ):  #Convert Decimal Degrees to DD:MM:SS
     
     if degrees > 0:
          degrees = int(degreesDec)
     else:
          degrees = int(degreesDec)
          degreesDec = abs(degreesDec)
     minutes = int(degreesDec*60%60)
     seconds = round(degreesDec*3600%3600%60,4)

     #Custom Output Settings
     if default == "DMS":
          return degrees, minutes, seconds
     elif default == "MS":
          return minutes + 60*degrees, seconds
     elif default == "S":
          return seconds + minutes*60 + degrees*3600

def toRadians( deg ):  #Convert Degrees to Radians
     return deg * pi / 180

def toDegrees( rad ):  #Convert Radians to Degrees
     return rad * 180 / pi

def metersToAU( meters ):  #Convert Meters to Astronomical Units
    return 6.68458134e-12*meters

def AUToMeters( AU ):  #Convert Astronomical Units to Meters
    return 149598000000*AU

def convertToDays( number, default = "Y", default2 = "D" ):
     if default == "Y":
          if default2 == "D":
               return number*365.242199

def listToVector( listIn ):  #Convert a List to a vector
     vectorOut = vector(listIn[0],listIn[1],listIn[2])
     return vectorOut

#- Julian Date -#

def isLeap( year ):

    if year % 400 == 0:
        return True
    elif year % 100 == 0:
        return False
    elif year % 4 == 0:
        return True
    else:
        return False

def JDInput():
    year = input("Please enter year (2000,2001,...): ")
    month = input("Please enter month (1,2,3,...,12): ")
    day = input("Please enter day (1,2,3,...): ")
    hours = input("Please enter hours (0,1,2,...23): ")
    minutes = input("Please enter minutes (0,1,...,59): ")
    seconds = input("Please enter seconds (0,1,...,59): ")

    print JD( year, month, day, hours, minutes, seconds )

def JD( year, month, day, hours, minutes, seconds ):

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

    JDate = JDN + (hours-12)/24. + minutes/1440. + seconds/86400.
    return JDate

def getLST( JD ):
     difference = JD - 2455743.791667
     LST = LST2011July + (24*1.00273791*difference)%24
     return LST

def newR( RVector, times ):
     LST = []
     for i in times:
          LST.append( getLST( i ) )
     vectors = []
     for i in range(len(LST)):
          LSTi = toRadians(convertHMStoDegDec( LST[i], 0, 0 ))
          g = rEarth*vector( cos(toRadians(latitude))*cos(LSTi), cos(toRadians(latitude))*sin(LSTi), sin(toRadians(latitude)) )
          vectors.append( RVector[i] + g )
     return vectors
     
#- Input Functions -#

def readInput():  #Read Input Data
     #Get File Path
     filePath = raw_input("Please enter the file name (without .txt): ")
     fileIn = open( filePath + ".txt", "r" )
     print ""

     ra = []
     dec = []
     time = []
     obs = 0

     for line in fileIn:
          obs = obs + 1
          column = 0
          newLine = line.split()

          ra.append( convertHMStoDegDec( float(newLine[0]), float(newLine[1]), float(newLine[2]) ) )
          dec.append( convertDMStoDegDec( float(newLine[3]), float(newLine[4]), float(newLine[5]) ) )
          time.append( float(newLine[6]) )

     return ra, dec, time, obs

def getBigRVector( time, ephem ):  #Get the Big R Vector
     RVector = []
     for element in time:
          RVector.append(vector(ephem.position(element, 10, 2)))
     return RVector

#- Linear Algebra Functions -#

def tripleScalarProduct( vector1, vector2, vector3 ):  #Perform a triple scalar product
     return dotProduct( crossProduct( vector1, vector2 ), vector3 )

def dotProduct( vector1, vector2 ):  #Calculate the dot product of a given vector
    if len(vector1) != len(vector2):
        return "Vectors are not same dimensions, please enter the same dimension vectors and try again.\nGoodbye!"
    
    dotValue = 0.0
    for i in range( 0, len(vector1) ):
        dotValue = dotValue + vector1[i]*vector2[i]

    return dotValue

def crossProduct( vector1, vector2 ):  #Calculate the cross product of a given vector, user choice to display
    if len(vector1) != 3 or len(vector2) != 3:
        return "Vectors must be three dimensional and the same number of dimensions, please try again.\nGoodbye!"

    matrix = zeros(( 3,3 ))
    matrix[0] = (1,1,1)
    matrix[1] = vector1
    matrix[2] = vector2

    vectorNew = ( determinant( minorArray( matrix, 0 ) ), -determinant( minorArray( matrix, 1 ) ), determinant( minorArray( matrix, 2 ) ) )

    return vectorNew

def minorArray( inArray, column ):  #Get Minor Array
    front = inArray[1:,0:column]
    back = inArray[1:,column+1:]
    newarray = append(front,back,axis = 1)
    return newarray

def determinant( inArray ):  #Get the Determinant
    det = 0
    for i in range(0,len(inArray)):
        if len(inArray) == 2:
            det = inArray[0,0]*inArray[1,1] - inArray[1,0]*inArray[0,1]
            return det
        else:
            det = det + inArray[0,i]*(-1)**i * determinant( minorArray( inArray, i ) )
    
    return det

#- Compute Functions -#

def getPHats( ra, dec ):  #Get the rho hat vectors
     phat = []

     for alpha, delta in zip(ra, dec):
          phat.append( vector( cos(toRadians(alpha))*cos(toRadians(delta)), \
                               sin(toRadians(alpha))*cos(toRadians(delta)), sin(toRadians(delta)) ) )

     return phat

def getTaus( time ):  #Get modified time variables
     taus = []
     for i in range( 0, len(time) ):
          if i != len(time)/2:
               if i >  len(time)/2:
                    value = k*( time[i] - time[len(time)/2] )
                    taus.append( value )
               else:
                    value = k*( time[i] - time[len(time)/2] )
                    taus.append( value )

     return taus

def tripleScalarProductFunction( RVector, phat ):  #Compute the Triple Scalar Product
     
     '''
     R1xp2dotp3, R2xp2dotp3, R3xp2dotp3, p1xR1dotp3, p1xR2dotp3, p1xR3dotp3,
     p2xR1dotp1, p2xR2dotp1, p2xR3dotp1, p1xp2dotp3, p2xp3dotp1
     '''
     
     tripleProducts = []
     tripleProducts.append(tripleScalarProduct( RVector[0], phat[1], phat[2] ))
     tripleProducts.append(tripleScalarProduct( RVector[1], phat[1], phat[2] ))
     tripleProducts.append(tripleScalarProduct( RVector[2], phat[1], phat[2] ))
     tripleProducts.append(tripleScalarProduct( phat[0], RVector[0], phat[2] ))
     tripleProducts.append(tripleScalarProduct( phat[0], RVector[1], phat[2] ))
     tripleProducts.append(tripleScalarProduct( phat[0], RVector[2], phat[2] ))
     tripleProducts.append(tripleScalarProduct( phat[1], RVector[0], phat[0] ))
     tripleProducts.append(tripleScalarProduct( phat[1], RVector[1], phat[0] ))
     tripleProducts.append(tripleScalarProduct( phat[1], RVector[2], phat[0] ))
     tripleProducts.append(tripleScalarProduct( phat[0], phat[1], phat[2] ))
     tripleProducts.append(tripleScalarProduct( phat[1], phat[2], phat[0] ))
                           
     return tripleProducts

def calculateRVectors( RVector, p, phat ):  #Recalculate the r vectors
     rvector = []
     for i in range( 0, len(phat) ):
          rvector.append( p[i]*phat[i] - RVector[i] )
     return rvector

def calculateFandGSmall( rvector, taus ):  #Perform the small f and g series
     f = []
     g = []
     middle = len(rvector)/2
     for i in range(0, len(taus)):
          fValue = 1. - (1./(2*((mag(rvector[middle]))**3)))*((taus[i])**2)
          gValue = taus[i] - (1./(6*((mag(rvector[middle]))**3)))*((taus[i])**2)
          f.append(fValue)
          g.append(gValue)
     return f, g

def calculateFandGLarge( rvector, rdotvector, taus ):  #Perform a larger f and g series
     f = []
     g = []
     middle = len(rvector)/2
     rMid = rvector[middle]
     rDMid = rdotvector
     for t in taus:
          fValue = (1. - (1./(2*((mag(rMid))**3)))*(t**2) + ((dot( rMid, rDMid ))/(2.0*(mag(rMid)**5)))*t**3 \
                   + (1./24.)*(((3.0)/((mag(rMid))**3))*(dot(rMid,rDMid)/(mag(rMid)**2)-(1./(mag(rMid))**3)) \
                   - ((15.*(dot(rMid,rDMid))**2)/(mag(rMid)**7))+(1./(mag(rMid)**6)))*(t**4))
          gValue = t - (1./(6*((mag(rMid)**3))))*(t**3)+((dot(rMid,rDMid))/(4.*(mag(rMid))**5))*(t**4)
          f.append(fValue)
          g.append(gValue)
     return f, g

def recalculateA( f, g ):  #Recalculate the a1 and a3 values
     base = (f[0]*g[1]-f[1]*g[0])*1.0
     a1 = (g[1])/base
     a3 = -(g[0])/base
     return a1, a3

def correctTime( p, taus, time ):  #Light time correction function
    pList = []
    for i in range(len(p)):
        pList.append(p[i])
    timeSubtract = [ (element/c)/86400.0 for element in pList ]
    timeArray = []
    for i in range(len(time)):
        timeArray.append(time[i]-timeSubtract[i])
    taus = getTaus( timeArray )
    return taus, timeArray

def calculateRDotVector( f, g, rvector ):  #Calculate the r dot vector
     base = g[0]*f[1]-g[1]*f[0]
     rdotvector = ((f[1]*1.0)/(base))*rvector[0] - ((f[0]*1.0)/(base))*rvector[2]
     return rdotvector

#- Orbital Elements Function -#

def orbitalElements( rvector, rdotvector, timeArray ):  #Compute the orbital elements
     #Convert to ectors
     rvector = listToVector( rvector[len(rvector)/2] )
     rdotvector = listToVector( rdotvector )

     #Rotate to ecliptic coordinates
     rvector = rotate(rvector, angle = toRadians(-epsilon), axis = (1,0,0) )
     rdotvector = rotate(rdotvector, angle = toRadians(-epsilon), axis = (1,0,0) )

     #Initial values
     mu = 1
     t = timeArray[len(timeArray)/2]
     magr = mag( rvector )
     magrdot = mag( rdotvector )

     #Compute the semi-major axis (a)
     a = 1/((2/magr)-((magrdot**2)/mu))

     #Get the h vector
     hvector = cross( rvector, rdotvector )
     magh = mag( hvector )
     rdotcrossh = listToVector( cross( rdotvector, hvector ) )

     #Compute the eccentricity (e)
     evector = ((rdotcrossh)/mu) - (rvector/magr)
     e = mag( evector )

     #Get the k vector
     kvector = vector(0,0,1)
     hvectordotk = dot( kvector, hvector )

     #Compute the inclination (i)
     i = toDegrees( acos( hvectordotk/magh ) )

     #Get the n vector
     nvector = listToVector( cross( kvector, hvector ) )
     magn = mag( nvector )

     #Get the i vector
     ivector = (1,0,0)
     idotn = dot( ivector, nvector )

     #Get the longitude of ascending node (omega)
     omega = toDegrees( acos( idotn/magn ) )
     if nvector.y < 0:
         omega = 360 - omega

     #Get the argument of pericenter (littleOmega)
     edotn = dot( evector, nvector )
     littleOmega = toDegrees( acos( ( dot( evector, nvector ) )/((e)*(magn)) ) )
     if evector.z < 0:
         littleOmega = 360 - littleOmega

     #Compute the eccentric anomaly (E)
     E = toDegrees(acos((1-(magr/a))/e))
     if dotProduct(rvector, rdotvector) < 0:
         E = 360 - E

     #Compute the mean anomaly (M)
     M = toDegrees((toRadians(E) - e*sin(toRadians(E))))

     #Compute the Period (P)
     P = sqrt(a**3)
     P = convertToDays(P)

     #Compute the time of pericenter passage (T)
     T = t - toRadians(M)/(k*sqrt((mu/(a**3))))

     #Output Results
     print "Semimajor Axis: ", a
     print "Eccentricity: ", e
     print "Inclination: ", i
     print "Longitude of Ascending Node: ", omega
     print "Argument of Perihelion: ", littleOmega
     print "Mean anomaly: ", M
     print "Period: ", P
     print "Time of last perihelion passage: ", T

     #User Option
     choice = raw_input("\nWould you like to export the Orbital Element Data (Y/N)?: ")
     if choice == "Y" or choice == "y" or choice == "Yes" or choice == "yes": 
          filePath = raw_input("Please enter the file name (without .txt): ")
          fileOut = open( filePath + ".txt", "w" )

          fileOut.write("Orbital Determination Results")
          fileOut.write("\n-----------------------------")

          fileOut.write("\n\nSemimajor Axis: " + str(a) )
          fileOut.write("\nEccentricity: " + str(e) )
          fileOut.write("\nInclination: " + str(i) )
          fileOut.write("\nLongitude of Ascending Node: " + str(omega) )
          fileOut.write("\nArgument of Perihelion: " + str(littleOmega) )
          fileOut.write("\nMean anomaly: " + str(M) )
          fileOut.write("\nPeriod: " + str(P) )
          fileOut.write("\nTime of last perihelion passage: " + str(T) )

          fileOut.close()
          
          
