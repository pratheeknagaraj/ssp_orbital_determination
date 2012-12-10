'''
Pratheek Nagaraj
July 18, 2011
Observation Data

This program will generate the results of the observation data
'''

# -- Import Packages -- #
from visual import *
from numpy import *
from math import *
import pyfits
import numdisplay

#################################### Functions ########################################

def convertDegDectoHMS( degrees, default = "HMS" ):
     
     if degrees > 0:
          hours = int(degrees)
     else:
          hours = int(degrees)
          degrees = abs(degrees)
     minutes = int(degrees*60%60)
     seconds = round(degrees*3600%3600%60,6)

     if default == "HMS":
          return hours, minutes, seconds
     elif default == "MS":
          return minutes + 60*hours, seconds
     elif default == "S":
          return seconds + minutes*60 + hours*3600

def convertDegDectoDMS( degreesDec, default = "DMS" ):
     
     if degreesDec > 0:
          degrees = int(degreesDec)
     else:
          degrees = int(degreesDec)
          degreesDec = abs(degreesDec)
     minutes = int(degreesDec*60%60)
     seconds = round(degreesDec*3600%3600%60,6)

     if default == "DMS":
          return degrees, minutes, seconds
     elif default == "MS":
          return minutes + 60*degrees, seconds
     elif default == "S":
          return seconds + minutes*60 + degrees*3600

def centroidFunc( ):

    # -- Display the Object Image -- #
    #Open the object image and display using DS9
    numdisplay.display(objectImage)

    image = noise( objectImage )

    choice = raw_input("Would you like to create an output centroid file (Y/N)? ")
    outputChoice = False
    fileOut = None
    if choice == "y" or choice == "Y" or choice == "Yes" or choice == "yes":
         
         fileName3 = raw_input("Please enter file name of centroid outputs: ")
         fileOut = open(fileName3 + ".txt", "w" )
         
         fileOut.write( "Centroid Report\n\n" )
         outputChoice = True

    choice2 = raw_input("Would you like to import search data from a file (Y/N)? ")
    if choice2 == "y" or choice2 == "Y" or choice2 == "Yes" or choice2 == "yes":
         centroidProcess2( image, outputChoice, fileOut )
    else:
         centroidProcess1( image, outputChoice, fileOut, True )

    if outputChoice == True:
         fileOut.close()

def centroidProcess2( image, outputChoice, fileOut ):
     filePath2 = raw_input("Please enter the file name (without .txt): ")
     fileIn = open( filePath2 + ".txt", "r" )
     print ""

     xlist = []
     ylist = []
     radii = []

     for line in fileIn:
          newLine = line.split()
           
          xlist.append(float(newLine[0]))
          ylist.append(float(newLine[1]))
          radii.append(float(newLine[2]))

     for i in range(0, len(xlist) ):
         brightLoc = image[ ylist[i]-radii[i]-1:ylist[i]+radii[i], xlist[i]-radii[i]-1:xlist[i]+radii[i] ]

         # -- Find Centroid -- #
         #Call Centroid Function to locate the given centroid'
         centroidArray = calcCentroid( brightLoc )
         centroid = ( xlist[i] + centroidArray[0], ylist[i] + centroidArray[1] )
         uncertainty = centroidArray[2]

         if outputChoice == True:
              if i == 0:
                   fileOut.write( "Main Object: " )
              fileOut.write( str(xlist[i] + centroidArray[0]) + " " + str(ylist[i] + centroidArray[1]) + " " + str(uncertainty) + "\n")
              if i != len(xlist) - 1:
                   fileOut.write( "Reference Object: " )

def centroidProcess1( image, outputChoice, fileOut, start = False ): ########## User Interface

    if start == True:
         print "Main Object"
         if outputChoice == True:
              fileOut.write( "Main Object: " )
     
    # -- Bright Pixel -- #
    #Get bright pixel location
    xPos = input("Please enter the x value for the location of the object pixel: ")
    yPos = input("Please enter the y value for the location of the object pixel: ")
    #radius corresponds to the numer of points to be used in the calculation i.e. 1->9 and 3->49"
    radius = input("Please enter the radius for the location of the object pixel: ")

    brightLoc = image[ yPos-radius-1:yPos+radius, xPos-radius-1:xPos+radius ]

    # -- Find Centroid -- #
    #Call Centroid Function to locate the given centroid'
    centroidArray = calcCentroid( brightLoc )
    centroid = ( xPos + centroidArray[0], yPos + centroidArray[1] )
    uncertainty = centroidArray[2]

    # -- Output -- #
    print "The centroid of the object is: " + str(centroid)
    print "The uncertainty of the calculation is: " + str(uncertainty)

    #Give choice for a visual'
    choice = raw_input("Would you like to see a visual (Y/N)? ")
    if choice == "y" or choice == "Y" or choice == "Yes" or choice == "yes":
        visualFunc( brightLoc, centroidArray[0], centroidArray[1], uncertainty )

    if outputChoice == True:
         fileOut.write( str(xPos + centroidArray[0]) + " " + str(yPos + centroidArray[1]) + " " + str(uncertainty) + "\n")
    
    # -- Iterate -- #
    choice2 = raw_input("Would you like to do another centroid calculation (Y/N)? ")
    if choice2 == "y" or choice2 == "Y" or choice2 == "Yes" or choice2 == "yes":
        print "\nReference Object"
        fileOut.write( "Reference Object: " )
        centroidProcess1( image )
    
def calcCentroid( array ):  ########## Mathematical Computation
    # -- Total Sum -- #
    #Use Python command to sum all elements in the array
    sumArray = sum(array);
    
    # -- X Coordinate -- #
    #Create the sum of the X Weighted Values'
    xSum = 0.0

    #Loop through the array and weight the X Values 
    weights = range( -len(array)/2+1, len(array)/2+1 )
    for element in array:
        pos = 0
        for point in element:
            xSum = xSum + weights[pos] * point
            pos = pos + 1

    #Find the X Coordinate by dividing by the total sum
    xCoor = xSum / sumArray
    xCoor = round(xCoor,3)

    # -- Y Coordinate -- #
    #Create the sum of the Y Weighted Values'
    ySum = 0.0

    #Loop through the array and weight the Y Values
    weight = -(len(array)/2)
    for element in array:
        ySum = ySum + weight*sum(element)
        weight = weight + 1

    #Find the Y Coordinate by dividing by the total sum
    yCoor = ySum / sumArray
    yCoor = round(yCoor,3)

    uncertainty = 0.0
    for element in array:
        for element2 in element:
            uncertainty = uncertainty + 2 * sqrt(element2)

    uncertainty = uncertainty / sumArray
    uncertainty = round( uncertainty, 3 )

    return xCoor, yCoor, uncertainty

def noise( normalizedObject ):  ##########  Noise Floor Portion
    # -- Noise Floor -- #
    #Get boundaties and adjust image'
    x1Bound = input("Please enter left x bound for the noise floor: ")
    x2Bound = input("Please enter right x bound for the noise floor: ")
    y1Bound = input("Please enter lower y bound for the noise floor: ")
    y2Bound = input("Please enter upper y bound for the noise floor: ")
    print "---------------------------------------------"

    subImage = normalizedObject[ y1Bound - 1: y2Bound, x1Bound - 1: x2Bound ]
    meanValue = subImage.mean()

    #New image modified with floor
    image = normalizedObject - meanValue
    image[where(image < 0)]=0
    return image

def visualFunc( array, x, y, unc ):  ########## Visualization
    # -- Display a visual of the centroid calculation -- #
    #Create brightness circles with nested for loop

    maxValue = array.max()
    print maxValue
    
    posY = -len(array)/2
    for row in array:
        posX = -len(array)/2
        for element in row:
            sphere( radius = element/(maxValue*2), pos = (posX, posY), color = color.blue )
            posX = posX + 1
        posY = posY + 1

    #Create circle for centroid
    sphere(radius = unc, pos = ( x, y, 1 ), color = color.red )

def minorArray( inArray, column ):
    #Get Minor Array
    'Split to make minor array based on column'
    front = inArray[1:,0:column]
    back = inArray[1:,column+1:]
    newarray = append(front,back,axis = 1)
    return newarray

def determinant( inArray ):
    #Get the Determinant
    det = 0
    'Loop through array for determinant'
    for i in range(0,len(inArray)):
        if len(inArray) == 2:
            'Minor is a 2x2'
            det = inArray[0,0]*inArray[1,1] - inArray[1,0]*inArray[0,1]
            return det
        else:
            'Recursively get minor'
            det = det + inArray[0,i]*(-1)**i * determinant( minorArray( inArray, i ) )
    
    return det

def cramerSlice( inArray, column ):
    #Get Minor Array
    'Split to make minor array based on column'
    front = inArray[:,0:column]
    mid = inArray[:,len(inArray):len(inArray)+1]
    back = inArray[:,column+1:]
    'Append arrays onto each other'
    frontmid = append(front,mid,axis = 1)
    newarray = append(frontmid,back,axis = 1)
    finalarray = newarray[0:len(inArray),0:len(inArray)]
    return newarray

def cramer( inArray1 , inArray2, GaussUse ):
    
    inArray=array([])
    inArray=zeros((len(inArray2),len(inArray2)+1))
    for i in range( 0,len(inArray2) ):
        for j in range( 0, len(inArray2) ):
            inArray[i,j] = inArray1[i][j]
        inArray[i,len(inArray2)]=inArray2[i]
    #Input system of equation to do Cramer's Rule

    if GaussUse or len(inArray2) > 6:
        inArray = Gauss( inArray )
        solutionArray = GaussSolve( inArray )
    else:
        Dbase = 0.0
        'Get D for Cramers'
        Dbase = determinant( cramerSlice( inArray, len(inArray) ) )
        if Dbase == 0:
            return "The coefficient array is singular"
        
        solutionArray = array([])
        solutionArray = zeros(len(inArray))
        'Get D of sub for Cramers'
        for i in range(0,len(inArray)):
            Dsub = determinant( cramerSlice( inArray, i ) )
            solutionArray[i] = ( Dsub/Dbase )
    
    return solutionArray

def Gauss( inArray ):

    for i in range(0,len(inArray)):
        pivotCol = i
        diagElement = 0.0
        diagElement = inArray[pivotCol,pivotCol]
        for j in range(i,len(inArray)+1):
            inArray[pivotCol,j] = inArray[pivotCol,j]/diagElement
        diagElement = inArray[pivotCol,pivotCol]
        for k in range(pivotCol+1,len(inArray)):
            point = 0.0
            point = inArray[k,pivotCol]
            multiple = point/diagElement
            for l in range(pivotCol,len(inArray)+1):
                inArray[k,l] = inArray[k,l] - multiple * inArray[pivotCol,l]
    return inArray

def GaussSolve( inArray ):

    solutions = []
    for i in range(0,len(inArray)):
        solutions.append(1)
        
    for i in arange(len(inArray)-1, -1, -1 ):
        value = inArray[i][len(inArray)]
        for j in range(i+1, len(inArray)):
            value = value - inArray[i][j]*solutions[j]
        solutions[i] = value
        
    return solutions

def LSPR():

    choice = raw_input("Would you like to import text data (Y/N)?: ")
    
    if choice == "Y" or choice == "y" or choice == "Yes" or choice == "yes":
        filePath = raw_input("Please enter the file name (without .txt): ")
        fileIn = open( filePath + ".txt", "r" )
        print ""

        ra = []
        dec = []
        xlist = []
        ylist = []
        column = 0

        for line in fileIn:
            column = 0
            newLine = line.split()
            
            data = []

            if newLine and newLine[0] != "//":

                newLine = [i for i in newLine if i != '|']
                newLine = [i for i in newLine if i != '||']

                if len(newLine) != 3:
                    for i in range(len(newLine)-9,len(newLine)):
                        data.append( float( newLine[i] ) )
                        
                    xlist.append( data[0] )
                    ylist.append( data[1] )

                    if data[3] < 0:
                        ra.append( data[3] - data[4]/60.0 - data[5]/3600.0 )
                    else:
                        ra.append( data[3] + data[4]/60.0 + data[5]/3600.0 )
                        
                    if data[6] < 0:
                        dec.append( data[6] - data[7]/60.0 - data[8]/3600.0 )
                    else:
                        dec.append( data[6] + data[7]/60.0 + data[8]/3600.0 )
                elif len(newLine) == 3:
                    for element in newLine:
                        data.append( float( element ) )
                    xObject = data[0]
                    yObject = data[1]
                    
        num = len(ra)

               
    else:
        num = input("Please enter the number of reference stars: ")
        print "\n- Please Input Data -\n"
        ra = []
        dec = []
        xlist = []
        ylist = []

        print "\n- Asteroid Inputs -\n"
        xObject = float(input("Please enter the x coordinate of the object: " ))
        yObject = float(input("Please enter the y coordinate of the object: " ))

        for i in range(0,num):
            xlist.append( input("\nPlease enter the x coordinate of reference star number " + str(i+1) + ": " ) )
            ylist.append( input("Please enter the y coordinate of reference star number " + str(i+1) + ": " ) )
            
            rah = input("Please enter the RA of reference star number " + str(i+1) + " (hours): ")
            ram = input("Please enter the RA of reference star number " + str(i+1) + " (minutes): ")
            ras = input("Please enter the RA of reference star number " + str(i+1) + " (seconds): ")

            decd = input("Please enter the DEC of reference star number " + str(i+1) + " (degrees): ")
            decm = input("Please enter the DEC of reference star number " + str(i+1) + " (minutes): ")
            decs = input("Please enter the DEC of reference star number " + str(i+1) + " (seconds): ")
            
            if i != num - 1:
                print "\n- Next Reference Star -"

            if rah < 0:
                ra.append( rah - ram/60.0 - ras/3600.0 )
            else:
                ra.append( rah + ram/60.0 + ras/3600.0 )
                
            if dec < 0:
                dec.append( decd - decm/60.0 - decs/3600.0 )
            else:
                dec.append( decd + decm/60.0 + decs/3600.0 )

    xsum = sum(xlist)
    xsqsum = sum([x**2 for x in xlist])
    xysum = sum([x*y for x, y in zip(xlist,ylist)])
    ysum = sum(ylist)
    ysqsum = sum([y**2 for y in ylist])
    asum = sum(ra)
    axsum = sum([a*x for a, x in zip(ra,xlist)])
    aysum = sum([a*y for a, y in zip(ra,ylist)])
    dsum = sum(dec)
    dxsum = sum([d*x for d, x in zip(dec,xlist)])
    dysum = sum([d*y for d, y in zip(dec,ylist)])
 
    coeff1 = [[num, xsum, ysum],[xsum, xsqsum, xysum],[ysum,xysum,ysqsum]]
    const1 = [asum, axsum, aysum]
    coeff2 = [[num, xsum, ysum],[xsum, xsqsum, xysum],[ysum,xysum,ysqsum]]
    const2 = [dsum, dxsum, dysum]

    sol1 = cramer( coeff1, const1, False )
    sol2 = cramer( coeff2, const2, False )

    RAObject = sol1[0] + sol1[1]*xObject + sol1[2]*yObject
    DecObject = sol2[0] + sol2[1]*xObject + sol2[2]*yObject

    if RAObject > 0:
        RAObj_h = int(RAObject)
    else:
        RAObj_h = int(RAObject)
        RAObject = abs(RAObject)
    RAObj_m = int(RAObject*60%60)
    RAObj_s = round(RAObject*3600%3600%60,4)

    if DecObject > 0:
        DecObj_d = int(DecObject)
    else:
        DecObj_d = int(DecObject)
        DecObject = abs(DecObject)
    DecObj_m = int(DecObject*60%60)
    DecObj_s = round(DecObject*3600%3600%60,4)
    
    print "RA of Object: ", RAObj_h, RAObj_m, RAObj_s
    print "Dec of Object: ", DecObj_d, DecObj_m, DecObj_s

    RARes = []
    for i in range(num):
        RARes.append( ( ra[i] - ( sol1[0] + sol1[1]*xlist[i] + sol1[2]*ylist[i] ) ) )
        
    DecRes = []
    for i in range(num):
        DecRes.append( ( dec[i] - ( sol2[0] + sol2[1]*xlist[i] + sol2[2]*ylist[i] ) ) )

    print "\n- Star Residuals -\n"
    print "Star #  RA (S)     Dec (S)"
    for i in range(num):
        print "Star " + str(i+1) + ":%10.6f" % convertDegDectoHMS(RARes[i], "S") + \
              "%10.6f" % convertDegDectoDMS(DecRes[i], "S")

    print "\n- Standard Deviation -\n"

    RAsq = sum([a**2 for a in RARes])
    Decsq = sum([d**2 for d in DecRes])

    RAstd = round(pow(RAsq / ( num - 3 ), 0.5), 8)
    Decstd = round(pow(Decsq / ( num - 3 ), 0.5), 8)

    print "RA Standard Deviation: ", RAstd
    print "Dec Standard Deviation: ", Decstd

    covariance = (sum([a*d for a, d in zip(RARes,DecRes)]))/num
    print "\n\nCovariance: ", covariance

    choice2 = raw_input("\nWould you like an output file (Y/N)?: ")
    
    if choice2 == "Y" or choice2 == "y" or choice2 == "Yes" or choice2 == "yes":
        filePath = raw_input("Please enter the output file name (without .txt): ")
        fileOut = open( filePath + ".txt", "w" )

        fileOut.write("LSPR Results")
        fileOut.write("\n------------")
                      
        fileOut.write("\nRA of Object: " + str(RAObj_h) + " " + str(RAObj_m) + " " + str(RAObj_s) )
        fileOut.write("\nDec of Object: "+ str(DecObj_d) + " " + str(DecObj_m) + " " + str(DecObj_s) )

        fileOut.write("\n\n- Star Residuals -\n")
        fileOut.write("\nStar #  RA (S)     Dec (S)")
        for i in range(num):
            fileOut.write("\nStar " + str(i+1) + ":%10.6f" % convertDegDectoHMS(RARes[i],"S") + \
                          "%10.6f" % convertDegDectoDMS(DecRes[i],"S") )

        fileOut.write("\n\n- Standard Deviation -\n")
        fileOut.write("\nRA Standard Deviation: " + str(RAstd) )
        fileOut.write("\nDec Standard Deviation: " + str(Decstd) )

        fileOut.write("\n\nCovariance: " + str(covariance))

        fileOut.close()

##################################### Main Frame ######################################### 

# -- Input -- #
#Open the object image'
global objectImage

choiceFirst = raw_input("Do you want to calculate centroids (Y/N)? " )

if choiceFirst == "y" or choiceFirst == "Y" or choiceFirst == "Yes" or choiceFirst == "yes":

    objectImage = pyfits.getdata("2011_07_24_13.fit")
    centroidFunc( )
    
choiceSecond = raw_input("Do you want to perform LSPR (Y/N)? " )

if choiceSecond == "y" or choiceSecond == "Y" or choiceSecond == "Yes" or choiceSecond == "yes":
     LSPR()
else:
     print "Goodbye!"
