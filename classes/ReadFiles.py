"""used for opening various file types. """
import os

def openReturnLines(fullPath):
    openFile = open(fullPath)
    lineList = openFile.readlines()
    openFile.close()
    return lineList

def ReadSpectraCalcLineList(localPath, molecule, isotopeDepth):

    rows = openReturnLines(os.getcwd() + localPath)
    #create a text list of isotopes for comparison
    i = 1
    isotopeList = []
    while i <= isotopeDepth:
        isotopeList.append(str(i))
        i += 1

    #   column names and locations
    id = 0
    isotope = 1
    wavenumber = 2
    intensity = 3
    einsteinA = 4
    airHalfwidth = 5
    selfHalfwidth = 6
    lowerEnergy = 7
    tempExponent = 8
    pressureShift = 9

    lineListDict = {}
    #   the csv files are comma-seperated, seperate into cells
    for row in rows:
        cell = row.split(',')
        if cell[id] == str(molecule) and cell[isotope] in isotopeList:
            #  create a dictionary of the contents without the wavenumber
            contents = {}
            contents['isotope'] = int(cell[isotope])
            contents['intensity'] = float(cell[intensity])
            contents['einsteinA'] = float(cell[einsteinA])
            contents['airHalfWidth'] = float(cell[airHalfwidth])
            contents['selfHalfWidth'] = float(cell[selfHalfwidth])
            contents['lowerEnergy'] = float(cell[lowerEnergy])
            contents['tempExponent'] = float(cell[tempExponent])
            contents['pressureShift'] = float(cell[pressureShift])

         #  add key : value pair to lineListDict, with the wavenumber being the key for quicker lookup. The value will
         #  be a second dictionary to aid in readability of the contents.
            lineListDict[float(cell[wavenumber])] = contents
    if lineListDict == {}:
        lineListDict = False
    return lineListDict

def readQ(ID, isotopeDepth):
    isoDict = {}
    isotope = 1
    while isotope <= isotopeDepth:
        fileName = 'q' + str(ID) + '-' + str(isotope) + '.csv'
        fullPath = os.getcwd() + '/HITRAN/' + fileName
        lines = openReturnLines(fullPath)
        qDict = {}
        for line in lines:
            cell = line.split(',')
            qDict[int(cell[0])] = float(cell[1])
        isoDict[int(isotope)] = qDict
        isotope += 1

    return isoDict