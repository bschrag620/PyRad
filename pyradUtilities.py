import os
from time import gmtime, strftime
import urllib.request as urlrequest

cwd = os.getcwd()
dataDir = '%s/data' % cwd
curvesDir = '%s/curves' % dataDir
molParamsFile = '%s/molparams.txt' % dataDir
debugFilePath = '%s/logger.txt' % cwd

debugFile = open(debugFilePath, 'w')
debugFile.write(strftime("%Y-%m-%d %H:%M:%S\n", gmtime()))


def outputToLog(text):
    debugFile = open(debugFilePath, 'a')
    debugFile.write('%s\n' % text)
    debugFile.close()


def setupDir():
    print('Verifying data stucture...', end='', flush=True)
    directoryList = [dataDir, curvesDir]
    fileList = []
    directoryCheck = True
    fileCheck = True
    for moleculeID in HITRAN_GLOBAL_ISO:
        for localIsoID in HITRAN_GLOBAL_ISO[moleculeID]:
            globalIsoPath = '%s/%s' % (dataDir, HITRAN_GLOBAL_ISO[moleculeID][localIsoID])
            isotopeFilePath = '%s/param.pyr' % globalIsoPath
            directoryList.append(globalIsoPath)
            fileList.append(isotopeFilePath)
    for directory in directoryList:
        if not os.path.isdir(directory):
            os.makedirs(directory)
    print('directories checked...', end='', flush=True)
    for file in fileList:
        if not os.path.isfile(file):
            if not os.path.isfile(molParamsFile):
                downloadMolParam()
            getMolParamsFromHitranFile()
    print('files checked.')
    return directoryCheck and fileCheck


def openReturnLines(fullPath):
    if not os.path.isfile(fullPath):
        return []
    openFile = open(fullPath)
    lineList = openFile.readlines()
    openFile.close()
    return lineList


def writeDictListToFile(dictionary, fullPath):
    openFile = open(fullPath, 'w')
    for key in dictionary:
        text = '%s,%s\n' % (key, ','.join(str(item) for item in dictionary[key]))
        openFile.write(text)
    openFile.close()


def writeListToFile():
    pass


def getCurves(curveType):
    curveDict = {}
    if curveType == 'voigt':
        curveFilePath = '%s/voigt.pyr' % curvesDir
    rows = openReturnLines(curveFilePath)
    if rows:
        for row in rows:
            cells = row.split(',')
            key = cells.pop(0)
            curveDict[key] = cells
    return curveDict


def getMolParamsFromHitranFile():
    outputToLog('getMolParams')
    rows = openReturnLines(molParamsFile)
    outputToLog('row length=%s' % len(rows))
    isotopeInfo = {}
    i = 0
    while i < len(rows) -1:
        outputToLog('i=%s' % i)
        cells = rows[i].split()
        outputToLog('cells length=%s' % len(cells))
        outputToLog('cells value=%s' % cells)
        while cells[0].lower() in MOLECULE_ID and i < len(rows) - 1:
            moleculeShortName = cells[0].lower()
            moleculeID = int(cells[1].replace(')', '').replace('(', ''))
            i += 1
            isotopeNumber = 1
            cells = rows[i].split()
            while cells[0].lower() not in MOLECULE_ID:
                IsoN = int(cells[0])
                abundance = float(cells[1])
                q296 = float(cells[2])
                gj = int(cells[3])
                molMass = float(cells[4])
                globalID = HITRAN_GLOBAL_ISO[moleculeID][isotopeNumber]
                infoList = [moleculeShortName, moleculeID, IsoN, abundance, q296, gj, molMass]
                isotopeInfo[globalID] = infoList
                writeDictListToFile({globalID: isotopeInfo[globalID]},'%s/%s/params.pyr' % (dataDir, globalID))
                isotopeNumber += 1
                i += 1
                if i == len(rows):
                    break
                cells = rows[i].split()
        i += 1
    return isotopeInfo


def gatherData(globalIsoId, rangeMin, rangeMax):
    rangeList = []
    segment = int(rangeMin / 100) * 100
    info = {}
    while segment < rangeMax:
        rangeList.append(segment)
        segment += 100
    for segment in rangeList:
        globalIsoDir = '%s/%s/' % (dataDir, globalIsoId)
        if not os.path.isdir(globalIsoDir):
            os.makedirs(globalIsoDir)
        filePath = '%s/%s/%s.pyr' % (dataDir, globalIsoId, segment)
        if not os.path.isfile(filePath):
            print('File not found for %s, range %s-%s. Downloading from HITRAN...' % (globalIsoId, segment, segment + 100))
            downloadHitran(filePath, globalIsoId, segment, segment + 100)
        info.update(readHitranOnlineFile(filePath, rangeMin, rangeMax))
        getQData(globalIsoId)
    return info


def getQData(isotope):
    qPath = cwd + '/data/%s/' % isotope
    filePath = qPath + 'q%s.txt' % isotope
    if not os.path.isfile(filePath):
        downloadQData(isotope)
    return readQFile(isotope)


def downloadMolParam():
    print('Missing molecule parameters file. Downloading from http://hitran.org')
    url = 'http://hitran.org/media/molparam.txt'
    try:
        request = urlrequest.urlopen(url)
    except urlrequest.HTTPError:
        print('Can not retrieve file at %s. Exiting.' % url)
        return False
    except urlrequest.URLError:
        print('Can not connect. Exiting')
        return False
    chunkSize = 1024 * 64
    openFile = open(molParamsFile, 'w')
    while True:
        chunk = request.read(chunkSize)
        if not chunk:
            print('Molecule parameters successfully downloaded.')
            break
        openFile.write(chunk.decode('utf-8'))
    openFile.close()


def downloadQData(isotope):
    url = 'http://hitran.org/data/Q/q%s.txt' % str(isotope)
    path = cwd + '/data/%s/q%s.txt' % (isotope, isotope)
    request = urlrequest.urlopen(url)
    print('Downloading Q table from %s' % url)
    openFile = open(path, 'w')
    chunkSize = 1024 * 64
    while True:
        chunk = request.read(chunkSize)
        if not chunk:
            break
        openFile.write(chunk.decode('utf-8'))
    openFile.close()


def downloadHitran(path, globalID, waveMin, waveMax):
    params = 'molec_id,local_iso_id,nu,sw,a,elower,gamma_air,gamma_self,delta_air,n_air'
    url = 'http://hitran.org/lbl/api?iso_ids_list=' + \
          str(globalID) + '&numin=' + str(waveMin) + \
          '&numax=' + str(waveMax) + \
          '&fixwidth=0&sep=[comma]' + \
          '&request_params=%s' % params
    try:
        request = urlrequest.urlopen(url)
    except urlrequest.HTTPError:
        print('Can not retrieve data for given parameters.')
        request = False
    except urlrequest.URLError:
        print('Can not connect to %s' % str(url))
        request = False
    print('Connected to http://hitran.org, beginning download.')
    dirLength = len(cwd)
    if request:
        i = 0
        openFile = open(path, 'w')
        chunkSize = 1024 * 64
        while True:
            i += 1
            chunk = request.read(chunkSize)
            if not chunk:
                break
            openFile.write(chunk.decode('utf-8'))
            print('%s downloaded and written to %s%s' % (chunkSize, path[dirLength:], '.' * i), end='\r', flush=True)
        print('\n', end='\r')
        openFile.close()


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


def readHitranOnlineFile(path, waveMin, waveMax):
    isotope = 1
    wavenumber = 2
    intensity = 3
    einsteinA = 4
    airHalfwidth = 6
    selfHalfwidth = 7
    lowerEnergy = 5
    tempExponent = 9
    pressureShift = 8
    lineListDict = {}
    rows = openReturnLines(path)
    for row in rows:
        cell = row.split(',')
        if waveMin < float(cell[wavenumber]):
            if float(cell[wavenumber]) < waveMax:
                contents = {'isotope': int(cell[isotope]),
                            'intensity': float(cell[intensity]),
                            'einsteinA': float(cell[einsteinA]),
                            'airHalfWidth': float(cell[airHalfwidth]),
                            'selfHalfWidth': float(cell[selfHalfwidth]),
                            'lowerEnergy': float(cell[lowerEnergy]),
                            'tempExponent': float(cell[tempExponent]),
                            'pressureShift': float(cell[pressureShift])}
                lineListDict[float(cell[wavenumber])] = contents
    return lineListDict


def readQFile(isotope):
    path = cwd + '/data/%s/q%s.txt' % (isotope, isotope)
    if not os.path.isfile:
        downloadQData(isotope)
    openFile = open(path, 'r')
    rows = openFile.readlines()
    qDict = {}
    for row in rows:
        cell = row.split()
        qDict[int(cell[0])] = float(cell[1])
    return qDict


def readMolParams(globalIso):
    filePath = '%s/%s/params.pyr' % (dataDir, globalIso)
    params = openReturnLines(filePath)
    row = params[0]
    cells = row.split(',')
    globalIso = int(cells[0])
    shortName = cells[1]
    moleculeNum = int(cells[2])
    isoN = int(cells[3])
    abundance = float(cells[4])
    q296 = float(cells[5])
    gj = int(cells[6])
    molMass = float(cells[7])
    return [globalIso, shortName, moleculeNum, isoN, abundance, q296, gj, molMass]


VERSION = '1.1'
titleLine = "***********************              PyRad              ***********************"
messageGap = int((len(titleLine) - len(VERSION) - 1) / 2)
GREETING = "%s\n" \
           "%sv%s%s\n" \
           "An open-source, amateur attempt at a radiative transfer model for an atmosphere\n\n" \
           "\tAll line lists are downloaded from HITRAN\n" \
           "\tAll information for calculation of lineshapes comes \n" \
           "\tfrom spectralcalc.com and hitran.org \n" \
           "\tFor a more complete option, check out HAPI from HITRAN.\n" \
           "\tThanks to those you have made the information available \n" \
           "\tand so easily accessible.\n\n" \
           "*******************************************************************************" \
           % (titleLine, ' ' * messageGap, VERSION, ' ' * (len(titleLine) - messageGap))

HITRAN_GLOBAL_ISO = {1: {1: 1,
                         2: 2,
                         3: 3,
                         4: 4,
                         5: 5,
                         6: 6,
                         7: 129},
                     2: {1: 7,
                         2: 8,
                         3: 9,
                         4: 10,
                         5: 11,
                         6: 12,
                         7: 13,
                         8: 14,
                         9: 121,
                         10: 15,
                         11: 120,
                         12: 122},
                     3: {1: 16,
                         2: 17,
                         3: 18,
                         4: 19,
                         5: 20},
                     4: {1: 21,
                         2: 22,
                         3: 23,
                         4: 24,
                         5: 25, },
                     5: {1: 26,
                         2: 27,
                         3: 28,
                         4: 29,
                         5: 30,
                         6: 31},
                     6: {1: 32,
                         2: 33,
                         3: 34,
                         4: 35},
                     7: {1: 36,
                         2: 37,
                         3: 38},
                     8: {1: 39,
                         2: 40,
                         3: 41},
                     9: {1: 42,
                         2: 43},
                     10: {1: 44},
                     11: {1: 45,
                          2: 46},
                     12: {1: 47,
                          2: 117},
                     13: {1: 48,
                          2: 49,
                          3: 50},
                     14: {1: 51,
                          2: 110},
                     15: {1: 52,
                          2: 53,
                          3: 107,
                          4: 108},
                     16: {1: 19,
                          2: 11,
                          3: 111,
                          4: 112},
                     17: {1: 56,
                          2: 113},
                     18: {1: 57,
                          2: 58},
                     19: {1: 59,
                          2: 60,
                          3: 61,
                          4: 62,
                          5: 63},
                     20: {1: 64,
                          2: 65,
                          3: 66},
                     21: {1: 67,
                          2: 68},
                     22: {1: 69,
                          2: 118},
                     23: {1: 70,
                          2: 71,
                          3: 72},
                     24: {1: 73,
                          2: 74},
                     25: {1: 75},
                     26: {1: 76,
                          2: 77,
                          3: 105},
                     27: {1: 78,
                          2: 106},
                     28: {1: 79},
                     29: {1: 80,
                          2: 119},
                     30: {1: 126},
                     31: {1: 81,
                          2: 82,
                          3: 83},
                     32: {1: 84},
                     33: {1: 85},
                     34: {1: 86},
                     35: {1: 127,
                          2: 128},
                     36: {1: 87},
                     37: {1: 88,
                          2: 89},
                     38: {1: 90,
                          2: 91},
                     39: {1: 92},
                     40: {1: 93,
                          2: 94},
                     41: {1: 95},
                     42: {1: 96},
                     43: {1: 116},
                     44: {1: 109},
                     45: {1: 103,
                          2: 115},
                     46: {1: 97,
                          2: 98,
                          3: 99,
                          4: 100},
                     47: {1: 114},
                     48: {1: 123},
                     49: {1: 124,
                          2: 125}}


MOLECULE_ID = {'h2o': 1, 'co2': 2, 'o3': 3, 'n2o': 4, 'co': 5,
               'ch4': 6, 'o2': 7, 'no': 8, 'so2': 9,
               'no2': 10, 'nh3': 11, 'hno3': 12, 'oh': 13,
               'hf': 14, 'hcl': 15, 'hbr': 16, 'hi': 17,
               'clo': 18, 'ocs': 19, 'h2co': 20, 'hocl': 21,
               'n2': 22, 'hcn': 23, 'ch3cl': 24, 'h2o2': 25,
               'c2h2': 26, 'c2h6': 27, 'ph3': 28, 'cof2': 29,
               'sf6': 30, 'h2s': 31, 'hcooh': 32, 'ho2': 33,
               'o': 34, 'clono2': 35, 'no+': 36, 'hobr': 37,
               'c2h4': 38, 'ch3oh': 39, 'ch3br': 40, 'ch3cn': 41,
               'cf4': 42, 'c4h2': 43, 'hc3n': 44, 'h2': 45,
               'cs': 46, 'so3': 47, 'c2n2': 48, 'cocl2': 49}


print(GREETING)
setupDir()