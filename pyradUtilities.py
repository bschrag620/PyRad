import os
import sys
import urllib.request as urlrequest
import urllib.error as urlexception
import datetime
import numpy as np


cwd = os.getcwd()
lineSep = os.linesep
dataDir = '%s/data' % cwd
curvesDir = '%s/curves' % dataDir
molParamsFile = '%s/molparams.txt' % dataDir
debuggerFilePath = '%s/logger.txt' % cwd
now = datetime.datetime.now()
debuggerFile = open(debuggerFilePath, 'wb')
debuggerFile.write(bytes('%s\n' % now.strftime("%Y-%m-%d %H:%M:%S"), 'utf-8'))
debuggerFile.close()


def magentaText(text):
    return '%s%s%s' % (TEXT_COLORS['regularMagenta'], text, TEXT_COLORS['colorEnd'])


def cyanText(text):
    return '%s%s%s' % (TEXT_COLORS['regularCyan'], text, TEXT_COLORS['colorEnd'])


def limeText(text):
    return '%s%s%s' % (TEXT_COLORS['regularLime'], text, TEXT_COLORS['colorEnd'])


def underlineCyan(text):
    return '%s%s%s' % (TEXT_COLORS['underlineCyan'], text, TEXT_COLORS['colorEnd'])


def underlineMagenta(text):
    return '%s%s%s' % (TEXT_COLORS['underlineMagenta'], text, TEXT_COLORS['colorEnd'])


def underlineLime(text):
    return '%s%s%s' % (TEXT_COLORS['underlineLime'], text, TEXT_COLORS['colorEnd'])


def underlineWhite(text):
    return '%s%s%s' % (TEXT_COLORS['underlineWhite'], text, TEXT_COLORS['colorEnd'])


def logToFile(text):
    debugFile = open(debuggerFilePath, 'a')
    debugFile.write('%s\n' % text)
    debugFile.close()


def setupDir():
    print('Verifying data structure...', end='')
    sys.stdout.flush()
    directoryList = [dataDir, curvesDir]
    fileList = []
    directoryCheck = True
    fileCheck = True
    for moleculeID in HITRAN_GLOBAL_ISO:
        for localIsoID in HITRAN_GLOBAL_ISO[moleculeID]:
            globalIsoPath = '%s/%s' % (dataDir, HITRAN_GLOBAL_ISO[moleculeID][localIsoID])
            isotopeFilePath = '%s/params.pyr' % globalIsoPath
            directoryList.append(globalIsoPath)
            fileList.append(isotopeFilePath)
    for directory in directoryList:
        if not os.path.isdir(directory):
            os.makedirs(directory)
    print('directories checked...', end='')
    sys.stdout.flush()
    for file in fileList:
        if not os.path.isfile(file):
            if not os.path.isfile(molParamsFile):
                downloadMolParam()
            getMolParamsFromHitranFile()
    print('files checked.')
    return directoryCheck and fileCheck


def openReturnLines(fullPath):
    if not os.path.isfile(fullPath):
        return False
    openFile = open(fullPath)
    lineList = openFile.readlines()
    openFile.close()
    if not lineList or NULL_TAG in lineList[0]:
        return False
    while lineList[0][0] == '#' and len(lineList) > 1:
        lineList.pop(0)
    return lineList

  
def writeDictListToFile(dictionary, fullPath, comments=None, mode='wb'):
    openFile = open(fullPath, mode)
    if comments:
        openFile.write(comments.encode('utf-8'))
    for key in dictionary:
        text = '%s,%s%s' % (key, ','.join(str(item) for item in dictionary[key]), lineSep)
        openFile.write(text.encode('utf-8'))
    openFile.close()


def getCurves(curveType, res):
    curveDict = {}
    resDirectory = '%s/res%s' % (curvesDir, res)
    print('Retrieving %s curves...' % curveType, end='', flush=True)
    if not os.path.isdir(resDirectory):
        os.mkdir(resDirectory)
    if curveType == 'voigt':
        curveFilePath = '%s/voigt.pyr' % resDirectory
    elif curveType == 'lorentz':
        curveFilePath = '%s/lorentz.pyr' % resDirectory
    elif curveType == 'gaussian':
        curveFilePath = '%s/gaussian.pyr' % resDirectory
    rows = openReturnLines(curveFilePath)
    if rows:
        for row in rows:
            cells = row.strip().split(',')
            key = cells.pop(0)
            for i in range(0, len(cells) - 1):
                if cells[i]:
                    try:
                        cells[i] = float(cells[i])
                    except ValueError:
                        print(cells[i])
            cells.pop()
            curveDict[key] = np.asarray(cells)
    print('%s built from cache.' % len(curveDict))
    return curveDict


def getMolParamsFromHitranFile():
    rows = openReturnLines(molParamsFile)
    isotopeInfo = {}
    haveMolecule = False
    for row in rows:
        logToFile('for loop row: %s' % row)
        cells = row.split()
        logToFile('for loop cells: %s' % cells)
        if cells[0].lower() in MOLECULE_ID:
            localIso = 0
            haveMolecule = True
            moleculeShortName = cells[0].lower()
            moleculeID = int(cells[1].replace(')', '').replace('(', ''))
        elif haveMolecule:
            localIso += 1
            IsoN = int(cells[0])
            abundance = float(cells[1])
            q296 = float(cells[2])
            gj = int(cells[3])
            molMass = float(cells[4])
            globalID = HITRAN_GLOBAL_ISO[moleculeID][localIso]
            infoList = [moleculeShortName, moleculeID, IsoN, abundance, q296, gj, molMass]
            isotopeInfo[globalID] = infoList
            writeDictListToFile({globalID: isotopeInfo[globalID]},
                                '%s/%s/params.pyr' % (dataDir, globalID),
                                MOLECULE_PARAM_COMMENTS)
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

# downloads a table of isotopes and data from hitran
def downloadMolParam():
    print('Missing molecule parameters file. Downloading from http://hitran.org')
    url = 'http://hitran.org/media/molparam.txt'
    try:
        request = urlrequest.urlopen(url)
    except urlexception.HTTPError:
        print('Can not retrieve file at %s. Exiting.' % url)
        return False
    except urlexception.URLError:
        print('Can not connect. Exiting')
        return False
    chunkSize = 1024 * 64
    openFile = open(molParamsFile, 'wb')
    while True:
        chunk = request.read(chunkSize)
        if not chunk:
            print('Molecule parameters successfully downloaded.')
            break
        openFile.write(chunk)
    openFile.close()

# downloads q table from hitran
def downloadQData(isotope):
    url = 'http://hitran.org/data/Q/q%s.txt' % str(isotope)
    path = cwd + '/data/%s/q%s.txt' % (isotope, isotope)
    request = urlrequest.urlopen(url)
    print('Downloading Q table from %s' % url)
    openFile = open(path, 'wb')
    chunkSize = 1024 * 64
    while True:
        chunk = request.read(chunkSize)
        if not chunk:
            break
        openFile.write(chunk)
    openFile.close()

# downloads data from hitran using a function found in hapi that's been slightly modified
def downloadHitran(path, globalID, waveMin, waveMax):
    params = 'molec_id,local_iso_id,nu,sw,a,elower,gamma_air,gamma_self,delta_air,n_air'
    url = 'http://hitran.org/lbl/api?iso_ids_list=' + \
          str(globalID) + '&numin=' + str(waveMin) + \
          '&numax=' + str(waveMax) + \
          '&fixwidth=0&sep=[comma]' + \
          '&request_params=%s' % params
    try:
        request = urlrequest.urlopen(url)
    except urlexception.HTTPError:
        print('Can not retrieve data for given parameters.')
        openFile = open(path, 'wb')
        openFile.write(bytes('%s no data available for %s, range %s-%s' %
                             (NULL_TAG, globalID, waveMin, waveMax), 'utf-8'))
        openFile.close()
        request = False
    except urlexception.URLError:
        print('Can not connect to %s' % str(url))
        request = False
    dirLength = len(cwd)
    if request:
        i = 0
        openFile = open(path, 'wb')
        chunkSize = 1024 * 64
        while True:
            i += 1
            chunk = request.read(chunkSize)
            if not chunk:
                break
            openFile.write(chunk)
            outputText = 'Downloading for isotope %s and writing to %s%s' \
                         % (globalID, path[dirLength:], '.' * i)
            print(outputText, end='\r', flush=True)
        print('%s%s downloaded.\n' % (outputText, i * chunkSize), end='\r')
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
    if not rows:
        return lineListDict
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


def writeCurveToFile(curveDict, curveName, res):
    resDirectory = '%s/res%s' % (curvesDir, res)
    resFile = '%s/%s.pyr' % (resDirectory, curveName)
    openFile = open(resFile, 'ab')
    for key in curveDict:
        openFile.write(bytes('%s,' % key, 'utf-8'))
        for value in curveDict[key]:
            openFile.write(bytes('%s,' % value, 'utf-8'))
        openFile.write(bytes('\n', 'utf-8'))
    openFile.close()


def displayAllMolecules():
    newLineIter = 0
    for molecule in MOLECULE_ID.keys():
        newLineIter += 1
        print('%s\t' % molecule, end='')
        if newLineIter == 7:
            newLineIter = 0
            print('\n')

RES_MULTIPLIER = 1
BASE_RESOLUTION = .01 * RES_MULTIPLIER

MOLECULE_PARAM_COMMENTS = "#\t#\t#\n" \
                          "# Molecule params for pyrad\n" \
                          "#\t#\t#\n"


if 'win' in sys.platform.lower():
    TEXT_COLORS = {'boldMagenta': '',
                   'boldLime': '',
                   'boldBlue': '',
                   'boldCyan': '',
                   'boldWhite': '',
                   'underlineWhite': '',
                   'underlineMagenta': '',
                   'underlineLime': '',
                   'underlineCyan': '',
                   'regularMagenta': '',
                   'regularLime': '',
                   'regularCyan': '',
                   'colorEnd': ''}
else:
    TEXT_COLORS =  {'boldMagenta': '\x1b[1;31;48m',
                    'boldLime': '\x1b[1;32;48m',
                    'boldBlue': '\x1b[1;34;48m',
                    'boldCyan': '\x1b[1;36;48m',
                    'boldWhite': '\x1b[1;30;48m',
                    'underlineWhite': '\x1b[4;30;48m',
                    'underlineMagenta': '\x1b[4;31;48m',
                    'underlineLime': '\x1b[4;32;48m',
                    'underlineCyan': '\x1b[4;36;48m',
                    'regularMagenta': '\x1b[0;31;48m',
                    'regularLime': '\x1b[0;32;48m',
                    'regularCyan': '\x1b[0;36;48m',
                    'colorEnd': '\x1b[0m'}

VERSION = '1.71'
titleLine = "%s***********************              PyRad v%s              ***********************%s" \
            % (TEXT_COLORS['underlineCyan'], VERSION, TEXT_COLORS['colorEnd'])
messageGap = int((len(titleLine) - len(VERSION) - 1) / 2)
GREETING = "%s\n\n" \
           "An open-source, amateur attempt at a radiative transfer model for an atmosphere\n\n" \
           "\tAll line lists are downloaded from HITRAN\n" \
           "\tAll information for calculation of lineshapes comes \n" \
           "\tfrom spectralcalc.com and hitran.org \n" \
           "\tFor a more complete option, check out HAPI from HITRAN.\n" \
           "\tThanks to those you have made the information available \n" \
           "\tand so easily accessible.\n\n" \
           "*******************************************************************************" \
           % (titleLine)


MOLECULE_PARAM_COMMENTS = "#\t#\t#\n" \
                          "# Molecule params for pyrad\n" \
                          "#\t#\t#\n"
NULL_TAG = '#/null/#'


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
                         11: 120},
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
