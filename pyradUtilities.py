import os
import urllib.request as urlrequest

cwd = os.getcwd()
dataDir = '%s/data' % cwd
curvesDir = '%s/curves' % dataDir


def setupDir():
    if not os.path.isdir(dataDir):
        os.makedirs(dataDir)
    if not os.path.isdir(curvesDir):
        os.makedirs(curvesDir)


def openReturnLines(fullPath):
    if not os.path.isfile(fullPath):
        return []
    openFile = open(fullPath)
    lineList = openFile.readlines()
    openFile.close()
    return lineList


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


GREETING = "***********************           PyRad ver 1           ***********************\n" \
           "An open-source, amateur attempt at a radiative transfer model for an atmosphere\n\n" \
           "\tAll line lists are downloaded from HITRAN\n" \
           "\tAll information for calculation of lineshapes comes \n" \
           "\tfrom spectralcalc.com and hitran.org \n" \
           "\tFor a more complete option, check out HAPI from HITRAN.\n" \
           "\tThanks to those you have made the information available \n" \
           "\tand so easily accessible.\n\n" \
           "*******************************************************************************"
print(GREETING)
setupDir()