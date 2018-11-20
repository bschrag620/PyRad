import os
import sys
import urllib.request as urlrequest
import urllib.error as urlexception
from datetime import datetime


cwd = os.getcwd()
lineSep = os.linesep
dataDir = '%s/data' % cwd
molParamsFile = '%s/molparams.txt' % dataDir
profileDir = '/%s/profiles' % dataDir
themeDir = '%s/themes' % cwd
debuggerFilePath = '%s/logger.txt' % cwd
now = datetime.now()
debuggerFile = open(debuggerFilePath, 'wb')
debuggerFile.write(bytes('%s\n' % now.strftime("%Y-%m-%d %H:%M:%S"), 'utf-8'))
debuggerFile.close()


class Theme:
    def __init__(self, value='dark'):
        self.theme = value
        self.faceColor = None
        self.textColor = None
        self.backingColor = None
        self.gridList = []
        self.colorList = []
        self.loadTheme(self.theme)

    def rgbTuple(self, value):
        rgb = value.split(',')
        return (float(rgb[0]) / 255, float(rgb[1]) / 255, float(rgb[2]) / 255)

    def loadTheme(self, value):
        fullPath = '%s/%s.pyr' % (themeDir, value)
        lines = openReturnLines(fullPath)
        self.colorList = []
        self.gridList = []
        for line in lines:
            cells = line.split(':')
            if cells[0] == 'colorList':
                self.colorList.append(self.rgbTuple(cells[1]))
            elif cells[0] == 'textColor':
                self.textColor = self.rgbTuple(cells[1])
            elif cells[0] == 'faceColor':
                self.faceColor = self.rgbTuple(cells[1])
            elif cells[0] == 'gridList':
                self.gridList.append(self.rgbTuple(cells[1]))
            elif cells[0] == 'backingColor':
                self.backingColor = self.rgbTuple(cells[1])

    @property
    def listOfThemes(self):
        themeList = []
        fileList = getPyrFileList(themeDir)
        for file in fileList:
            themeList.append(file.split('.')[0])
        return themeList


class Settings:
    def __init__(self, setting, userName='@bradschrag'):
        self.setting = setting
        self.reduceRes = True
        self.userName = userName

    def changeSetting(self, value):
        self.setting = value

    def changeResReduce(self, bool):
        self.reduceRes = bool


    @property
    def baseRes(self):
        if self.setting == 'low':
            return .01
        elif self.setting == 'mid':
            return .01
        elif self.setting == 'hi':
            return .01

    @property
    def lineIntensityCutoff(self):
        if self.setting == 'low':
            return 1E-21
        elif self.setting == 'mid':
            return 1E-24
        elif self.setting == 'hi':
            return 0.0

    @property
    def smoothing(self):
        if self.setting == 'low':
            return 1
        elif self.setting == 'mid':
            return 50
        elif self.setting == 'hi':
            return 100


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
    directoryList = [dataDir, profileDir, themeDir]
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
    # check theme
    filePath = '%s/dark.pyr' % themeDir
    if not os.path.isfile(filePath):
        writeTheme(filePath)
    print('files checked.')
    return directoryCheck and fileCheck


def writeTheme(fullPath):
    openFile = open(fullPath, 'wb')
    text = '# created by PyRad v%s on %s.\n' \
           '# values are (r, g, b) [0-255]\n' % (VERSION, now.strftime("%Y-%m-%d %H:%M:%S"))
    properties = THEME['dark']
    value = 'faceColor'
    r,g,b = properties[value]
    text += '%s: %s, %s, %s\n' % (value, r, g, b)
    value = 'gridColor'
    r, g, b = properties[value]
    text += '%s: %s, %s, %s\n' % (value, r, g, b)
    value = 'textColor'
    r, g, b = properties[value]
    text += '%s: %s, %s, %s\n' % (value, r, g, b)
    for color in properties['colorList']:
        r,g,b = color
        text += '%s: %s, %s, %s\n' % ('colorList', r,g,b)
    openFile.write(text.encode('utf-8'))
    openFile.close()
    return


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


def writePlanetProfile(name, layer, processingTime, moleculeList, moleculeSpecific=False):
    folderPath = '%s/%s' % (profileDir, name)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
    filePath = '%s/%s.pyr' % (folderPath, layer.name)
    now = datetime.now()
    openFile = open(filePath, 'wb')
    text = '# created by PyRad v%s on %s.\n' % (VERSION, now.strftime("%Y-%m-%d %H:%M:%S"))
    openFile.write(text.encode('utf-8'))
    text = ('time: %ssecs\n'
            'depth: %s\n'
            'T: %s\n'
            'P: %s\n'
            'rangeMin: %s\n'
            'rangeMax: %s\n'
            'height: %s\n'
            'name: %s\n'
            'molecule list: %s\n'
            '# layer abs coef\n'
            'layer absCoef: %s\n'
            % (int(processingTime), layer.depth, int(layer.T), layer.P, layer.rangeMin, layer.rangeMax, layer.height,
               layer.name, ','.join(moleculeList), ','.join(map(str, layer.absCoef.tolist()))))
    openFile.write(text.encode('utf-8'))
    if moleculeSpecific:
        for molecule in layer:
            text1 = '# %s abs coef\n' % molecule.name
            text2 = '%s absCoef: %s\n' % (molecule.name, ','.join(map(str, molecule.absCoef.tolist())))
            text3 = '%s concentration: %s\n' % (molecule.name, molecule.concentration)
            openFile.write(text1.encode('utf-8'))
            openFile.write(text2.encode('utf-8'))
            openFile.write(text3.encode('utf-8'))
    openFile.close()
    print('\t\t\t\t\t\t ||> %s.pyr' % layer.name, end='\r', flush=True)
    return


def writePlanetTransmission(name, height, values, direction, number, mode='wb'):
    folderPath = '%s/%s' % (profileDir, name)
    filePath = '%s/%s.pyr' % (folderPath, 'trans looking %s-%s' % (direction, number))
    openfile = open(filePath, mode)
    if mode == 'wb':
        text = '# PyRad v%s transmission file\n' \
           '# layer height for this file is %s\n' % (VERSION, height)
        openfile.write(text.encode('utf-8'))
        text = '# general layer data is listed below this line. Other data contained within this file is\n' \
               '# the reduced transmission for the layer and each molecule that is part of the atmosphere.\n' \
               '# Effective emissivity for is the avg emissivity, normalized is the avg with zeroes removed.\n#\n'
        openfile.write(text.encode('utf-8'))
    for item in values:
        try:
            iter(values[item])
            value = ','.join(map(str, values[item]))
        except TypeError:
            value = values[item]
        text = '%s: %s\n' % (item, value)
        openfile.write(text.encode('utf-8'))
    openfile.write('#\n'.encode('utf-8'))
    openfile.close()
    return


def clearAbsData(folderPath):
    fullPath = '%s/%s' % (profileDir, folderPath)
    fileList = os.listdir(fullPath)
    for file in fileList:
        if 'layer' in file:
            filePath = '%s/%s' % (fullPath, file)
            print('removing %s' % file, end='\r', flush=True)
            os.remove(filePath)
    return


def emptyFile(filePath):
    fullPath = '%s/%s' % (cwd, filePath)
    openfile = open(fullPath, 'wb')
    openfile.close()
    return


def getPyrFileList(folder=cwd):
    fileList = os.listdir(folder)
    profileFiles = []
    for file in fileList:
        if '.pyr' in file:
            profileFiles.append(file)
    return profileFiles


def getCompletedProfileList():
    fileList = os.listdir(profileDir)
    dirList = []
    completedList = []
    for file in fileList:
        if os.path.isdir('%s/%s' % (profileDir, file)):
            dirList.append('%s/%s' % (profileDir, file))
    for dir in dirList:
        if os.path.isfile('%s/profileComplete.pyr' % dir):
            completedList.append(dir.split('/')[-1])
    return completedList


def getCompletedTransmissionList():
    fileList = os.listdir(profileDir)
    dirList = []
    completedList = []
    for file in fileList:
        if os.path.isdir('%s/%s' % (profileDir, file)):
            dirList.append('%s/%s' % (profileDir, file))
    for dir in dirList:
        if os.path.isfile('%s/transmissionComplete.pyr' % dir):
            completedList.append(dir.split('/')[-1])
    return completedList


def readCompleteProfile(folderPath):
    fullPath = '%s/%s/profileComplete.pyr' % (profileDir, folderPath)
    lines = openReturnLines(fullPath)
    values = {}
    for line in lines:
        if line[0] == '#':
            pass
        else:
            cells = line.split(':')
            values[cells[0].strip()] = cells[1].strip()
    return values


def readCompleteTransmission(folderPath):
    fullPath = '%s/%s/transmissionComplete.pyr' % (profileDir, folderPath)
    lines = openReturnLines(fullPath)
    values = {}
    for line in lines:
        if line[0] == '#':
            pass
        else:
            cells = line.split(':')
            values[cells[0].strip()] = cells[1].strip()
    return values


def readTransmissionValues(fileName, folderPath):
    fullPath = '%s/%s/%s' % (profileDir, folderPath, fileName)
    lines = openReturnLines(fullPath)
    values = {}
    for line in lines:
        if line[0] == '#':
            pass
        else:
            cells = line.split(':')
            transmission = []
            for t in cells[1].split(','):
                transmission.append(float(t))
            values[cells[0]] = transmission
    return values


def molSpecProfileList():
    completedProfileDirList = getCompletedProfileList()
    completeList = []
    for completeDir in completedProfileDirList:
        folderPath = '%s/%s' % (profileDir, completeDir)
        fullPath = '%s/profileComplete.pyr' % folderPath
        lines = openReturnLines(fullPath)
        for line in lines:
            if line[0] == '#':
                pass
            elif line.split(':')[0] == 'molSpecific':
                if 'true' in line.split(':')[1].lower():
                    completeList.append(completeDir)
    return completeList


def checkPlanetProfile(name):
    folderPath = '%s/%s' % (profileDir, name)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
        return False
    else:
        fileName = 'profileComplete.pyr'
        filePath = '%s/%s' % (folderPath, fileName)
        if os.path.isfile(filePath):
            return True
        else:
            return False


def checkPlanetTransmission(name):
    folderPath = '%s/%s' % (profileDir, name)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
        return False
    else:
        fileName = 'transmissionComplete.pyr'
        filePath = '%s/%s' % (folderPath, fileName)
        if os.path.isfile(filePath):
            return True
        else:
            return False


def parseCustomProfile(name):
    filePath = '%s/%s.pyr' % (cwd, name)
    lines = openReturnLines(filePath)
    params = {'compositionRules': [],
              'temperatureRules': [],
              'molecules': [],
              'name': name,
              'molSpecific': False}
    for line in lines:
        if line[0] == '#' or line == '':
                pass
        elif 'composition' in line.lower() or 'lapse' in line.lower():
            tempDict = {}
            values = line.split(':')[1].split(',')
            tempDict['name'] = values[0].strip()
            tempDict['finalHeight'] = float(values[1])
            tempDict['finalValue'] = float(values[2])
            if len(values) > 3:
                tempDict['moleculeName'] = values[3].strip()
            if 'composition' in line.lower():
                params['compositionRules'].append(tempDict)
            elif 'lapse' in line.lower():
                params['temperatureRules'].append(tempDict)
        elif 'surface' in line.lower():
            tempDict = {'rangeMin': 0,
                        'rangeMax': 2000}
            values = line.split(':')[1].split(',')
            tempDict['surfacePressure'] = float(values[0])
            tempDict['surfaceTemperature'] = int(values[1])
            tempDict['maxHeight'] = float(values[2])
            tempDict['initialDepth'] = float(values[3])
            tempDict['gravity'] = float(values[4])
            if len(values) > 5:
                tempDict['rangeMin'] = int(values[5])
            if len(values) > 6:
                tempDict['rangeMax'] = int(values[6])
            params['initialValues'] = tempDict
        elif 'molecule' in line.lower():
            values = line.split(':')[1].split(',')
            try:
                concentration = float(values[1])
                tempDict = {'name': values[0].lower().strip(),
                            'concentration': concentration}
            except ValueError:
                print('invalid concentration %s in %s' % (values[1], filePath))
                exit()
            params['molecules'].append(tempDict)
        elif 'molspecific' in line.lower():
            if line.split(':')[1].lower().strip() == 'true':
                value = True
            else:
                value = False
            params['molSpecific'] = value
    return params


def profileLength(name):
    folderPath = '%s/%s' % (profileDir, name)
    filePath = '%s/profileComplete.pyr' % folderPath
    lines = openReturnLines(filePath)
    for line in lines:
        if line[0] == '#':
            pass
        else:
            values = line.split(':')
            if values[0].lower().strip() == 'completed':
                return int(values[1])


def profileWriteProgress(name, completed, expected, processingTime, moleculeSpecific, moleculeList, res):
    folderPath = '%s/%s' % (profileDir, name)
    fileName = 'profileProgress'
    filePath = '%s/%s.pyr' % (folderPath, fileName)
    text = '# profile for %s\n' \
           '# PyRad v%s\n' \
           '# total processing time: %ssecs\n' \
           'molSpecific: %s\n' \
           'molList: %s\n' \
           'res: %s\n' \
           'expected: %s\n' \
           'completed: %s' % (name, VERSION, int(processingTime),
                              moleculeSpecific, ','.join(moleculeList), res, expected, completed)
    openFile = open(filePath, 'wb')
    openFile.write(text.encode('utf-8'))
    openFile.close()
    return


def profileProgress(name):
    folderPath = '%s/%s' % (profileDir, name)
    fileName = 'profileProgress'
    filePath = '%s/%s.pyr' % (folderPath, fileName)
    startTime = 0
    if not os.path.isfile(filePath):
        return False, startTime
    lines = openReturnLines(filePath)
    for line in lines:
        if line[0] == '#':
            pass
        else:
            cells = line.split(':')
            if cells[0] == 'completed':
                completed = int(cells[1])
            elif cells[0] == 'time':
                startTime = int(cells[1].strip())
    return completed, startTime


def emptyProfileDirectory(name):
    folderPath = '%s/%s' % (profileDir, name)
    if not os.path.isdir(folderPath):
        return
    fileList = os.listdir(folderPath)
    for file in fileList:
        filePath = '%s/%s' % (folderPath, file)
        os.remove(filePath)
    os.removedirs(folderPath)
    return


def profileComplete(name):
    folderPath = '%s/%s' % (profileDir, name)
    fileName = 'profileComplete'
    filePath = '%s/%s.pyr' % (folderPath, fileName)
    if os.path.isfile(filePath):
        return True
    return False


def profileWriteTransmissionComplete(folderPath, heightList):
    folderPath = '%s/%s' % (profileDir, folderPath)
    fileName = 'profileComplete'
    filePath = '%s/%s.pyr' % (folderPath, fileName)
    completeProfileData = openReturnLines(filePath)
    transmissionPath = '%s/transmissionComplete.pyr' % folderPath
    openFile = open(transmissionPath, 'wb')
    for line in completeProfileData:
        openFile.write(line.encode('utf-8'))
    text = '\nheightList: %s' % (','.join(str(n) for n in heightList))
    openFile.write(text.encode('utf-8'))
    openFile.close()
    return


def profileWriteComplete(planet, completed, expected, processingTime, res, moleculeSpecific=False):
    folderPath = '%s/%s' % (profileDir, planet.folderPath)
    fileName = 'profileComplete'
    filePath = '%s/%s.pyr' % (folderPath, fileName)
    time = datetime.now()
    text = '# profile for %s completed on %s\n' \
           '# PyRad v%s\n' \
           '# total processing time: %ssecs\n' \
           'name: %s\n' \
           'molSpecific: %s\n' \
           'surfacePressure: %s\n' \
           'surfaceTemperature: %s\n' \
           'maxHeight: %s\n' \
           'initialDepth: %s\n' \
           'gravity: %s\n' \
           'rangeMin: %s\n' \
           'rangeMax: %s\n' \
           'molList: %s\n' \
           'res: %s\n' \
           'expected: %s\n' \
           'completed: %s' % (planet.folderPath, time.strftime("%Y-%m-%d %H:%M:%S"), VERSION, int(processingTime),
                              planet.folderPath, moleculeSpecific, planet.surfacePressure, planet.surfaceTemperature,
                              int(planet.maxHeight / 100000), planet.heightList[1], planet.gravity,
                              planet.rangeMin, planet.rangeMax, ','.join(planet.moleculeList), res, expected, completed)
    openFile = open(filePath, 'wb')
    openFile.write(text.encode('utf-8'))
    openFile.close()
    fileName = 'profileProgress'
    filePath = '%s/%s.pyr' % (folderPath, fileName)
    os.remove(filePath)
    return


def readPlanetProfileMolecule(name, layerNumber, length, moleculeName):
    folderPath = '%s/%s' % (profileDir, name)
    fileName = 'Layer %s:%s' % (layerNumber, length)
    filePath = '%s/%s.pyr' % (folderPath, fileName)
    layerDict = {'molecule list': []}
    lines = openReturnLines(filePath)
    if not lines:
        fileName = 'layer %s_%s' % (layerNumber, length)
        filePath = '%s/%s.pyr' % (folderPath, fileName)
        lines = openReturnLines(filePath)
        if not lines:
            print('cant find layer file %s' % filePath)
            exit()
    print('Reading %s profile from %s...                                ' % (moleculeName, fileName), end='\r',
          flush=True)
    for line in lines:
        keyValue = line.split(':')
        if line[0] == '#':
            pass
        elif keyValue[0].strip() == 'name':
            layerDict[keyValue[0]] = keyValue[1].strip()
        elif keyValue[0].strip() == '%s absCoef' % moleculeName:
            absList = keyValue[1].split(',')
            absCoefList = []
            for value in absList:
                absCoefList.append(float(value))
            layerDict['absCoef'] = absCoefList
        elif keyValue[0].strip() == 'T' or keyValue[0] == 'rangeMin' or keyValue[0] == 'rangeMax':
            layerDict[keyValue[0]] = int(keyValue[1])
        elif keyValue[0].strip() == 'P' or keyValue[0].strip() == 'height' or keyValue[0].strip() == 'depth' \
                or keyValue[0].strip() == '%s concentration' % moleculeName:
            layerDict[keyValue[0]] = float(keyValue[1])
        else:
            pass
    return layerDict


def readPlanetProfile(name, layerNumber, length):
    folderPath = '%s/%s' % (profileDir, name)
    fileName = 'Layer %s:%s' % (layerNumber, length)
    filePath = '%s/%s.pyr' % (folderPath, fileName)
    lines = openReturnLines(filePath)
    if not lines:
        fileName = 'layer %s_%s' % (layerNumber, length)
        filePath = '%s/%s.pyr' % (folderPath, fileName)
        lines = openReturnLines(filePath)
        if not lines:
            print('cant find layer file %s' % filePath)
            exit()
    print('Reading profile from %s...                                ' % fileName, end='\r', flush=True)
    layerDict = {'molecule list': []}
    for line in lines:
        keyValue = line.split(':')
        if line[0] == '#':
            pass
        elif keyValue[0].strip() == 'name':
            layerDict[keyValue[0]] = keyValue[1].strip()
        elif keyValue[0].strip() == 'layer absCoef':
            absList = keyValue[1].split(',')
            absCoefList = []
            for value in absList:
                absCoefList.append(float(value))
            layerDict['absCoef'] = absCoefList
        elif keyValue[0].strip() == 'T' or keyValue[0] == 'rangeMin' or keyValue[0] == 'rangeMax':
            layerDict[keyValue[0]] = int(keyValue[1])
        elif keyValue[0].strip() == 'P' or keyValue[0].strip() == 'height' or keyValue[0].strip() == 'depth' \
                or 'concentration' in keyValue[0]:
            layerDict[keyValue[0]] = float(keyValue[1])
        else:
            pass
    return layerDict


def writeDictListToFile(dictionary, fullPath, comments=None, mode='wb'):
    openFile = open(fullPath, mode)
    if comments:
        openFile.write(comments.encode('utf-8'))
    for key in dictionary:
        text = '%s,%s%s' % (key, ','.join(str(item) for item in dictionary[key]), lineSep)
        openFile.write(text.encode('utf-8'))
    openFile.close()


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
        dirLength = len(cwd)
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
    except urlexception.HTTPError:
        print('Can not retrieve data for given parameters.')
        openFile = open(path, 'wb')
        openFile.write(bytes('%s no data available for %s, range %s-%s' %
                             (NULL_TAG, globalID, waveMin, waveMax), 'utf-8'))
        openFile.close()
    except urlexception.URLError:
        print('Can not connect to %s' % str(url))


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

VERSION = '3.0'
titleLine = "%s******************              PyRad v%s              ******************%s" \
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

THEME = {'dark':    {'faceColor':   (29, 29, 29),
                     'colorList':    [(239, 244, 255),
                                      (0, 163, 0),
                                      (126,56, 120),
                                      (0, 171, 169),
                                      (185, 29,171),
                                      (30, 113, 69),
                                      (227, 162, 26)],
                     'textColor':      (239, 244, 255),
                     'gridColor':   (204, 204, 204)}}

print(GREETING)
setupDir()
