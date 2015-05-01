__author__ = 'Sanjarbek Hudaiberdiev'


import sys
sys.path.append('/users/hudaiber/Projects/SystemFiles/')
sys.path.append('/users/hudaiber/Projects/lib/BioPy')
import globalVariables as gv
import os
from ftplib import FTP
import urllib
from collections import deque


#Global variables

arGOGDataPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Archea', 'arCOG', 'ar14')
genomesDataPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Archea', 'genomes')
projectDataPath = os.path.join('../', 'data', 'Archea')
neigborhoodsPath = os.path.join(projectDataPath, 'neighborhoods')
targetArCOGgenesFile = os.path.join(projectDataPath,'target_arcog_genes.txt')

columnNum = 7  # arCOG ids are stored in the 7th column in ar14.arCOG.csv

flankingGeneWindow = 50

def getMatchingNeighborhoodsFromArcogArchive(className):

    arCOGarchiveFile = os.path.join(arGOGDataPath, 'ar14.arCOG.csv')
    arCOGnbrFile = os.path.join(projectDataPath, "arcog_nbr_%s.txt" % className)

    tmpf = open(os.path.join(projectDataPath, 'target_nrbs_%s.txt' % className), 'w')


    arch_lines = [l for l in open(arCOGarchiveFile).readlines() if 'arCOG' in l]

    for l in open(arCOGnbrFile).readlines():
        nbrhd = l.split()[1:]
        nbrhd_reverse = l.split()[1:]
        nbrhd.reverse()

        cnt = 0
        detected_block = []
        active_nbr = None

        for i in range(len(arch_lines)):
            l = arch_lines[i]

            if nbrhd[0] in l:
                active_nbr = nbrhd

            if nbrhd_reverse[0] in l:
                active_nbr = nbrhd_reverse

            if active_nbr:
                detected_block = arch_lines[i:i+len(active_nbr)]
                print detected_block
                detected_block = [tmpL.split(',')[6] for tmpL in detected_block]
                print detected_block
                sys.exit()
                print 'detected', detected_block
                print 'active  ', active_nbr
                if detected_block == active_nbr:
                    print detected_block
                    print active_nbr
                    sys.exit()
                active_nbr = None






    tmpf.close()


def getGenomicFiles():

    ftpLink = 'ftp.ncbi.nlm.nih.gov'
    ftpBacGenomes = 'http://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/'

    archeaDataPath = os.path.join(gv.LOCAL_DATA_FOLDER,'Archea')
    savePath=os.path.join(archeaDataPath,'genomes')

    tmpLines = open(targetArCOGgenesFile).readlines()
    gnmNames = set([l.split(',')[1] for l in tmpLines])

    ftp = FTP(ftpLink)
    ftp.login()
    ftp.cwd('genomes')
    ftp.cwd('Bacteria')
    cnt=0

    for gnmName in gnmNames:
        curDir=os.path.join(savePath, gnmName)
        if not os.path.exists(curDir):
            os.mkdir(curDir)

        curFiles = []
        try:
            curFiles = ftp.nlst(gnmName)
        except:
            print 'Couldn\'t find',gnmName
            continue

        curFiles = [l.split()[-1] for l in curFiles if 'ptt' in l]

        cnt+=1
        # print cnt,gnmName

        for curFile in curFiles:
            sourceFile=ftpBacGenomes+curFile
            targetFile=os.path.join(curDir, os.path.basename(curFile))
            urllib.urlretrieve(sourceFile, targetFile)


def writeToFile(targetPath, fname, dataLines, headerLines):

    with open(os.path.join(targetPath,'%s.nbr'%fname),'w') as f:
        [f.write(l) for l in headerLines]
        [f.write(l) for l in dataLines]

def prepareInfoHeader(gnmName, pttName, pttFirstLines):
    infoLines=[]

    infoLines.append("Genome name:%s\n"%gnmName)
    infoLines.append("ptt file:%s\n"%pttName)
    infoLines += pttFirstLines

    return infoLines


def extractNeighborhoodsClass(classFile):

    nbrDir = os.path.splitext(classFile)[0].split('_')[-1]
    tmpNeighborhoodsPath = os.path.join(neigborhoodsPath, nbrDir)

    arcogs = [l.strip() for l in open(classFile).readlines()]

    selectedLines = []
    with open(targetArCOGgenesFile) as depoFile:
        l = depoFile.readline()
        while l:
            if l.split(',')[columnNum-1] in arcogs:
                selectedLines.append(l)
            l = depoFile.readline()

    distinctGenomes = [l.split(',')[1] for l in selectedLines]
    distinctGenomes = set(distinctGenomes)

    for gnm in distinctGenomes:
        gnmLines = [l for l in selectedLines if l.split(',')[1] == gnm]

        gids = [l.split(',')[0] for l in gnmLines]

        gnmPath = os.path.join(genomesDataPath,gnm)
        pttFiles = [l for l in os.listdir(gnmPath) if l.endswith('ptt')]

        for pttFile in pttFiles:
            pttLines = open(os.path.join(gnmPath, pttFile)).readlines()
            pttHeaderLines = pttLines[:3]
            pttLines = pttLines[3:]

            leftFlank, rightFlank = [], []
            leftFlank = deque(leftFlank)

            for i in range(0, len(pttLines)-1):
                tmpLine = pttLines[i]
                tmpGid = tmpLine.split('\t')[3]

                if tmpGid in gids:
                    rightFlank = pttLines[i+1:i+flankingGeneWindow+1]
                    tmpLines = list(leftFlank) + ['**'+tmpLine] + rightFlank

                    headerLines = prepareInfoHeader(gnm, pttFile, pttHeaderLines)
                    writeToFile(tmpNeighborhoodsPath, tmpGid, tmpLines, headerLines)

                if len(leftFlank)==flankingGeneWindow:
                    leftFlank.popleft()
                leftFlank.append(tmpLine)


def extractNeighborhoods():

    arCogFiles=['selected_arcogs_integrase.txt',
                'selected_arcogs_recombinase.txt',
                'selected_arcogs_transposase.txt']

    for curArcogFile in arCogFiles:
        print curArcogFile
        curArcogFile = os.path.join(projectDataPath,curArcogFile)
        extractNeighborhoodsClass(curArcogFile)













