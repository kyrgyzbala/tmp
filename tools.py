__author__ = 'hudaiber'

import os
import sys
sys.path.append('/users/hudaiber/Projects/lib/BioPy')
from BioClasses import Gene

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


def extract_pty_file(genomes_path, pty_file_path):

    with open(pty_file_path) as f:
        cur_genome, cur_genome_path = '', ''
        cur_lines = []
        cur_chromosome = ''

        for l in f:
            if cur_genome == '' and cur_chromosome == '':
                cur_genome = l.split()[3]
                cur_chromosome = l.split()[4]
                cur_genome_path = os.path.join(ptyGenomesPath, cur_genome)
                if not os.path.exists(cur_genome_path):
                    os.mkdir(cur_genome_path)

            if cur_chromosome != l.split()[4]:
                cur_chromosome_file = open(os.path.join(cur_genome_path, "%s.pty" % cur_chromosome), 'w')
                [cur_chromosome_file.write(t) for t in cur_lines]
                cur_lines = []
                cur_chromosome = l.split()[4]
                cur_chromosome_file.close()

            cur_lines.append(l)
            if cur_genome != l.split()[3]:
                cur_genome = l.split()[3]
                cur_genome_path = os.path.join(ptyGenomesPath, cur_genome)
                if not os.path.exists(cur_genome_path):
                    os.mkdir(cur_genome_path)


def get_ptt_map(file_path):
    ptt_map = {}
    gnm_name = os.path.dirname(file_path).split('/')[-1]
    chromosome = os.path.splitext(os.path.basename(file_path))[0]
    lines = open(file_path).readlines()[3:]
    for l in lines:
        terms = l.strip().split()
        gid = terms[3]
        coordinates = terms[0]
        strand = terms[1]
        pfrom, pto = coordinates.split('..')
        cogid = terms[7]
        curGene = Gene(source=chromosome, gid=gid, pFrom=pfrom, pTo=pto, organism=gnm_name, strand=strand, cogid=cogid)
        ptt_map[gid] = curGene
    return ptt_map


def get_pty_map(folder_path):
    pty_map = {}
    for f in os.listdir(folder_path):
        if f.endswith('.pty'):
            for l in open(os.path.join(folder_path, f)):
                gid, coordinates, strand, genome, chromosome = l.strip().split()
                pfrom, pto = coordinates.split('..')
                curGene = Gene(source=chromosome, gid=gid, pFrom=pfrom, pTo=pto, organism=genome, strand=strand)
                pty_map[gid] = curGene
    return pty_map