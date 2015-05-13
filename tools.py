__author__ = 'hudaiber'

import os
import sys
sys.path.append('/users/hudaiber/Projects/lib/BioPy')
from BioClasses import Gene
from CogClasses import Block

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


def pty_to_genes_map(folder_path, source):
    pty_map = {}

    for l in open(os.path.join(folder_path, '%s.pty'%source)):
        gid, coordinates, strand, genome, chromosome = l.strip().split()
        pfrom, pto = coordinates.split('..')
        curGene = Gene(source=chromosome, gid=gid, pFrom=pfrom, pTo=pto, organism=genome, strand=strand)
        pty_map[gid] = curGene
    return pty_map


def ptt_to_genes_map(folder_path, source):

    ptt_map = {}
    gnm_name = os.path.dirname(folder_path).split('/')[-1]
    file_path = os.path.join(folder_path, '%s.ptt'%source)
    lines = open(file_path).readlines()[3:]
    for l in lines:
        terms = l.strip().split()
        gid = terms[3]
        coordinates = terms[0]
        strand = terms[1]
        pfrom, pto = coordinates.split('..')
        cogid = terms[7]
        curGene = Gene(source=source, gid=gid, pFrom=pfrom, pTo=pto, organism=gnm_name, strand=strand, cogid=cogid)
        ptt_map[gid] = curGene
    return ptt_map


def get_ptt_map(ptt_path, pty_path, organism, source):

    """Somewhat deceiving function name. It combines ptt and pty into one map"""

    ptt_path = os.path.join(ptt_path, organism)
    pty_path = os.path.join(pty_path, organism)
    out_map={}
    ptt_map, pty_map = {}, {}

    if not os.path.exists(ptt_path):
        print 'Directory not found:', ptt_path
    if not os.path.exists(pty_path):
        print 'Directory not found:', pty_path

    try:
        ptt_map = ptt_to_genes_map(ptt_path, source)
    except Exception as e:
        print 'Could not load ptt for:', ptt_path, source
        print 'Exception:', e
        print

    try:
        pty_map = pty_to_genes_map(pty_path, source)
    except Exception as e:
        print 'Could not load pty for:', pty_path, source
        print 'Exception:', e
        print

    gids = set(ptt_map.keys()+pty_map.keys())
    for gid in gids:
        out_map[gid] = ptt_map[gid] if ptt_map.has_key(gid) else pty_map[gid]

    return out_map


def get_class_arcogs(class_name):
    fname = os.path.join(projectDataPath, 'selected_arcogs', 'arcogs_%s.txt'%class_name)
    return [l.strip() for l in open(fname).readlines()]


def get_cog_def(fpath):

    return {l.split()[0]:' '.join(l.split()[1:]) for l in fpath}


def get_arcog_def(fpath):
    def_map = {}
    for l in fpath:
        terms = l.split()
        k = terms[0]
        v = ' '.join(terms[3:]) if terms[2]=='-' else ' '.join(terms[2:])
        def_map[k]=v
    return def_map


def merge_to_regions(gene_list):
    out_list = []
    for gid in set([g.gid for g in gene_list]):
        cur_seqs = [g for g in gene_list if g.gid == gid]
        cur_seq = cur_seqs[0]
        if len(cur_seqs) > 1:
            cur_arcog_ids = set([g.arcogid for g in cur_seqs])
            if len(cur_arcog_ids) > 1:
                cur_seq.arcogid = ' '.join(cur_arcog_ids)
        out_list.append(cur_seq)

    return sorted(out_list)


def merge_to_blocks(regions, arcog_hits, max_merge_length):

    raw_blocks = []
    cur_block=[]

    for i in range(len(regions)):
        cur_gene = regions[i]

        if not cur_block:
            cur_block.append(cur_gene)
            continue

        if abs(cur_block[-1].pTo - cur_gene.pFrom) > 1000:
            raw_blocks.append(cur_block)
            cur_block = [cur_gene]
        else:
            cur_block.append(cur_gene)

    # If the two raw_blocks are within 100 genes, by their further genes, then merge these two into one,
    # filling the gap between.

    gids = [g.gid for g in arcog_hits]

    merged_blocks = []
    cur_block = []

    for i in range(len(raw_blocks)):
        cur_raw_block = raw_blocks[i]
        if not cur_block:
            cur_block = cur_raw_block
            continue

        first_outer = gids.index(cur_block[0].gid)
        last_outer = gids.index(cur_raw_block[-1].gid)

        if last_outer - first_outer < max_merge_length:
            cur_block += cur_raw_block
        else:
            merged_blocks.append(cur_block)
            cur_block = cur_raw_block

        if i == len(raw_blocks)-1:
            merged_blocks.append(cur_block)

    return merged_blocks
