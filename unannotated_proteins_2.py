__author__ = 'Sanjarbek Hudaiberdiev'

import os
from ftplib import FTP
import urllib
import sys
sys.path.append('/users/hudaiber/Projects/lib/BioPy')
from BioClasses import Gene
from copy import deepcopy

from Bio import SeqIO

def download_from_ftp():

    genomes = [l.strip() for l in open(os.path.join(folder_path, 'genome_names.txt')).readlines()]

    ftpLink = 'ftp.ncbi.nlm.nih.gov'
    ftpBacGenomes = 'http://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/'

    save_path = os.path.join(folder_path, 'genomes')

    ftp = FTP(ftpLink)
    ftp.login()
    ftp.cwd('genomes')
    ftp.cwd('Bacteria')

    for genome in genomes:

        try:
            cur_files = ftp.nlst(genome)
        except:
            print 'Couldn\'t find', genome
            continue

        cur_dir = os.path.join(save_path, genome)

        if not os.path.exists(cur_dir):
            os.mkdir(cur_dir)

        cur_files = [l.split()[-1] for l in cur_files if 'fna' in l]

        for cur_file in cur_files:
            sourceFile=ftpBacGenomes+cur_file
            targetFile=os.path.join(cur_dir, os.path.basename(cur_file))
            urllib.urlretrieve(sourceFile, targetFile)


def frame_translate(seq, frame):
    out_seq = deepcopy(seq)
    out_seq.id = 'Frame: ' + str(frame)
    out_seq.description = ""

    if frame>0:
        out_seq.seq = out_seq.seq[frame-1:].translate()
    else:
        out_seq.seq = out_seq.seq.reverse_complement()[abs(frame)-1:].translate()
    return out_seq


def translate(seq):
    frames = [1, 2, 3, -1, -2, -3]
    seqs = [frame_translate(seq, frame) for frame in frames]
    return seqs


def pty2list(pty_file_path):

    pty_list = []

    with open(pty_file_path) as f:
        for l in f:            
            gid, coordinates, strand, genome, chromosome = l.strip().split()
            pfrom, pto = coordinates.split('..')

            curGene = Gene(source=chromosome, gid=gid, pFrom=pfrom, pTo=pto, organism=genome, strand=strand)
            pty_list.append(curGene)

    return pty_list


def list2pty(pty_list, fname):

    out_fmt = "%s\t%d..%d\t%s\t%s\t%s\n"

    with open(fname, 'w') as f:
        for g in pty_list:
            out_line = out_fmt % (g.gid, g.pFrom, g.pTo, g.strand, g.organism, g.src)
            f.write(out_line)


def annotate_proteins():

    prots_file_path    = '/Users/hudaiber/data/Pty/unannotated_proteins/unannotataed_prots_2.fa'
    pty_source_path    = '/Users/hudaiber/data/Pty/genomes/'
    prots_genomes_path = '/Users/hudaiber/data/Pty/unannotated_proteins/genomes/'

    prots = SeqIO.parse(prots_file_path, 'fasta')

    for seq in prots:
        
        curGenome = seq.description.split('|')[-1]        
        search_seq = seq.seq.tostring().replace('X','*')

        print seq.id, '\t', curGenome, '\t', seq.seq
        found = False

        fna_files = os.listdir(os.path.join(prots_genomes_path,curGenome))
        fna_files = [l for l in fna_files if l.endswith('.fna')]

        for f in fna_files:

            if found:
                continue

            dna_file = os.path.join(prots_genomes_path, curGenome, f)
            dna_seq = SeqIO.read(dna_file,'fasta')
            six_frame_tr = translate(dna_seq)
            chromosome = f.split('.')[0]
            

            for frame_tr in six_frame_tr[:3]:
                pos = frame_tr.seq.find(search_seq)                
                if pos!=-1:
                    reading_frame = int(frame_tr.id.split()[1])
                    pFrom = pos*3 + reading_frame-1
                    pTo = (pos+len(seq.seq))*3 + reading_frame - 1
                    strand = '+'
                    found = True

            for frame_tr in six_frame_tr[3:]:
                pos = frame_tr.seq.find(search_seq)
                dna_len = len(dna_seq.seq)
                if pos!=-1:
                    reading_frame = int(frame_tr.id.split()[1])
                    pFrom = pos*3 + abs(reading_frame) - 1
                    pTo = (pos+len(seq.seq))*3 + abs(reading_frame) - 1
                    pFrom, pTo = dna_len - pTo, dna_len - pFrom
                    strand = '-'
                    found = True

            if found:

                cur_seq = Gene(source=chromosome, gid=seq.id.split('|')[1], pFrom=pFrom, pTo=pTo, organism=curGenome, strand=strand)                
                pty_file_source = os.path.join(pty_source_path, curGenome, '%s.pty'%chromosome)
                pty_lines = pty2list(pty_file_source)
                pty_lines.append(cur_seq)
                pty_lines.sort()
                list2pty(pty_lines, pty_file_source)

        if not found:
            print 'Not found:', seq.id



if __name__=='__main__':
    folder_path = '/Users/hudaiber/data/Pty/unannotated_proteins'
    fasta_file = os.path.join(folder_path, 'unannotataed_prots_2.fa')

    # download_from_ftp()
    annotate_proteins()
