__author__ = 'Sanjarbek Hudaiberdiev'


import sys
sys.path.append('/users/hudaiber/Projects/SystemFiles/')
sys.path.append('/users/hudaiber/Projects/lib/BioPy')

from BioClasses import Gene

import globalVariables as gv
import os
from CogClasses import *
import cPickle as pickle


#Global variables

arGOGDataPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Archea', 'arCOG')
genomesDataPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Archea', 'genomes')
projectDataPath = os.path.join('../', 'data', 'Archea')
ptyFilePath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Pty', 'Prok1402_3added.pty')
ptyGenomesPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Pty', 'genomes')

FLANK_LENGTH = 50


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


def prepare_line(l):
    
    out_line_fmt = "%s\t%s\t%s\t%s\t%d\t%d\t%s\n"
    active_genomes_path = ptyGenomesPath

    terms = l.split(',')
    cur_gid = terms[0]
    cur_genome = terms[1]
    cur_arcog = terms[6]

    cur_genome_path = os.path.join(active_genomes_path, cur_genome)
    pty_map = get_pty_map(cur_genome_path)
    
    gene = pty_map[cur_gid]
    out_line = out_line_fmt % (cur_genome, gene.src, cur_gid, gene.strand, gene.pFrom, gene.pTo, cur_arcog)

    return out_line


def map_all_arcog_hits():
    """arcog_hits => ar14.arCOG.csv.
    Add coordinate information to the arcog_hits"""

    arcog_hit_file_path = os.path.join(arGOGDataPath, 'ar14.arCOG.csv')
    new_arcog_hit_file_path = os.path.join(projectDataPath, 'map_arcogs_neighborhoods', 'mapped_ar14.arCOG.csv')

    cnt=0
    with open(arcog_hit_file_path) as infile:
        with open(new_arcog_hit_file_path,'w') as outfile:
            for l in infile:
                if 'arCOG' not in l:
                    continue

                try:
                    new_line = prepare_line(l)
                    outfile.write(new_line)
                except:
                    print 'Couldn\'t find:', l.strip()

                cnt+=1

                if cnt % 1000 == 0 and cnt>0:
                    print cnt


if __name__=='__main__':

    pass
    #
    #
    # integrase_arcogs   = [l.strip() for l in open(os.path.join(projectDataPath, 'selected_arcogs', 'arcogs_integrase.txt'  )).readlines()]
    # recombinase_arcogs = [l.strip() for l in open(os.path.join(projectDataPath, 'selected_arcogs', 'arcogs_recombinase.txt')).readlines()]
    # transposase_arcogs = [l.strip() for l in open(os.path.join(projectDataPath, 'selected_arcogs', 'arcogs_transposase.txt')).readlines()]
    #
    # integrases, transposases, recombinases = 0, 0, 0
    #
    # for arcog_hit in selected_arcog_hits:
    #     if arcog_hit.arcogid in integrase_arcogs:
    #         arcog_hit.set_product_enzyme('integrase')
    #         integrases += 1
    #
    #     if arcog_hit.arcogid in transposase_arcogs:
    #         arcog_hit.set_product_enzyme('transposase')
    #         transposases += 1
    #
    #     if arcog_hit.arcogid in recombinase_arcogs:
    #         arcog_hit.set_product_enzyme('recombinase')
    #         recombinases += 1
    #
    #
    # pickle.dump(selected_neighborhoods, open('selected_nbrhoods.p', 'w'))
    # pickle.dump(selected_arcog_hits, open('selected_arcog_hits.p', 'w'))

    # selected_arcog_hits = pickle.load(open('selected_arcog_hits.p'))
    # selected_neighborhoods = pickle.load(open('selected_nbrhoods.p'))

    # orgname = selected_arcog_hits[0].organism
    # org_hits = [g for g in selected_arcog_hits if g.organism==orgname]
    # org_hits.sort()
    # for g in org_hits[:2]:
    #     print g
    # sys.exit()

    # for n in selected_neighborhoods:
    #     print n.classname,  n.cogs
    #     print n.organisms

    # process_neighborhoods_of_arcog_hits(selected_neighborhoods, selected_arcog_hits)








