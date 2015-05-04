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


def get_ptt_map(folder_path):

    ptt_map = {}
    gnm_name = folder_path.split('/')[-1]
    for f in os.listdir(folder_path):
        if f.endswith('.ptt'):
            chromosome = os.path.splitext(f)[0]
            lines = open(os.path.join(folder_path, f)).readlines()[3:]
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

def map_arcog_with_pty(arcog_file_path, mapped_arcog_file):

    out_line_fmt = "%s\t%s\t%s\t%s\t%d\t%d\t%s\n"
    active_genomes_path = ptyGenomesPath

    with open(arcog_file_path) as f:
        cur_genome, cur_genome_path = '', ''

        with open(mapped_arcog_file,'w') as outf:

            for l in f:
                if cur_genome == '':
                    cur_genome = l.split(',')[0]
                    cur_genome_path = os.path.join(active_genomes_path, cur_genome)
                    ptt_map = get_pty_map(cur_genome_path)

                if cur_genome != l.split(',')[0]:

                    genes = [gene for gene in ptt_map.values() if hasattr(gene, 'arcogid')]
                    genes.sort()

                    for gene in genes:
                        out_line = out_line_fmt % (gene.organism, gene.src, gene.gid, gene.strand, gene.pFrom, gene.pTo, gene.arcogid)
                        outf.write(out_line)

                    cur_genome = l.split(',')[0]
                    cur_genome_path = os.path.join(active_genomes_path, cur_genome)
                    ptt_map = get_pty_map(cur_genome_path)

                terms = l.split(',')
                curGid, curArcog = terms[1], terms[5]

                if ptt_map.has_key(curGid):
                    ptt_map[curGid].set_arcog(curArcog)
                else:
                    terms = l.split(',')
                    curOrganism = '--'+terms[0]
                    curGene = Gene(source='none', gid=curGid, pFrom=0, pTo=0, organism=curOrganism, strand='none')
                    curGene.set_arcog(curArcog)
                    ptt_map[curGid] = curGene

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

        

def map_neighborhood_related_arcog_hits():

    """arcog_hits => ar14.arCOG.csv.
    This function gets only these hits from the big archive (indicated above),
    which are hit by arcogs included in the neighborhoods"""

    class_neighborhood_files = ['arcog_nbr_integrase.txt', 'arcog_nbr_recombinase.txt', 'arcog_nbr_transposase.txt']
    class_neighborhood_files = [os.path.join(projectDataPath, 'neighborhoods', f) for f in class_neighborhood_files]
    
    arcogs = []
    for fname in class_neighborhood_files:
        for l in open(fname).readlines():
            arcogs += l.split()[1:]

    arcogs = set(arcogs)
    
    arcog_hit_file_path = os.path.join(arGOGDataPath, 'ar14.arCOG.csv')
    new_arcog_hit_file_path = os.path.join(projectDataPath, 'map_arcogs_neighborhoods', 'filtered_mapped_ar14.arCOG.csv')

    cnt=0
    with open(arcog_hit_file_path) as infile:
        with open(new_arcog_hit_file_path,'w') as outfile:
            for l in infile:
                if not 'arCOG' in l:
                    continue

                if l.split(',')[6] in arcogs:

                    try:
                        new_line = prepare_line(l)
                    except:
                        print 'Couldn\'t find:', l.strip()
                    outfile.write(new_line)
                    cnt+=1

                if cnt % 1000 == 0 and cnt>0:
                    print cnt

def map_all_arcog_hits():
    """arcog_hits => ar14.arCOG.csv.
    Add coordinate information to the arcog_hits"""

    arcog_hit_file_path = os.path.join(arGOGDataPath, 'ar14.arCOG.csv')
    new_arcog_hit_file_path = os.path.join(projectDataPath, 'map_arcogs_neighborhoods', 'mapped_ar14.arCOG.csv')

    cnt=0
    with open(arcog_hit_file_path) as infile:
        with open(new_arcog_hit_file_path,'w') as outfile:
            for l in infile:
                if not 'arCOG' in l:
                    continue

                try:
                    new_line = prepare_line(l)
                    outfile.write(new_line)
                except:
                    print 'Couldn\'t find:', l.strip()

                cnt+=1

                if cnt % 1000 == 0 and cnt>0:
                    print cnt


def get_genomes_of_neighborhoods():

    class_neighborhood_files = ['arcog_nbr_integrase.txt', 'arcog_nbr_recombinase.txt', 'arcog_nbr_transposase.txt']
    class_neighborhood_files = [os.path.join(projectDataPath, 'neighborhoods', f) for f in class_neighborhood_files]

    arcog_hit_file = os.path.join(projectDataPath, 'map_arcogs_neighborhoods', 'mapped_sorted_ar14.arCOG.csv')

    neighborhoods = []
    arcog_hits = []

    for f in class_neighborhood_files:
        clsname = os.path.basename(f).split('.')[0].split('_')[-1]
        neighborhoods += [COG_neighborhood(l.split()[1:], clsname=clsname) for l in open(f).readlines()]

    for l in open(arcog_hit_file):
        terms = l.split()
        arcog_hits.append(Gene(source=terms[1], gid=terms[2], pFrom=terms[4], pTo=terms[5], organism=terms[0], strand=terms[3], arcogid=terms[6]))

    organisms = set(g.organism for g in arcog_hits)

    for organism in organisms:
        org_arcog_hits = [g for g in arcog_hits if g.organism==organism]
        sources = set(g.src for g in org_arcog_hits)

        for source in sources:
            source_arcog_hits = [g for g in org_arcog_hits if g.src == source]
            source_arcog_hits.sort()
            source_arcog_ids = [g.arcogid for g in source_arcog_hits]
            
            for nbr in neighborhoods:
                if nbr.set_cogs.issubset(source_arcog_ids):
                    if organism not in nbr.organisms:
                        nbr.organisms.append(organism)

    selected_neighborhoods = [nbr for nbr in neighborhoods if len(nbr.organisms)>1]

    selected_organisms = []
    for nbr in selected_neighborhoods:
        selected_organisms += nbr.organisms

    selected_arcog_hits = []
    for organism in set(selected_organisms):
        selected_arcog_hits += [g for g in arcog_hits if g.organism==organism]

    return selected_neighborhoods, selected_arcog_hits

def process_neighborhoods_of_arcog_hits(neighborhoods, arcog_hits):

    classes = ['integrase', 'transposase', 'recombinase']

    nbrhds_out_dir = os.path.join(projectDataPath, 'neighborhoods', 'neighborhoods')

    cnt = 1
    block_fmt = '%s,%s,%s,%s'
    for nbr in neighborhoods:

        out_columns = []

        for organism in nbr.organisms:
            print 'Organism name', organism

            ptt_path = os.path.join(genomesDataPath, organism)

            # ptt_map = get_ptt_map(ptt_path)
            #
            # for k,v in ptt_map.items():
            #     if not k in org_arcog_hit_gis:
            #         org_arcog_hits.append(v)


            org_arcog_hits = [g for g in arcog_hits if g.organism==organism]
            sources = set([g.src for g in org_arcog_hits])

            for source in sources:
                source_arcog_hits = [g for g in org_arcog_hits if g.src==source]
                source_arcog_hits.sort()

                left_flank, nbr_block, right_flank = [], [], []

                nbr_size = len(nbr.cogs)
                tmp_cnt=0

                for i in range(len(source_arcog_hits)):
                    if source_arcog_hits[i].arcogid in nbr.cogs:
                        nbr_block.append(source_arcog_hits[i])

                if len(nbr_block)<len(nbr.cogs):
                    nbr_block = []


                if nbr_block:
                    print source
                    for g in nbr_block:
                        print g, g.arcogid
                    print

            # print 'Nbr_block size:', len(nbr_block)
            # for g in nbr_block:
            #     print g
            # print
            org_block = []
        print

        nbr_file_name = 'nbr_%d_%d.csv' % (len(nbr.cogs), cnt)



if __name__=='__main__':

    # get_neighborhood_related_arcog_hits()
    # map_all_arcog_hits()
    # extract_pty_file(ptyGenomesPath, ptyFilePath)

    # map_arcog_with_pty(arcogfile, mapped_arCog_file)
    # selected_neighborhoods, selected_arcog_hits = get_genomes_of_neighborhoods()
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

    selected_arcog_hits = pickle.load(open('selected_arcog_hits.p'))
    selected_neighborhoods = pickle.load(open('selected_nbrhoods.p'))

    # orgname = selected_arcog_hits[0].organism
    # org_hits = [g for g in selected_arcog_hits if g.organism==orgname]
    # org_hits.sort()
    # for g in org_hits[:2]:
    #     print g
    # sys.exit()

    # for n in selected_neighborhoods:
    #     print n.classname,  n.cogs
    #     print n.organisms

    process_neighborhoods_of_arcog_hits(selected_neighborhoods, selected_arcog_hits)








