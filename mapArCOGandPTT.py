__author__ = 'Sanjarbek Hudaiberdiev'


import sys
sys.path.append('/users/hudaiber/Projects/SystemFiles/')
sys.path.append('/users/hudaiber/Projects/lib/BioPy')

from BioClasses import Gene

import globalVariables as gv
import os


#Global variables

arGOGDataPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Archea', 'arCOG')
genomesDataPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Archea', 'genomes')
projectDataPath = os.path.join('../', 'data', 'Archea')
ptyFilePath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Pty', 'Prok1402_3added.pty')
ptyGenomesPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Pty', 'genomes')

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

                curGene = Gene(source=chromosome, gid=gid, pFrom=pfrom, pTo=pto, organism=gnm_name, strand=strand)
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



def get_neighborhood_related_arcog_hits():

    """arcog_hits => ar14.arCOG.csv.
    This function gets only these hits from the big archive (indicated above),
    which are hit by arcogs included in the neighborhoods"""

    class_neighborhood_files = ['arcog_nbr_integrase.txt','arcog_nbr_recombinase.txt','arcog_nbr_transposase.txt']
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
                        print 'Couldn\'t find:',l.strip()
                    outfile.write(new_line)
                    cnt+=1

                if cnt % 1000 == 0 and cnt>0:
                    print cnt

if __name__=='__main__':

    get_neighborhood_related_arcog_hits()
    # extract_pty_file(ptyGenomesPath, ptyFilePath)

    # map_arcog_with_pty(arcogfile, mapped_arCog_file)


