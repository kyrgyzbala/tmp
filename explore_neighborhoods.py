__author__ = 'Sanjarbek Hudaiberdiev'

import sys
sys.path.append('/users/hudaiber/Projects/SystemFiles/')
sys.path.append('/users/hudaiber/Projects/lib/BioPy')
from BioClasses import Gene
import globalVariables as gv
import os
import CogClasses as cc
import cPickle as pickle
import tools as t
from collections import Counter
import xlsxwriter as x

#Global variables

genomesDataPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Archea', 'genomes')
projectDataPath = os.path.join('../', 'data', 'Archea')
ptyGenomesPath = os.path.join(gv.LOCAL_DATA_FOLDER, 'Pty', 'genomes')

# projectDataPath = '/home/sanjarbek/NCBI_projects/NewSystems/data/Archea/'

FLANK_LENGTH = 50


def arcog_hits_to_genes():
    arcog_hit_file = os.path.join(projectDataPath, 'map_arcogs_neighborhoods', 'mapped_sorted_ar14.arCOG.csv')
    arcog_hits = []
    for l in open(arcog_hit_file):
        terms = l.split()
        arcog_hits.append(Gene(source=terms[1], gid=terms[2], pFrom=terms[4], pTo=terms[5], organism=terms[0], strand=terms[3], arcogid=terms[6]))
    return arcog_hits


def get_class_arcogs(class_name):
    fname = os.path.join(projectDataPath, 'selected_arcogs', 'arcogs_%s.txt'%class_name)
    return [l.strip() for l in open(fname).readlines()]


def get_cog_def():

    return {l.split()[0]:' '.join(l.split()[1:]) for l in open(os.path.join(projectDataPath,'additional','COGdef.txt'))}


def get_arcog_def():
    def_map = {}
    for l in open(os.path.join(projectDataPath, 'additional', 'arCOGdef.txt')):
        terms = l.split()
        k = terms[0]
        v = ' '.join(terms[3:]) if terms[2]=='-' else ' '.join(terms[2:])
        def_map[k]=v
    return def_map


def get_selected_neighborhoods():
    # 'arcog_nbr_recombinase.txt' is not included, because no neighborhoods found to have recombinase arcog ids
    class_neighborhood_files = ['arcog_nbr_integrase.txt', 'arcog_nbr_transposase.txt']
    class_neighborhood_files = [os.path.join(projectDataPath, 'neighborhoods', f) for f in class_neighborhood_files]

    neighborhoods = []

    for f in class_neighborhood_files:
        neighborhoods += tuple(sorted(l.split()[1:]) for l in open(f).readlines())
    neighborhoods = sorted(neighborhoods)
    neighborhoods = [neighborhoods[i] for i in range(len(neighborhoods)) if i==0 or neighborhoods[i]!=neighborhoods[i-1]]
    return neighborhoods


def pack_neighborhoods(neighborhood_arcogid_lists):

    """    For each list of arcogs in 'neighborhoods', form COG_neighborhood objects with:
           - Organisms that it was found in
           - For each organism; for each source, retreive arcog_hit lines as Gene objects"""

    neighborhoods=[]
    for arcogid_list in neighborhood_arcogid_lists:
        neighborhoods.append(cc.COG_neighborhood(coglist=arcogid_list))

    arcog_hits = arcog_hits_to_genes()
    organisms = set(g.organism for g in arcog_hits)

    for organism in organisms:        
        org_arcog_hits = [g for g in arcog_hits if g.organism == organism]
        source_names = set(g.src for g in org_arcog_hits)

        for source in source_names:
            source_arcog_hits = [g for g in org_arcog_hits if g.src == source]
            source_arcog_hits.sort()
            source_arcog_ids = set([g.arcogid for g in source_arcog_hits])

            for nbr in neighborhoods:
                if nbr.set_cogs < source_arcog_ids:
                    nbr.organisms.add(organism)

    selected_neighborhoods = [nbr for nbr in neighborhoods if len(nbr.organisms)>1]

    return selected_neighborhoods


def get_match_all_members(arcog_hits, arcog_ids):

    """From arcog_hits, get regions where all of the arcogs in arcog_ids exists"""
    blocks=[]
    block_size = 2*len(arcog_ids)

    for i in range(0, len(arcog_hits)-block_size+1, block_size):
        block = arcog_hits[i:i+block_size]
        block_arcogs = [g.arcogid for g in block]

        if sum([1 for a in arcog_ids if a in block_arcogs]) == len(arcog_ids):
            blocks.append(block)

    return blocks


def get_matching_blocks(arcog_hits, arcog_ids):
    matching_blocks = get_match_all_members(arcog_hits, arcog_ids)
    matching_blocks = merge_nearby_blocks(matching_blocks)
    return matching_blocks


def label_class_arcogs(blocks):
    integrase_arcogs = get_class_arcogs('integrase')
    transposase_arcogs = get_class_arcogs('transposase')
    recombinase_arcogs = get_class_arcogs('recombinase')

    for block in blocks:
        for g in block:
            if g.arcogid in recombinase_arcogs:
                g.arcogid = 'r_'+g.arcogid
            if g.arcogid in integrase_arcogs:
                g.arcogid = 'i_'+g.arcogid
            if g.arcogid in transposase_arcogs:
                g.arcogid = 't_'+g.arcogid
    return blocks


def add_flanks(nbr_blocks, ptt_path):

    out_blocks = []
    ptt_map = t.get_ptt_map(ptt_path)

    for block in nbr_blocks:

        if len(block) != len(set(g.gid for g in block)):
            new_block=[]
            block_gis = [g.gid for g in block]
            block_gis_count = Counter(block_gis)
            for gi, cnt in block_gis_count.items():
                if cnt > 1:
                    gi_seqs = [g for g in block if g.gid==gi]
                    cur_seq = gi_seqs[0]
                    cur_seq.arcogid = " ".join([g.arcogid for g in gi_seqs])
                    cur_seq.cogid = ptt_map[gi].cogid
                else:
                    cur_seq = [g for g in block if g.gid==gi][0]
                    cur_seq.cogid = ptt_map[gi].cogid
                new_block.append(cur_seq)
            block = new_block
        block.sort()
        first, last = block[0], block[-1]

        source_genes = ptt_map.values()
        source_genes.sort()

        for i in range(len(source_genes)):
            if source_genes[i].gid == first.gid:
                left_flank = source_genes[i-FLANK_LENGTH:i]
            if source_genes[i].gid == last.gid:
                right_flank = source_genes[i:i+FLANK_LENGTH]
                break
        block = left_flank+block+right_flank
        out_blocks.append(block)

    return out_blocks


def align_neighborhoods_flanking_regions(neighborhoods):

    nbrs_out_dir = os.path.join(projectDataPath, 'neighborhoods', 'neighborhoods')
    cnt = 1
    cog_def_map = get_cog_def()
    arcog_def_map = get_arcog_def()
    arcog_hits = arcog_hits_to_genes()

    nbr_to_report = []

    for nbr in neighborhoods:
        cnt+=1
        out_organisms, out_sources, out_columns, out_columns_depo = [], [], [], []
        cur_rep_nbr = cc.REP_neighborhood(nbr.cogs)

        for organism in nbr.organisms:

            cur_rep_organism = cc.REP_organism(organism)
            org_arcog_hits = [g for g in arcog_hits if g.organism==organism]
            sources = set([g.src for g in org_arcog_hits])

            for source in sources:
                source_arcog_hits = [g for g in org_arcog_hits if g.src==source]
                source_arcog_hits.sort()
                ptt_path = os.path.join(genomesDataPath, organism, '%s.ptt'%source)
                nbr_blocks = get_matching_blocks(source_arcog_hits, nbr.cogs)

                if nbr_blocks:
                    cur_rep_source = cc.REP_source(source)
                    nbr_blocks = label_class_arcogs(nbr_blocks)
                    cur_rep_source.blocks = add_flanks(nbr_blocks, ptt_path)
                    cur_rep_organism.sources.append(cur_rep_source)

            if cur_rep_organism.sources:
                cur_rep_nbr.organisms.append(cur_rep_organism)

        if cur_rep_nbr.organisms:
            nbr_to_report.append(cur_rep_nbr)
    return nbr_to_report


def write_to_file(nbr, nbr_rep_file):

    cog_def = get_cog_def()

    arcog_def = get_arcog_def()

    workbook = x.Workbook(nbr_rep_file)
    worksheet = workbook.add_worksheet()

    row_len = 5
    column_names = ['From','To','Strand','COG','arCOG']

    # Writing header section

    title_format = workbook.add_format()
    title_format.set_font_size(14)
    title_format.set_bold()
    title_format.set_align('center')

    header_format = workbook.add_format()
    header_format.set_font_size(12)
    header_format.set_bold()
    header_format.set_align('center')

    top_border = 0

    worksheet.merge_range(0, 0, 0, 10, 'Neighborhood: '+ ' '.join(nbr.cogs), title_format)
    top_border += 1

    for arcogid in nbr.cogs:
        worksheet.merge_range(top_border, 0, top_border, 10, '%s:  %s'%(arcogid, arcog_def[arcogid]))
        top_border += 1

    left_border = 0

    for org in nbr.organisms:
        cur_top_border = top_border
        org_len = (row_len+1)*sum([1 for source in org.sources for b in source.blocks])
        worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + org_len-2, org.name, header_format)
        cur_top_border += 1

        for source in org.sources:
            for block in source.blocks:
                worksheet.merge_range(cur_top_border, left_border, cur_top_border, left_border + row_len-1, source.name, header_format)
                worksheet.write_row(cur_top_border+1, left_border, column_names, header_format)
                left_border += row_len+1

    top_border += 4
    # Writing the data

    left_border = 0

    for org in nbr.organisms:
        for source in org.sources:
            for block in source.blocks:
                for i in range(len(block)):
                    cur_seq = block[i]
                    if cur_seq.cogid==None:
                        cur_seq.cogid='-'
                    cur_cogid = cur_seq.cogid if len(cur_seq.cogid)==7 else cur_seq.cogid[:-1]
                    cog_desc = cog_def[cur_cogid] if cog_def.has_key(cur_cogid) else '-'

                    arcog_desc = ' '
                    if cur_seq.arcogid:
                        if len(cur_seq.arcogid.split())>1:
                            for arcogid in cur_seq.arcogid.split():
                                arcog_desc += arcog_def[arcogid[2:]] if '_' in arcogid else arcog_def[arcogid] + ' || '
                        else:
                            arcog_desc = arcog_def[cur_seq.arcogid[2:]] if '_' in cur_seq.arcogid else arcog_def[cur_seq.arcogid]
                    
                    data_raw = [cur_seq.pFrom, cur_seq.pTo, cur_seq.strand, cog_desc, arcog_desc,' ']
                    worksheet.write_row(i+top_border, left_border, data_raw)            
                left_border += row_len+1


    # worksheet.write('B1', 'Cell B1', format)


    workbook.close()


def write_to_files(nbr_to_report, path):

    cnt=1
    ffmt='neighborhoods_%d_%d.xlsx'
    for nbr in nbr_to_report:
        nbr_rep_file = os.path.join(path, ffmt%(len(nbr.cogs),cnt))
        cnt+=1
        write_to_file(nbr, nbr_rep_file)


if __name__=='__main__':

    # selected_neighborhood_arcogid_lists = get_selected_neighborhoods()
    # selected_neighborhoods = pack_neighborhoods(selected_neighborhood_arcogid_lists)
    #
    # pickle.dump(selected_neighborhoods, open('selected_nbrhoods.p', 'w'))
    # selected_neighborhoods = pickle.load(open('selected_nbrhoods.p'))
    #
    # rep_neighborhoods = align_neighborhoods_flanking_regions(selected_neighborhoods)
    # pickle.dump(rep_neighborhoods, open('rep_neighborhoods.p','w'))
    rep_neighborhoods = pickle.load(open('rep_neighborhoods.p'))

    report_files_path = os.path.join(projectDataPath, 'report', 'xls')
    write_to_files(rep_neighborhoods, report_files_path)