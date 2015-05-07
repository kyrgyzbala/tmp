__author__ = 'Sanjarbek Hudaiberdiev'


import sys
sys.path.append('/users/hudaiber/Projects/SystemFiles/')
sys.path.append('/users/hudaiber/Projects/lib/BioPy')
import os

project_data_path = os.path.join('../', 'data', 'Archea')


def extract_class_nbr(nbr_path, arcogs_path, nbr_file, class_name):

    class_arcog_ids = open(os.path.join(arcogs_path,'arcogs_%s.txt'%class_name)).readlines()
    class_arcog_ids = [l.strip() for l in class_arcog_ids]

    picked_nbrs = []

    with open(os.path.join(nbr_path,'arcog_nbr_%s.txt'%class_name),'w') as outf:
        for l in open(nbr_file):
            terms = l.split()
            for arcog_id in class_arcog_ids:
                if arcog_id in terms:
                    outf.write(l)
                    break


def class_neighborhood_files():

    classes = ['integrase', 'recombinase', 'transposase']

    search_arcog_files =['selected_arcogs_integrase.txt',
                         'selected_arcogs_recombinase.txt',
                         'selected_arcogs_transposase.txt']

    nbr_source_file = os.path.join(project_data_path, 'neighborhoods', 'arcog_nbr.txt')
    nbr_path = os.path.join(project_data_path, 'neighborhoods')
    selected_arcogs_path = os.path.join(project_data_path, 'selected_arcogs')

    for i in range(3):
        extract_class_nbr(nbr_path, selected_arcogs_path, nbr_source_file, classes[i])


if __name__=='__main__':
    class_neighborhood_files()