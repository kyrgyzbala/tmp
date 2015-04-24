__author__ = 'Sanjarbek Hudaiberdiev'


import sys
sys.path.append('/users/hudaiber/Projects/SystemFiles/')
sys.path.append('/users/hudaiber/Projects/lib/BioPy')
import globalVariables as gv
import os


if __name__=='__main__':
    selected_files = ['selected_arcogs_integrase.txt',
                      'selected_arcogs_recombinase.txt',
                      'selected_arcogs_transposase.txt']

    nbr_files = ['arcog_nbr_integrase.txt',
                   'arcog_nbr_recombinase.txt',
                   'arcog_nbr_transposase.txt']

    projectDataPath = os.path.join('../', 'data', 'Archea')
    arcog_nbr_total = os.path.join(projectDataPath,'arcog_nbr.txt')

    for i in range(3):
        tmp_selected_file = os.path.join(projectDataPath,selected_files[i])
        tmp_nbr_file = os.path.join(projectDataPath,nbr_files[i])

        with open(tmp_nbr_file,'w') as f:
            selected_arcogs = open(tmp_selected_file).readlines()
            selected_arcogs = [l.strip() for l in selected_arcogs]

            for l in open(arcog_nbr_total).readlines():
                terms = l.split()[1:]
                for term in terms:
                    if term in selected_arcogs:
                        f.write(l)
                        print term
                        continue