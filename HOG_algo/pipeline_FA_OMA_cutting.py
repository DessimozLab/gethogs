__author__ = 'traincm'

import subprocess

cutting_param = None
pathfolder_output = None
pathfile_cutting = None
import os

taxons=["Homininae","Murinae"]




def main_pipeline_FA(cut):
    global cutting_param
    global pathfile_cutting
    global pathfolder_output
    cutting_param = str(cut)
    pathfile_cutting="../Result/OMA_cutting_"+cutting_param+"/OMA_HOGS_with_cutparam_"+cutting_param+"_oid.orthoxml"
    pathfolder_output="../Result/OMA_cutting_"+cutting_param+"/FA/"
    for taxon in taxons:
        if taxon =='Homininae':
            os.system("cd ../family-analyzer; bin/familyanalyzer --xreftag omaId "+ pathfile_cutting +" "+ taxon +" HUMAN PANTR GORGO >"  +pathfolder_output+ "Homininae.txt" )
            os.system("cd ../family-analyzer; bin/familyanalyzer --xreftag omaId "+ pathfile_cutting +" "+ taxon +" HUMAN PANTR >"  +pathfolder_output+ "Homininae_HP.txt" )
        elif taxon == "Murinae":
            os.system("cd ../family-analyzer; bin/familyanalyzer --xreftag omaId "+ pathfile_cutting +" "+ taxon +" MOUSE RATNO >"  +pathfolder_output+ "Murinae.txt" )




        '''
        elif taxon == "Euarchontoglires":
            write_output_FA_into_file("Euarchontoglires.txt", ["HUMAN", "PANTR", "GORGO", "MOUSE", "RATNO"], taxon)
            write_output_FA_into_file("Euarchontoglires_GHP.txt", ["HUMAN", "PANTR", "GORGO"], taxon)
            write_output_FA_into_file("Euarchontoglires_HP.txt", ["HUMAN", "PANTR"], taxon)
            write_output_FA_into_file("Euarchontoglires_MR.txt", ["MOUSE", "RATNO"], taxon)
        elif taxon == "Boreoeutheria":
            write_output_FA_into_file("Boreoeutheria.txt", ["CANFA", "HUMAN", "PANTR", "GORGO", "MOUSE", "RATNO"], taxon)
            write_output_FA_into_file("Boreoeutheria_GHP.txt", ["HUMAN", "PANTR", "GORGO"], taxon)
            write_output_FA_into_file("Boreoeutherias_HP.txt", ["HUMAN", "PANTR"], taxon)
            write_output_FA_into_file("Boreoeutherias_MR.txt", ["MOUSE", "RATNO"], taxon)

'''
