__author__ = 'traincm'

import subprocess
import classes as classes
import os

taxons=["Homininae","Murinae"]

def main_pipeline_FA(settings):
    for taxon in taxons:
        if taxon =='Homininae':
            os.system("cd ../family-analyzer; ./bin/familyanalyzer --xreftag protId "+ "../Result/"+settings.dir_name_param+"/"+ settings.xml_name_param +" "+ "HUMANPANTR" +" HUMAN PANTR  >"  +"../Result/"+ settings.dir_name_param+"/FA/"+ "Homininae.txt" )
        elif taxon == "Murinae":
            os.system("cd ../family-analyzer; ./bin/familyanalyzer --xreftag protId "+ "../Result/"+settings.dir_name_param+"/"+ settings.xml_name_param +" "+ "MOUSERATNO" +" MOUSE RATNO >"  +"../Result/"+ settings.dir_name_param+"/FA/"+ "Murinae.txt" )




