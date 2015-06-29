__author__ = 'traincm'

import subprocess
import classes as classes
import os

taxons=["Homininae","Murinae"]

def main_pipeline_FA():
    for taxon in taxons:
        if taxon =='Homininae':
            os.system("cd ../family-analyzer; ./bin/familyanalyzer --xreftag protId "+ "../Result/"+classes.Settings.dir_name_param+"/"+ classes.Settings.xml_name_param +" "+ "HUMANPANTR" +" HUMAN PANTR  >"  +"../Result/"+ classes.Settings.dir_name_param+"/FA/"+ "Homininae.txt" )
        elif taxon == "Murinae":
            os.system("cd ../family-analyzer; ./bin/familyanalyzer --xreftag protId "+ "../Result/"+classes.Settings.dir_name_param+"/"+ classes.Settings.xml_name_param +" "+ "MOUSERATNO" +" MOUSE RATNO >"  +"../Result/"+ classes.Settings.dir_name_param+"/FA/"+ "Murinae.txt" )

