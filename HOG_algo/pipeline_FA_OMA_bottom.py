__author__ = 'traincm'

import subprocess
import classes as classes
import os

taxons=["Homininae","Murinae"]



def write_output_FA_into_file(fn, taxon_arg, taxon):
    cmd = ["../family-analyser/bin/familyanalyzer", "--xreftag", "omaId", "../Result/"+classes.Settings.dir_name_param+ classes.Settings.xml_name_param, taxon]
    for arg in taxon_arg:
        cmd.append(arg)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
    out, err = p.communicate()
    with open("../Result/"+ classes.Settings.dir_name_param+"/FA/"+ fn,"w+") as f:
        for line in out:
            f.write(str(line))
        f.close()

def main_pipeline_FA():
    for taxon in taxons:
        if taxon =='Homininae':
            os.system("cd ../family-analyzer; bin/familyanalyzer --xreftag protId "+ "../Result/"+classes.Settings.dir_name_param+"/"+ classes.Settings.xml_name_param +" "+ "HUMANPANTR" +" HUMAN PANTR  >"  +"../Result/"+ classes.Settings.dir_name_param+"/FA/"+ "Homininae.txt" )
        elif taxon == "Murinae":
            os.system("cd ../family-analyzer; bin/familyanalyzer --xreftag protId "+ "../Result/"+classes.Settings.dir_name_param+"/"+ classes.Settings.xml_name_param +" "+ "MOUSERATNO" +" MOUSE RATNO >"  +"../Result/"+ classes.Settings.dir_name_param+"/FA/"+ "Murinae.txt" )

