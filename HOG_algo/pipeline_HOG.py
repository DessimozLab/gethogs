__author__ = 'traincm'
import sys, getopt


from test_hogs import *
import classes as classes
import os
#import pipeline_FA_OMA_cutting
import pipeline_FA_OMA_bottom
import threading
import main as ma


def pipeline(array_param, cluster):
    for param in array_param:
        if cluster:
            os.system("bsub python -c 'import pipeline_HOG; pipeline_HOG.launch_job("+ str(param)+","+str(cluster)+")'")
        else:
            os.system(" python -c 'import pipeline_HOG; pipeline_HOG.launch_job("+ str(param)+","+str(cluster)+")'")





def launch_job(param,cluster):
    print("new thread")
    set = classes.Settings()
    set.dir_name_param = "OMA_bottom_" + str(param)
    set.param_merge = param
    print(set.dir_name_param)
    if not os.path.isdir("./Result/" + set.dir_name_param + "/"):
        os.mkdir("../Result/" + set.dir_name_param + "/")
        os.mkdir("../Result/" + set.dir_name_param + "/FA")
        os.mkdir("../Result/" + set.dir_name_param + "/MOBA")
    classes.reset_uniqueId()
    set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_tiny.xml'
    ma.main(set, 'tiny')
    launch_test(set)
    pipeline_FA_OMA_bottom.main_pipeline_FA(set)


def main(argv):
    array_param = []
    rerun=False
    cluster=False
    try:
        opts, args = getopt.getopt(argv,"hp:r:c" )
    except getopt.GetoptError:
        print('Usage: pipeline_HOG.py [-c] [-r] -p PARAMETER')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h"):
            print('Usage: pipeline_HOG.py [-c] [-r] -p PARAMETER')
            sys.exit()
        elif opt in ("-p"):
            g = arg.split(",")
            for e in g:
                array_param.append(int(e))
        elif opt in ("-r"):
            rerun = True
        elif opt in ("-r"):
            cluster = True
    if rerun:
        print("rerun")
        #pipeline_FA_OMA_cutting.main_pipeline_FA("0")
        #pipeline_FA_OMA_cutting.main_pipeline_FA("0_65")


    pipeline(array_param,cluster)


if __name__ == "__main__":
   main(sys.argv[1:])




