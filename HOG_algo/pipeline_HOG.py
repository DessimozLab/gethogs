__author__ = 'traincm'
import sys, getopt


from test_hogs import *
import classes as classes
import os
import pipeline_FA_OMA_cutting
import pipeline_FA_OMA_bottom
import threading
import main as ma


def pipeline(array_param, cluster,dataset,FA):
    print(dataset)
    for param in array_param:
        if cluster:
            arg= str(param)+',"'+str(dataset)+'"'
            cmd= "bsub python -c 'import pipeline_HOG; pipeline_HOG.launch_job("+ arg +")'"
            os.system(cmd)
        else:
            arg= str(param)+',"'+str(dataset)+'","'+str(FA)+'"'
            cmd= "python -c 'import pipeline_HOG; pipeline_HOG.launch_job("+ arg +")'"
            os.system(cmd)





def launch_job(param, dataset,FA):
    set = classes.Settings()
    set.dir_name_param = "OMA_bottom_" + str(param)
    set.param_merge = param
    print(set.dir_name_param)
    if not os.path.isdir("./Result/" + set.dir_name_param + "/"):
        os.mkdir("../Result/" + set.dir_name_param + "/")
        os.mkdir("../Result/" + set.dir_name_param + "/FA")
        os.mkdir("../Result/" + set.dir_name_param + "/MOBA")
    classes.reset_uniqueId()
    if dataset == 'big':
        set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_big.xml'
    elif dataset == 'tiny':
        set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_tiny.xml'
    elif dataset == 'huge':
        set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_huge.xml'
    ma.main(set, dataset)
    #launch_test(set,dataset)
    if FA == "False":
        return
    pipeline_FA_OMA_bottom.main_pipeline_FA(set)


def main(argv):
    array_param = []
    rerun=False
    cluster=False
    dataset = None
    FA = True
    try:
        opts, args = getopt.getopt(argv,"hp:d:rcf" )
    except getopt.GetoptError:
        print('Usage: pipeline_HOG.py [options] -p [PARAMETER...] -d [DATASET] ')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h"):
            print('\n  Python script for HOGs pipeline inferences \n'
            '\nUsage: pipeline_HOG.py [options] -p [PARAMETER...] -d [DATASET] \n'
            '\nArguments: \n'
            '   -p  List of integers separated by ",", corresponding to parameters used for each pipeline-runs. \n'
            '   -d  Dataset used for inferences, possibilities: "tiny", "big", "huge" (coming soon).  \n'
            '\nOptions: \n'
            '   -c  If specify, use bsub cmd to start a job in a LSF-cluster. Otherwise use bash python cmd. \n'
            '   -r  If specify, re-run family analyzer on OMA current HOGs xml file (top-down method), otherwise use previous family analyzer output.  \n'
            '   -f  If specify, avoid running family analyzer on new HOGs xml file, by default family analyzer is runned.  \n')
            sys.exit()
        elif opt in ("-p"):
            g = arg.split(",")
            for e in g:
                array_param.append(int(e))
        elif opt in ("-r"):
            rerun = True
        elif opt in ("-c"):
            cluster = True
        elif opt in ("-f"):
            FA = False
        elif opt in ("-d"):
            if arg != "tiny" and arg!= "big" and arg!= "huge" :
                print("wrong dataset")
                sys.exit()
            dataset = arg
    if rerun:
        print("rerun")
        #pipeline_FA_OMA_cutting.main_pipeline_FA("0")
        #pipeline_FA_OMA_cutting.main_pipeline_FA("0_65")


    pipeline(array_param,cluster,dataset,FA)


if __name__ == "__main__":
   main(sys.argv[1:])




