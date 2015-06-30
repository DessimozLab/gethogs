__author__ = 'traincm'
import sys, getopt


from test_hogs import *
import classes as classes
import os
import pipeline_FA_OMA_cutting
import pipeline_FA_OMA_bottom
import threading
import main as ma


def pipeline(array_param, cluster,dataset):
    print(dataset)
    for param in array_param:
        if cluster:
            arg= str(param)+',"'+str(dataset)+'"'
            cmd= "bsub python -c 'import pipeline_HOG; pipeline_HOG.launch_job("+ arg +")'"
            os.system(cmd)
        else:
            arg= str(param)+',"'+str(dataset)+'"'
            cmd= "python -c 'import pipeline_HOG; pipeline_HOG.launch_job("+ arg +")'"
            os.system(cmd)





def launch_job(param, dataset):
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
    ma.main(set, dataset)
    launch_test(set,dataset)
    pipeline_FA_OMA_bottom.main_pipeline_FA(set)


def main(argv):
    array_param = []
    rerun=False
    cluster=False
    dataset = None
    try:
        opts, args = getopt.getopt(argv,"hp:d:rc" )
    except getopt.GetoptError:
        print('Usage: pipeline_HOG.py [-c] [-r] -p PARAMETER -d DATASET')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h"):
            print('Usage: pipeline_HOG.py [-c] [-r] -p PARAMETER -d DATASET')
            sys.exit()
        elif opt in ("-p"):
            g = arg.split(",")
            for e in g:
                array_param.append(int(e))
        elif opt in ("-r"):
            rerun = True
        elif opt in ("-c"):
            cluster = True
        elif opt in ("-d"):
            if arg != "tiny" and arg!= "big":
                print("wrong dataset")
                sys.exit()
            dataset = arg
    if rerun:
        print("rerun")
        #pipeline_FA_OMA_cutting.main_pipeline_FA("0")
        #pipeline_FA_OMA_cutting.main_pipeline_FA("0_65")


    pipeline(array_param,cluster,dataset)


if __name__ == "__main__":
   main(sys.argv[1:])




