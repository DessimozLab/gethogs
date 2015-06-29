__author__ = 'traincm'


from test_hogs import *
import classes as classes
import os
import pipeline_FA_OMA_cutting
import pipeline_FA_OMA_bottom


def pipeline(array_param):
    for param in array_param:
        classes.Settings.dir_name_param = "OMA_bottom_" + str(param)
        classes.Settings.param_merge = param
        if not os.path.isdir("./Result/" + classes.Settings.dir_name_param + "/"):
            os.mkdir("../Result/" + classes.Settings.dir_name_param + "/")
            os.mkdir("../Result/" + classes.Settings.dir_name_param + "/FA")
            os.mkdir("../Result/" + classes.Settings.dir_name_param + "/MOBA")
        launch_test()
        pipeline_FA_OMA_bottom.main_pipeline_FA()


#pipeline_FA_OMA_cutting.main_pipeline_FA("0")
#pipeline_FA_OMA_cutting.main_pipeline_FA("0_65")

pipeline([1, 8])




