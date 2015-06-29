__author__ = 'traincm'
import unittest
from main import main
import filecmp
import gc
import classes as classes

class HOG_Integration_Test(unittest.TestCase):
    def test_hogs_big(self):
        classes.reset_uniqueId()
        resfn = '/tmp/testrun'+ str(classes.param_merge) +'.xml'
        main(resfn, 'big')
        self.assertTrue(filecmp.cmp(resfn, 'Result/OMA_HOG_bottom_none_copy.xml'))
        gc.collect()


class HOG_Integration_Test_Tiny(unittest.TestCase):
    def test_hogs_tiny(self):
        classes.reset_uniqueId()
        classes.Settings.xml_name_param = "OMA_HOGS_bottom_"+ str(classes.Settings.param_merge) +'_tiny.xml'
        main(classes.Settings.xml_name_param, 'tiny')
        self.assertTrue(filecmp.cmp("../Result/" + classes.Settings.dir_name_param + "/" + classes.Settings.xml_name_param, '../Result/test_reference/OMA_HOG_bottom_tiny.xml'))
        gc.collect()


def launch_test():
    #test_normal = unittest.TestLoader().loadTestsFromTestCase(HOG_Integration_Test)
    #unittest.TextTestRunner().run(test_normal)

    test_tiny = unittest.TestLoader().loadTestsFromTestCase(HOG_Integration_Test_Tiny)
    unittest.TextTestRunner().run(test_tiny)
