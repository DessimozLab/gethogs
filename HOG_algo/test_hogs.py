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
    settings = None
    def test_hogs_tiny(self):
        self.assertTrue(filecmp.cmp("../Result/" + self.settings.dir_name_param + "/" + self.settings.xml_name_param, '../Result/test_reference/OMA_HOG_bottom_tiny.xml'))
        gc.collect()


def launch_test(set):
    #test_normal = unittest.TestLoader().loadTestsFromTestCase(HOG_Integration_Test)
    #unittest.TextTestRunner().run(test_normal)

    HOG_Integration_Test_Tiny.settings = set
    test_tiny = unittest.TestLoader().loadTestsFromTestCase(HOG_Integration_Test_Tiny)
    unittest.TextTestRunner().run(test_tiny)
