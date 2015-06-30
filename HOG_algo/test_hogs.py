__author__ = 'traincm'
import unittest
import filecmp
import gc


class HOG_Integration_Test(unittest.TestCase):
    settings = None
    dataset = None
    def test_hogs_tiny(self):
        self.assertTrue(filecmp.cmp("../Result/" + self.settings.dir_name_param + "/" + self.settings.xml_name_param, '../Result/test_reference/OMA_HOG_bottom_'+ self.dataset +'.xml'))
        gc.collect()

def launch_test(set,dataset):
    if dataset == 'big':
        set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_big.xml'
    elif dataset == 'tiny':
        HOG_Integration_Test.settings = set
        HOG_Integration_Test.dataset = dataset
        test_tiny = unittest.TestLoader().loadTestsFromTestCase(HOG_Integration_Test)
        unittest.TextTestRunner().run(test_tiny)


