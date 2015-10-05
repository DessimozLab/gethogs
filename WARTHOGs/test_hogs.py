__author__ = 'traincm'
import unittest
import filecmp
import gc


class HOG_Integration_Test_tiny(unittest.TestCase):


    def test_hogs_tiny(self):
        self.assertTrue(filecmp.cmp("../Results/OMA_bottom_0/OMA_HOGS_bottom_0_tiny.xml", '../Results/test_reference/OMA_HOGS_bottom_0_tiny.xml'))
        gc.collect()



class HOG_Integration_Test_big(unittest.TestCase):

    def test_hogs_big(self):
        self.assertTrue(filecmp.cmp("../Results/OMA_bottom_0/OMA_HOGS_bottom_0_big.xml", '../Results/test_reference/OMA_HOGS_bottom_0_big.xml'))
        gc.collect()


class HOG_Integration_Test_huge(unittest.TestCase):

    def test_hogs_huge(self):
        self.assertTrue(filecmp.cmp("../Results/OMA_bottom_0/OMA_HOGS_bottom_0_huge.xml", '../Results/test_reference/OMA_HOGS_bottom_0_huge.xml'))
        gc.collect()




def launch_test(set,dataset):
    if dataset == 'big':
        test_big = unittest.TestLoader().loadTestsFromTestCase(HOG_Integration_Test_big)
        unittest.TextTestRunner().run(test_big)

    elif dataset == 'tiny':
        test_tiny = unittest.TestLoader().loadTestsFromTestCase(HOG_Integration_Test_tiny)
        unittest.TextTestRunner().run(test_tiny)

    elif dataset == 'huge':
        test_huge = unittest.TestLoader().loadTestsFromTestCase(HOG_Integration_Test_huge)
        unittest.TextTestRunner().run(test_huge)

    elif dataset == 'insane':
        print("***** IMPOSSIBRU")


