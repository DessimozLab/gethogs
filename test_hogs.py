__author__ = 'traincm'
import unittest
from main import main
import filecmp

class HOG_Integration_Test(unittest.TestCase):
    def test_hogs_big(self):
        resfn = '/tmp/testrun.xml'
        main(resfn, 'big')
        self.assertTrue(filecmp.cmp(resfn, 'Result/OMA_HOG_bottom_none_copy.xml'))

    def test_hogs_tiny(self):
        resfn = '/tmp/testrun_tiny.xml'
        main(resfn, 'tiny')
        self.assertTrue(filecmp.cmp(resfn, 'Result/OMA_HOG_bottom_tiny.xml'))
