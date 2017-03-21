#! /usr/bin/env python

import unittest
import sys
sys.path.append('libs')
import functions

class readTest(unittest.TestCase):
    def test_readData():
        assertEqual(a.readData("./tests/calibrated_IBD cohort sample H6-run 1_0_E24_1.xy")[0][0],1299.11)

if __name__ == "__main__":
    unittest.main(verbosity=2)
