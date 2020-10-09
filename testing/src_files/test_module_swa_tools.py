#!/usr/bin/env python3

import unittest
import os
from datetime import datetime
from Bio import AlignIO

import sys
sys.path.insert(0, "../../sliding_window/src_files")
import swa_tools


class test_shannon_entropy(unittest.TestCase):
    """
        testing the shannon_entropy function
    """
    def test_function_exists(self):
        """
            uses dir([obj]). When the object is a module, 
            dir() returns the names of the module's attributes 
        """
        self.assertIn("shannon_entropy", dir(swa_tools))

    def test_entropy_output(self):
        """
            check that the output makes sense by comparing 
            the output of two sequences passed to shannon_entropy() 
            with values calculated by hand.
        """
        # high=1.366
        high="ATCGATCGAT"
        # low=0
        low="AAAAAAAAAA"
        self.assertAlmostEqual(swa_tools.shannon_entropy(high), 1.366, 3)
        self.assertEqual(swa_tools.shannon_entropy(low), 0)

class test_mean_entropy(unittest.TestCase):
    """
        testing the mean_entropy function
    """
    checkruntime=datetime.now().strftime('%m_%d_%H_%M')
    
    def test_function_exists(self):
        """
            uses dir([obj]). When the object is a module, 
            dir() returns the names of the module's attributes 
        """
        self.assertIn("mean_entropy", dir(swa_tools))
        
    def test_mean_output(self):
        """
            check that a mock file returns a calculated mean entropy of ~.4 
            I calculate the entropy of mock_clustal=.39526, but 
            the function returns .40890
        """     
        alignment = AlignIO.read("../resource_files/mock_clustal.aln","clustal")
        self.assertEqual(alignment.get_alignment_length(), 60)
        self.assertAlmostEqual(swa_tools.mean_entropy(alignment, "testrun"), .4, 1)

    def test_mean_output_file(self):
        """
            mean_output writes a file, check that there really is one
        """
        get_files=os.listdir("../processed_data")
        self.assertIn("{}_testrun_full_entropy.csv".format(self.checkruntime), get_files)

if __name__ == '__main__':
    unittest.main()