import unittest
import os
from datetime import datetime
from Bio import AlignIO
from Bio import SeqIO

import sys
sys.path.insert(0, "../../digests/src_files")
import digest
        
class test_simple_cuts(unittest.TestCase):

    def test_1_function_exists(self):
        """
            uses dir([obj]). When the object is a module, 
            dir() returns the names of the module's attributes 
        """
        self.assertIn("simple_cuts", dir(digest))
        
    def test_find_site(self):
        """
            check that Restriction returns the right position
            
            ecori cuts G^AATTC
            
            Restriction counts along the sequence until the cut site is to the left of the current bp, so the right answer is position 24
            
            dummy file returns this:
            {
            'dummy_seq1': {'EcoRI': [24], 'BamHI': []}, 
            'dummy_seq2': {'EcoRI': [], 'BamHI': []}, 
            'dummy_seq3': {'EcoRI': [], 'BamHI': [14]}
            }
        """
        cleave=digest.simple_cuts("../resource_files/mock_sequence.fna", ['EcoRI','BamHI'], False)
        self.assertEqual(cleave['dummy_seq1']['EcoRI'][0], 24)
        self.assertEqual(cleave['dummy_seq3']['BamHI'][0], 14)
        
        # empty list evaluates to false
        self.assertFalse(cleave['dummy_seq2']['EcoRI'])
        self.assertFalse(cleave['dummy_seq1']['BamHI'])
        
class test_strains_cut(unittest.TestCase):
    """
    
    """
    
    def test_1_function_exists(self):
        """
            uses dir([obj]). When the object is a module, 
            dir() returns the names of the module's attributes 
        """
        self.assertIn("strains_cut", dir(digest)) 
        
    def test_return_dict(self):
        """
            dummy sequence file returns this:
            {
                'EcoRI': ['dummy_seq1', 'dummy_seq4'], 
                'BamHI': ['dummy_seq3', 'dummy_seq4']
                }
        """
        enzyme_dict=digest.strains_cut("../../testing/resource_files/mock_sequence.fna", ['EcoRI','BamHI'], False)
        
        self.assertEqual(enzyme_dict['EcoRI'][0], 'dummy_seq1')
        self.assertEqual(len(enzyme_dict['EcoRI']), 2)
        # put one seq with both cut sites
        self.assertEqual(enzyme_dict['EcoRI'][1], enzyme_dict['BamHI'][1])
        
class test_percent_utility(unittest.TestCase):
    """
    
    """
    
    def test_1_function_exists(self):
        """
            uses dir([obj]). When the object is a module, 
            dir() returns the names of the module's attributes 
        """
        self.assertIn("percent_utility", dir(digest))     

    def test_number(self):
        """
            dummy file returns:
            {
                'EcoRI': 0.5, 
                'BamHI': 0.5
                }
        """
        
        get_dict=digest.percent_utility("../resource_files/mock_sequence.fna", ['EcoRI','BamHI'], .1, False)
       
        self.assertEqual(get_dict['EcoRI'], .5)
        
        get_dict2=digest.percent_utility("../resource_files/mock_sequence.fna", ['EcoRI','BamHI'], .6, False)
        
        self.assertFalse(get_dict2)
        
class test_common_enzymes(unittest.TestCase):
    """
    
    """
    
    def test_1_function_exists(self):
        """
            uses dir([obj]). When the object is a module, 
            dir() returns the names of the module's attributes 
        """
        self.assertIn("common_enzymes", dir(digest))      
        
    def test_inersection(self):
        """
        
        """
        mock_list=["../resource_files/mock_sequence.fna"]
        get_intersect=digest.common_enzymes(mock_list, ['EcoRI','BamHI','HindIII'], .4, is_ribo=True)

        self.assertEqual(get_intersect, ['BamHI', 'EcoRI'])
        
class test_calc_fragments(unittest.TestCase):
    """
        got totally confused with the first attempt at this
    """
    
    def test_1_function_exists(self):
        """
            uses dir([obj]). When the object is a module, 
            dir() returns the names of the module's attributes 
        """
        self.assertIn("calc_fragments", dir(digest))
        
    def test_check_digest(self):
        """
            check the recognition sites in the dummy fasta
        """
        enzymes=['EcoRI','BamHI']
        cleave=digest.simple_cuts("../resource_files/mock_sequence.fna", enzymes)
        
        #print(cleave)
        
        for e in enzymes:
            for k,v in cleave.items():
                if e=="EcoRI" and k=="dummy_seq1":
                    #print(v[e])
                    self.assertEqual(v[e][0], 24)
                if e=="BamHI" and k=="dummy_seq3":
                    self.assertEqual(v[e][0], 14)

    def test_build_dict(self):
        """
            check the length of each dummy record is actually 40 bp
        """        
        
        get_seq=[x for x in SeqIO.parse("../resource_files/mock_sequence.fna", "fasta")]
        seq_length={x.id:len(x.seq) for x in get_seq}
        for k,v in seq_length.items():
            self.assertEqual(v, 40)
        
    def test_build_range(self):
        """
            check that the calculated fragments from dummy records make sense
        """        
        eco_fragments=digest.calc_fragments([24], 40)
        self.assertEqual(eco_fragments, [24, 16])
        bam_fragments=digest.calc_fragments([14], 40)
        self.assertEqual(bam_fragments, [14, 26])
        

    def test_total_fragments(self):
        """
            alltogether now, check that the returned fragments 
            from all sequences with the recognition 
            sites of an enzyme make sense
        """        
        
        complete=digest.digest_fragments("../resource_files/mock_sequence.fna", ['EcoRI','BamHI'], verbose=True)
        
        self.assertEqual(complete['EcoRI'], {24: 2, 16: 2})
        self.assertEqual(complete['BamHI'], {14: 2, 26: 2})
        
        
        
        
        


if __name__ == '__main__':
    unittest.main()