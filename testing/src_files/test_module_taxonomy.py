import unittest
import os
from datetime import datetime
from Bio import AlignIO
from Bio import SeqIO
import subprocess
import shutil
import filecmp

import sys
sys.path.insert(0, "../../taxonomy/src_files")
import prep_id_file


class test_prep_id_file(unittest.TestCase):
    """
        testing the shannon_entropy function
    """
    def test_function_exists(self):
        """
            uses dir([obj]). When the object is a module, 
            dir() returns the names of the module's attributes 
        """
        self.assertIn("get_fasta_headers", dir(prep_id_file))
        
    def test_get_headers(self):
        
        get_sequence="../resource_files/mock_fna_for_idfile.fna"
        get_headers=prep_id_file.get_fasta_headers(get_sequence)
        self.assertEqual(get_headers[0], "JQ765433.1.1505 Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus vallismortis")

    def test_split_headers(self):
        
        get_sequence="../resource_files/mock_fna_for_idfile.fna"
        
        get_headers=prep_id_file.split_description(get_sequence)
        self.assertEqual(get_headers['JQ765433.1.1505'], ['Bacteria', 'Firmicutes', 'Bacilli', 'Bacillales', 'Bacillaceae', 'Bacillus', 'Bacillus$vallismortis'])  

    def test_divert_taxons(self):
        """
        
        """
        pass
        


        
class test_file_writers(unittest.TestCase):
    """
        I'm using the addCleanup() and tearDown() methods here 
        to control the file clutter. Even though the functions 
        being tested are living in /taxonomy and write files 
        to /taxonomy/processed_files, the test functions are 
        writing files to the directories under /testing. 
        That's a lot easier than I thought it would be, and 
        it's probably due to the relative path names 
        in the functions. Happy accident.
    """
    def test_OLDtag_taxons(self):
        """
        
        """
        get_sequence="../resource_files/mock_fna_for_idfile.fna"
        runtime=datetime.now().strftime('%B%d%y_%H_%M')
        tagged_taxon="../processed_files/old_tagged_taxons_{}.txt".format(runtime)
        get_taxon_dict=prep_id_file.split_description(get_sequence)
        useable_headers=prep_id_file.OLDtag_taxons(get_taxon_dict)
        
        with open(tagged_taxon, 'r') as f:
            check_tags=f.readlines()
        
        self.assertEqual(len(check_tags),2)
        self.assertEqual(check_tags[0].strip(), "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus$vallismortis")
        self.assertEqual(check_tags[1].strip(), "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Betaproteobacteriales;f__Burkholderiaceae;g__Burkholderia-Caballeronia-Paraburkholderia;s__Burkholderia$sp.$BRHS/P92")
        
        self.addCleanup(os.remove, tagged_taxon)
        
    def test_tag_taxons(self):
        """
        
        """
        get_sequence="../resource_files/mock_fna_for_idfile.fna"
        runtime=datetime.now().strftime('%B%d%y_%H_%M')
        tagged_taxon="../processed_files/tagged_taxons_{}.txt".format(runtime)
        get_taxon_dict=prep_id_file.split_description(get_sequence)
        useable_headers=prep_id_file.tag_taxons(get_taxon_dict)
        
        
        with open(tagged_taxon, 'r') as f:
            check_tags=f.readlines()
        
        self.assertEqual(len(check_tags),16)
        self.assertEqual(check_tags[0].strip(), "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus$vallismortis$num0")
        self.assertEqual(check_tags[2].strip(), "d__Bacteria;p__Bacteroidetes;c__Ignavibacteria;o__OPB56;f__uncultured$marine$bacterium$num3")        
        
   
        
    def test_weed_fasta(self):
        """
            
        """
        runtime=datetime.now().strftime('%B%d%y_%H_%M')
        tagged_taxon="../processed_files/tagged_taxons_{}.txt".format(runtime)
        subset_records="../processed_files/subset_records_{}.fna".format(runtime)
        get_sequence="../resource_files/mock_fna_for_idfile.fna"
        
        get_taxon_dict=prep_id_file.split_description(get_sequence)
        useable_headers=prep_id_file.tag_taxons(get_taxon_dict)
        
            
        prep_id_file.weed_fastas(get_sequence, useable_headers)
        
        check_subset={x.id:[x.description, x.seq] for x in SeqIO.parse(subset_records, "fasta")}
        
        self.assertEqual(check_subset["JQ765578.1.1444"][0].split()[1], "Root;Bacteria;Proteobacteria;Gammaproteobacteria;Betaproteobacteriales;Burkholderiaceae;Burkholderia-Caballeronia-Paraburkholderia;Burkholderia")
        
        self.assertEqual(check_subset["JQ765578.1.1444"][0].split()[4], "num1")
        
        self.assertEqual(check_subset["JQ765578.1.1444"][1], "AUCG")
        
        self.addCleanup(os.remove, tagged_taxon)
        self.addCleanup(os.remove, subset_records)
        

class test_validate_match(unittest.TestCase):

    def test_script_exists(self):
        """
            check the script exists
        """
        get_files=os.listdir("../../taxonomy/src_files")
        self.assertIn("validate_match_batch.py", get_files)
       

    def test_unclassified(self):
        """
            does validate_match catch the ambiguous words
        """

        checkit=subprocess.run(["python", "../../taxonomy/src_files/validate_match_batch.py", "-i", "../resource_files/validate_folder", "-m", "../resource_files/testing_bad_mapfile.csv"], capture_output=True, text=True)
        spl_folder=checkit.stdout.strip().split("/")[-2]
        spl_output="{}/{}".format(spl_folder, checkit.stdout.strip().split("/")[-1])

        with open("../processed_files/{}".format(spl_output), 'r') as f:
            get_lines=f.readlines()
            self.assertEqual(get_lines[0].strip(),"id,query,blca,confidence,match")
            self.assertEqual(get_lines[1].strip(),"FC000001.01.02,Pretendbacterium bacterium,unavailable,16.3265306122,0")                  # uncultured bacterium,species
            self.assertEqual(get_lines[2].strip(),"FC000002.01.02,Pretendbacterium bacterium2,unavailable,16.3265306122,0")                 # unavailable,species
            self.assertEqual(get_lines[3].strip(),"FC000003.01.02,Pretendbacterium bacterium3,Pretendbacterium metagenome,16.3265306122,0") # metagenome,species
            self.assertEqual(get_lines[4].strip(),"FC000004.01.02,Pretendbacterium bacterium4,unavailable,16.3265306122,0")                 # no genus taxon
            self.assertEqual(get_lines[5].strip(),"FC000005.01.02,Pretendbacterium bacterium5,unavailable,16.3265306122,0")                 # unidentified,species
            self.assertEqual(get_lines[6].strip(),"FC000006.01.02,Pretendbacterium bacterium6,unavailable,16.3265306122,0")                 # no species taxon
            self.assertEqual(get_lines[7].strip(),"FC000007.01.02,Pretendbacterium bacterium7,unavailable,16.3265306122,0")                 # no genus or species taxons
        
        print("removing ../processed_files/{}".format(spl_folder))
        shutil.rmtree("../processed_files/{}".format(spl_folder))

    def test_multiples(self):
        """
            does validate_match_batch handle multiple word taxons
        """

        checkit=subprocess.run(["python", "../../taxonomy/src_files/validate_match_batch.py", "-i", "../resource_files/validate_folder2", "-m", "../resource_files/testing_good_mapfile.csv"], capture_output=True, text=True)
        spl_folder=checkit.stdout.strip().split("/")[-2]
        spl_output="{}/{}".format(spl_folder, checkit.stdout.strip().split("/")[-1])
        
        with open("../processed_files/{}".format(spl_output), 'r') as f:
            get_lines=f.readlines()
            self.assertEqual(get_lines[0].strip(),"id,query,blca,confidence,match")
            self.assertEqual(get_lines[1].strip(),"FC000001.01.02,Pretendbacterium bacterium,Pretendbacterium bacterium,16.3265306122,1")   # regular match
            self.assertEqual(get_lines[2].strip(),"FC000002.01.02,Pretendbacterium bacterium2,Pretendbacterium bacterium2,16.3265306122,1") # species only
            self.assertEqual(get_lines[3].strip(),"FC000003.01.02,Pretendbacterium bacterium3,Pretendbacterium bacterium3,16.3265306122,1") # number after genus, species only
            self.assertEqual(get_lines[4].strip(),"FC000004.01.02,Pretendbacterium bacterium4 SK52 = DSM 20,Pretendbacterium bacterium4,16.3265306122,1") # strain number in reference
            self.assertEqual(get_lines[5].strip(),"FC000005.01.02,Pretendbacterium bacterium5 SK52 = DSM 20,Pretendbacterium bacterium5,16.3265306122,1") # strain number in reference, species only
            
        print("removing ../processed_files/{}".format(spl_folder))
        shutil.rmtree("../processed_files/{}".format(spl_folder))    

    def test_synonym(self):
        """
            does validate_match_batch handle synonyms
        """     
        pass  

class test_rdp_lineage_to_tax(unittest.TestCase):

    def test_script_exists(self):
        """
            check the script exists
        """
        get_files=os.listdir("../../taxonomy/src_files")
        self.assertIn("rdp_lineage_to_tax.py", get_files)
        
    def test_written(self):
        """
            check for files written, and some comparisons
        """

        checkit=subprocess.run(["python", "../../taxonomy/src_files/rdp_lineage_to_tax.py", "-i", "../resource_files/rdp_test_taxonomy.csv", "-o", "test_rdp_taxonomy"], capture_output=True, text=True)
        
        # is the folder there
        self.assertTrue(os.path.exists(os.path.exists("../processed_files/rdp_prep_taxonomy")))
        
        # there should be 2 files in there
        files_in_dir=os.listdir("../processed_files/rdp_prep_taxonomy")
        self.assertEqual(len(files_in_dir), 2)
        
        for x in files_in_dir:
            if x.split('.')[-1]=='txt':
                taxonomy_file=x
          
        # does the test match the provided actual output
        # rdp_team_taxonomy_check can be found on https://github.com/rdpstaff/classifier/issues/18
        self.assertTrue(filecmp.cmp("../resource_files/rdp_team_taxonomy_check.txt", "../processed_files/rdp_prep_taxonomy/{}".format(taxonomy_file)))
        
        shutil.rmtree("../processed_files/rdp_prep_taxonomy")
        
class test_write_qiime_train_db(unittest.TestCase):      

    def test_script_exists(self):
        """
            check the script exists
        """
        get_files=os.listdir("../../taxonomy/src_files")
        self.assertIn("write_qiime_train_db.py", get_files)
        
    def test_splits(self):
        """
            make sure splitting on the whitespace doesn't 
            cause more problems than it's supposed to.
        """
        
        checkit=subprocess.run(["python", "../../taxonomy/src_files/write_qiime_train_db.py", "-i", "../resource_files/spaces.fna", "-o", "../processed_files/test_spaces.fna"], capture_output=True, text=True)
        
        get_processed=os.listdir("../processed_files")
        self.assertIn("test_spaces.fna", get_processed)
        
        with open("../processed_files/test_spaces.fna", 'r') as f:
            lines=f.readlines()
            
        self.assertEqual(lines[0].strip(), ">GY203941.1.1493	Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella_7;unidentified")
        
        os.remove("../processed_files/test_spaces.fna")

if __name__ == '__main__':
    unittest.main()