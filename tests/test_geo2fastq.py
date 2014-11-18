from unittest import TestCase
import geo2fastq
from geo2fastq import Geo
import os
#import shutil
#import glob
#import re

config = geo2fastq.config.config()


class TestClass(TestCase):
    def test_Geo_search(self):
        gse = 'GSE14025'
        results = Geo.search(gse)
        assert(results[gse]['title'] == "A Hierarchy of H3K4me3 and H3K27me3 Acquisition in Spatial Gene Regulation in Xenopus Embryos")
        assert(results[gse]['samples'].has_key("GSM352204"))
        assert(results[gse]['samples'].has_key("GSM352202"))
        assert(results[gse]['samples'].has_key("GSM352203"))
        assert(results[gse]['samples'].has_key("GSM419463"))


    def test_Geo_get_sample_info(self):
        gse = 'GSE14025'
        g = Geo(gse)
        
        info = g.get_sample_info(gse)
        for i in info:
            assert(i.has_key('name'))

        
    def test_Geo_download(self):
        gse = 'GSE14025'
        g = Geo(gse)
        fname_generator = g.download() #takes about 20 minutes
        for sample, fname in fname_generator:
            for path in fname:
                assert(os.path.exists(path))
                
       
    def test_sra2fastq(self):
        gse = 'GSE14025'
        g = Geo(gse)
        fname_generator = g.download()
        for sample, fname in fname_generator:
            for path in fname:
                assert(os.path.exists(path))
        all_fqs = [fqs for fqs in g.sras2fastqs(keep_sra=True)]
        for fqs in all_fqs:
            for fq in fqs:
                assert(os.path.exists(fq))
        
      
      
#    def test_fastq2bam(self):
#        # depends on test_sra2fastq
#        gse = 'GSE14025'
#        g = Geo(gse)
#        sample = g.samples[g.samples.keys()[1]]
#        fqs = glob.glob(os.path.join(gse, "*{0}*.fq.gz".format(sample['gsm'])))
#        name = re.sub(r'[^a-zA-Z1-9_-]', "", sample['name'])
#        bam = os.path.join(gse, "{0}.{1}.bam".format(sample['gsm'], name))
#        aligner = config['aligner'].setdefault(sample['library'], config['aligner']['default'])
#        genome_dir = config['genome_dir']
#        genome = config['genome_build'][sample['tax_id']] 
#        algn_cmd = config['ALIGN_CMD']
#        #TODO: in future, test with some less exotic genome (now xenTro3beta)
#        geo2fastq.convert.fastq2bam(fqs, bam, genome, aligner, genome_dir, algn_cmd)
#        assert(os.path.exists(bam))
#        assert(os.path.exists("{0}.bai".format(bam)))
#        for fq in fqs:
#            os.unlink(fq)
#        
#        
#    def test_bam2bw(self):
#        # depends on test_fastq2bam
#        gse = 'GSE14025'
#        g = Geo(gse)
#        sample = g.samples[g.samples.keys()[1]]
#        bam = 'GSE14025/GSM352202.H3K4me3_ChIPSeq.bam' 
#        bw = bam.replace(".bam", ".bw")
#        
#        geo2fastq.convert.bam2bw(bam, bw, library=sample['library']) 
#        assert(os.path.exists(bw))
    
       



