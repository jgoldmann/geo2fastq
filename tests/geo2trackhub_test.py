from geo2fastq import Geo
import os
import shutil
from geo2fastq.convert import sra2fastq, fastq2bam, bam2bw
from geo2fastq.config import config
import glob
import re


config = config()


class TestClass:
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
        files_present_before = os.path.exists("./GSE14025")
        if files_present_before:
            shutil.rmtree("./GSE14025")
            
        fname_generator = g.download() #takes about 20 minutes
        for sample, fname in fname_generator:
            for path in fname:
                assert(os.path.exists(path))
                
        if not files_present_before:
            shutil.rmtree("./GSE14025")
       
       
    def test_sra2fastq(self):
        gse = 'GSE14025'
        g = Geo(gse)
        sample = g.samples[g.samples.keys()[1]]
        sra_fname = g._download_sample(sample).next()
        fqs = sra2fastq(sra_fname, 
                        sample['gsm'], 
                        gse)
        for fq in fqs:
            assert(os.path.exists(fq))
        # do not delete the output file at this stage
      
      
    def test_fastq2bam(self):
        gse = 'GSE14025'
        g = Geo(gse)
        sample = g.samples[g.samples.keys()[1]]
        fqs = glob.glob(os.path.join(gse, "*{0}*.fq.gz".format(sample['gsm'])))
        name = re.sub(r'[^a-zA-Z1-9_-]', "", sample['name'])
        bam = os.path.join(gse, "{0}.{1}.bam".format(sample['gsm'], name))
        aligner = config['aligner'].setdefault(sample['library'], config['aligner']['default'])
        genome_dir = config['genome_dir']
        genome = config['genome_build'][sample['tax_id']] 
        #TODO: in future, test with some less exotic genome (now xenTro3beta)
        fastq2bam(fqs, bam, genome, aligner, genome_dir)
        assert(os.path.exists(bam))
        assert(os.path.exists("{0}.bai".format(bam)))
        for fq in fqs:
            os.unlink(fq)
        
        
    def test_bam2bw(self):
        gse = 'GSE14025'
        g = Geo(gse)
        sample = g.samples[g.samples.keys()[1]]
        bam = 'GSE14025/GSM352202.H3K4me3_ChIPSeq.bam' 
        bw = bam.replace(".bam", ".bw")
        
        bam2bw(bam, bw, library=sample['library']) # not functional yet!
        assert(os.path.exists(bw))
    
       




