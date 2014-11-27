from unittest import TestCase
import geo2fastq
from geo2fastq import Geo
import os

config = geo2fastq.config.config()

# run tests with test debugging mode: nosetests --pdb tests/test_geo2fastq.py

class TestClass(TestCase):
    def test_Geo_search(self):
        gse = 'GSE14025'
        for gse_out, info in Geo.search(gse).items():
           #print "{0}\t{1}".format(gse, info['title'])
           assert(info['title'] == "A Hierarchy of H3K4me3 and H3K27me3 Acquisition in Spatial Gene Regulation in Xenopus Embryos")
           for gsm,title in info['samples'].items():
               #print "  {0}\t{1}".format(gsm, title)
               assert(gsm.startswith('GSM'))
           assert(info['samples'].has_key("GSM352204"))



    def test_Geo_get_sample_info(self):
        gse = 'GSE14025'
        g = Geo(gse)
        
        info = g.get_sample_info(gse)
        for i in info:
            assert(i.name != '')

        
    def test_Geo_download(self):
        gse = 'GSE14025'
        g = Geo(gse)
        g.download()
        all_sras = []
        for sample in g.samples.values():
            all_sras.append(sum(sample.sra_files.values(),[]))
        for sra in sum(all_sras,[]):
            assert(os.path.exists(sra))
                
       
    def test_sra2fastq(self):
        gse = 'GSE14025'
        g = Geo(gse)
        g.download()
        all_fastqs = []
        g.sras2fastqs(keep_sra=True)
        for sample in g.samples.values():
            all_fastqs.append(sum(sample.fastqs, []))
        for fastq in sum(all_fastqs, []):
            assert(os.path.exists(fastq))
        
      
      
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
    
       




