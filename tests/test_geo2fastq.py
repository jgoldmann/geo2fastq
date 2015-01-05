from unittest import TestCase
import geo2fastq
from geo2fastq import Geo
import os
import urllib

config = geo2fastq.config.config()
test_gse = "GSE54401"

# run tests with test debugging mode: nosetests --pdb tests/test_geo2fastq.py

class TestClass(TestCase):
    def test_Geo_search(self):
        for gse_out, info in Geo.search(test_gse).items():
           assert(info['title'] == "A Hierarchy of H3K4me3 and H3K27me3 Acquisition in Spatial Gene Regulation in Xenopus Embryos")
           for gsm,title in info['samples'].items():
               assert(gsm.startswith('GSM'))
           assert(info['samples'].has_key("GSM352204"))



    def test_Geo_get_sample_info(self):
        g = Geo(test_gse, config)
        
        info = g.get_sample_info(test_gse)
        for i in info:
            assert(i.name != '')

        
    def test_Geo_download(self):
        g = Geo(test_gse, config)
        g.download()
        all_sras = []
        for sample in g.samples.values():
            all_sras.append(sum(sample.sra_files.values(),[]))
        for sra in sum(all_sras,[]):
            assert(os.path.exists(sra))
                
       
    def test_sra2fastq(self):
        g = Geo(test_gse, config)
        g.download()
        all_fastqs = []
        g.sras2fastqs(keep_sra=True)
        for sample in g.samples.values():
            all_fastqs.append(sum(sample.fastqs, []))
        for fastq in sum(all_fastqs, []):
            assert(os.path.exists(fastq))
        
      
    def test_fastqs2bams(self):
        g = Geo('GSE14025', config) 
        g.download()
        g.sras2fastqs(keep_sra = True)
        bams = g.fastqs2bams(config)
        for bam in bams:
            assert(os.path.exists(bam))
        
        
    def test_bams2bws(self):
        g = Geo('GSE14025', config)
        g.download()
        g.sras2fastqs(keep_sra = True)
        g.fastqs2bams(config)
        bws = g.bams2bws()
        for bw in bws:
            assert(os.path.exists(bw))
            
    
    def test_generate_trackhub(self):
        g = Geo('GSE14025', config)
        g.download()
        g.sras2fastqs(keep_sra = True)
        g.fastqs2bams(config)
        g.bams2bws()
        url = g.generate_trackhub(config)
        url_basename = url[0:url.rfind('/')]+'/'
        hub_file = urllib.urlopen(url)
        # check for correct gse:
        gse = hub_file.readline().strip().split(' ')[1]
        assert(g.gse == gse)
        
        # check genomes file
        hublines = hub_file.readlines()
        hub_file.close()
        genomes_filename = hublines[2].strip().split(' ')[1]
        genomes_file = urllib.urlopen(url_basename + genomes_filename)
        genomeLines = genomes_file.readlines()
        genomes_file.close()
        
        # check trackDb file
        for genomeLine in genomeLines:
            if genomeLine.startswith('trackDb '):
                trackdb = genomeLine.strip().split(' ')[1]
        trackdb_filename = trackdb[trackdb.rfind('/')+1:len(trackdb)]
        trackdb_file = urllib.urlopen(url_basename + trackdb_filename)
        
        # check tracks
        dataUrls = []
        for line in trackdb_file.readline():
            if line.startswith("bigDataUrl "):
                dataUrls.apped(line.split(' ')[1])
        trackdb_file.close()
        for dataUrl in dataUrls:
            dataHandle = urllib.urlopen(dataUrl)
            dataHandle.close()
        
    
        
        




