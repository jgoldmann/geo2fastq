from geo2fastq import Geo
import os
import shutil

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
            
        fname_generator = g.download()
        for sample, fname in fname_generator:
            for path in fname:
                assert(os.path.exists(path))
                
        if not files_present_before:
            shutil.rmtree("./GSE14025")
        

