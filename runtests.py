from geo2fastq import Geo

class TestClass:
    def test_Geo_search():
        results = Geo.search('GSE14025')
        assert(results["GSE14025"] == "A Hierarchy of H3K4me3 and H3K27me3 Acquisition in Spatial Gene Regulation in Xenopus Embryos")
        assert(results.has_key("GSM352204"))
        assert(results.has_key("GSM352202"))
        assert(results.has_key("GSM352203"))
        assert(results.has_key("GSM419463"))
