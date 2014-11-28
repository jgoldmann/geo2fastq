#!/usr/bin/env python
from urllib2 import urlopen
from Bio import Entrez
import sys
import StringIO
import gzip
import os
import subprocess as sp
from sample import Sample
import pp
import pdb


GEOFTP_URLBASE = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/{0}/{0}_family.soft.gz"
FTP_ROOT = "ftp-trace.ncbi.nlm.nih.gov"


class Geo:
    def __init__(self, gse="", email=""):
        """ Create Geo object to search GEO and parse GEO experiment data
        
        :param gse: GEO accession or filehandle to GEO soft file
        :type accession: string or file
        """
        self.gse = gse
        self.email = email
        self.samples = {}
        if gse:
            try:
                callable(gse.read)
                gen = self.get_sample_info(fo=gse)
            except:
                self.gse = gse
                gen = self.get_sample_info(accession=self.gse) 
            for sample in gen:
                self.samples[sample.gsm] = sample
                

    @classmethod
    def search(self, term):
        """ Search NCBI GEO and return a dictionary of GEO accessions
        describing series (GSE) and samples (GSM).
        
        :param term: search term(s)
        :type term: string
        """
        Entrez.email = "s.vanheeringen@ncmls.ru.nl"  #TODO: parameterize
        handle = Entrez.esearch("gds", term)
        record = Entrez.read(handle)
        results = {}
        for gid in record['IdList']:
            handle = Entrez.esummary(db="gds", id=gid)
            r = Entrez.read(handle)[0]
            if r['entryType'] == 'GSE':
                gse = "GSE{0}".format(r['GSE'])
                results[gse] = {
                                'title': r['title'].encode('utf8'),
                                'samples':{}
                                }
                for sample in r['Samples']:
                    results[gse]['samples'][sample['Accession']] = sample['Title'].encode('utf8')
        return results
        
    
    def get_sample_info(self, accession=None, fo=None):
        """ Accepts a single accession or a list of accessions and returns
        an generator of detailed sample information.
       
        Supply either an accession or a SOFT filehandle.    
     
        :param accession: GEO accession
        :type accession: string
        :params fo: filehandle to gzipped SOFT file
        :type fo: file
        """
        if not fo:
            if not accession:
                sys.stderr.write("Please supply either accession or fo")
                sys.exit(1)
            # Dowload GEO SOFT file
            soft = GEOFTP_URLBASE.format(accession)
            try:
                fh = urlopen(soft)
            except:
                sys.stderr.write("Could not open SOFT URL {0}\n".format(soft))
                sys.exit(1)
            fo = StringIO.StringIO(fh.read())
        
        # Parse gzipped SOFT file
        g = gzip.GzipFile(fileobj=fo)
        record = self._soft_read(g)
        self.gse = record['SERIES'].keys()[0]
        for gsm, data in record['SAMPLE'].items():
            sample = Sample(gsm, data)
            yield sample


    def _soft_read(self, fh):
        """ Parses a filehandle of a SOFT formatted file
        """
        soft = {}
        current = soft
        for line in fh:
            try:
                key, val = line[1:].strip().split(" = ")
                if line.startswith("^"):
                    soft.setdefault(key, {})[val] = {}
                    current =  soft[key][val]
                else:
                    current.setdefault(key, []).append(val)
            except:
                pass
        return soft

 
    def download(self, outdir='./'):
        outdir = os.path.join(outdir, self.gse)
        server = pp.Server(len(self.samples))
        for sample in self.samples.values():
            sra_files = server.submit(sample.download, (outdir,), modules=('geo2fastq.sample',))()
            if 'sra_files' not in dir(sample):  #make sure that the sra filenames arrive in the sample instance of the main process
                sample.sra_files = sra_files
        return []
        
        

        
    def check_sras(self):
        """Check the SRA files of all samples for sanity."""
        #check for vdb-validate
        try:
            sp.check_output('which vdb-validate', shell=True)
        except sp.CalledProcessError:
            print "Could not find 'vdb-validate'-tool in system path; cannot ceck sanity of downloaded files.\n"
            return []
        print 'Checking downloaded files for sanity...'
        server = pp.Server(len(self.samples))
        ok = [server.submit(sample.check_sras)() for sample in self.samples.values()]
        return ok

        
    
    def sras2fastqs(self, keep_sra, outdir='./'):
        """Convert the sra files of all samples to fastq files."""
        outdir = os.path.join(outdir, self.gse)
        #check for fastq-dump
        try:
            sp.check_output('which fastq-dump', shell=True)
        except sp.CalledProcessError:
            print "Could not find 'fastq-dump'-tool in system path; cannot ceck sanity of downloaded files.\n"
            return []
        #else, check the sample files
        print 'Converting files to fastq...'
        server = pp.Server(len(self.samples))
        for sample in self.samples.values():
            fastqs = server.submit(sample.sras2fastqs, (outdir, keep_sra))()
            sample.fastqs = fastqs
        return []
    
    
    def fastqs2bams(self, config, outdir, force=False):
        print 'Starting Alignments...'
        outdir = os.path.join(outdir, self.gse)
        server = pp.Server(len(self.samples))
        bams = []
        for sample in self.samples.values():
            bam = server.submit(sample.fastqs2bam, (config, outdir, force), modules=('re',))()
            sample.bam = bam
            bams.append(bam)
        print 'Finished with alining.'
        return bams


