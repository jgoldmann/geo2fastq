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
import trackhub
import trackhub.upload
import pdb


GEOFTP_URLBASE = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/{0}/{0}_family.soft.gz"
FTP_ROOT = "ftp-trace.ncbi.nlm.nih.gov"


class Geo:
    def __init__(self, gse="", config={}):
        """ Create Geo object  to search GEO and parse GEO experiment data
        
        :param gse: GEO accession or filehandle to GEO soft file
        :type accession: string or file
        :param config: a dictionary of configuration settings. must contain entry for "email".
        :type accession: dictionary
        """
        self.gse = gse
        self.email = config['email']
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
    def search(self, term, config):
        """ Search NCBI GEO and return a dictionary of GEO accessions
        describing series (GSE) and samples (GSM).
        
        :param term: search term(s)
        :type term: string
        """
        Entrez.email = config['email']
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
        self.record = record
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
    
    
    def fastqs2bams(self, config, outdir='./', force=False):
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

    
    def bams2bws(self,force=False):
        print 'Starting bigwig conversions...'
        server = pp.Server(len(self.samples))
        bws = []
        for sample in self.samples.values():
            bw = server.submit(sample.bam2bw, (force,), modules=('subprocess as sp',))() #parallel
            #bw = sample.bam2bw(force) #sequentially
            sample.bw = bw
            bws.append(bw)
        print 'Done with generating bigwigs.'
        return bws


    def generate_trackhub(self, config, outdir='./'):
        hub_path = os.path.join(os.getenv('USER'), self.gse)
        hub_path_local =  os.path.join(outdir, self.gse)
        hub = trackhub.Hub(hub = self.gse,
                           short_label = self.gse,
                           long_label =  self.record['SERIES'][self.gse]['Series_title'],
                           email = self.email,
                           )
        hub.local_fn =  os.path.join(hub_path_local,                    '{0}.hub.txt'.format(self.gse))
        hub.remote_fn = os.path.join(config['HUB_LOCALBASE'], hub_path, '{0}.hub.txt'.format(self.gse))
        hub.url =       os.path.join(config['HUB_URLBASE'],   hub_path, '{0}.hub.txt'.format(self.gse))
        genomes_in_samples = set()
        tracks = []
        for sample in self.samples.values():
            #collect tracks
            track = sample.generate_track_object(config, self.gse)
            sample.track = track
            tracks.append(track)
            #collect genomes
            genomes_in_samples.add(config['genome_build'][sample.tax_id])
                
        trackdb = trackhub.TrackDb()
        trackdb.local_fn =  os.path.join(hub_path_local,                    '{0}.trackDb.txt'.format(self.gse))
        trackdb.remote_fn = os.path.join(config['HUB_LOCALBASE'], hub_path, '{0}.trackDb.txt'.format(self.gse))
        trackdb.url =       os.path.join(config['HUB_URLBASE'],   hub_path, '{0}.trackDb.txt'.format(self.gse))
        trackdb.add_tracks(tracks)
        #crude hack: add every track to every genome
        genomes = [trackhub.Genome(genome, trackdb) for genome in genomes_in_samples]
        genomesFile = trackhub.GenomesFile(genomes)
        genomesFile.local_fn  = os.path.join(hub_path_local,                    '{0}.genomes.txt'.format(self.gse))
        genomesFile.remote_fn = os.path.join(config['HUB_LOCALBASE'], hub_path, '{0}.genomes.txt'.format(self.gse))
        genomesFile.url       = os.path.join(config['HUB_URLBASE'],   hub_path, '{0}.genomes.txt'.format(self.gse))
        hub.add_genomes_file(genomesFile)
        
        #upload tracks
        hub.render()
        for track in trackdb.tracks:
            trackhub.upload.upload_track(track=track, host = 'localhost', user = os.getenv("USER"))
        trackhub.upload.upload_hub(hub = hub, host = 'localhost', user = os.getenv("USER"))
        print 'Generated a trackhub at: {0}'.format(hub.url)
        return hub.url


