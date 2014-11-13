#!/usr/bin/env python
from urllib2 import urlopen
from Bio import Entrez
import sys
import StringIO
import gzip
from ftplib import FTP
import re
import os
import subprocess as sp
import glob
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
                self.samples[sample['gsm']] = sample

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
            #for k,v in data.items():
            #    print k,v
            sample = {'gsm':gsm}
            sample['tax_id'] = int(data['Sample_taxid_ch1'][0])
            sample['sra'] = []
            sample['name'] = data['Sample_title'][0]
            sample['library'] = data['Sample_library_strategy'][0]
            sample['info'] = data['Sample_characteristics_ch1']
            for sra_link in [x for x in data['Sample_relation'] if x.startswith("SRA")]:
                sample['sra'].append(sra_link)
            yield sample
 
    def download_srx(self, srx, outdir):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        ftp = FTP(FTP_ROOT)
        ftp.login()
        rootdir = "/sra/sra-instant/reads/ByExp/sra/SRX/{0}/{1}".format(srx[:6], srx)
        ftp.cwd(rootdir)
        dirs  = []
        ftp.retrlines('LIST', callback=lambda x: dirs.append(x.split(" ")[-1]))
        for dirname in dirs:
            ftp.cwd(os.path.join(rootdir, dirname))
            fnames = []
            ftp.retrlines('LIST', callback=lambda x: fnames.append(x.split(" ")[-1]))
            for fname in fnames:
                local_name = os.path.join(outdir,fname)
                if not os.path.exists(local_name):
                    sys.stderr.write("Downloading {0}...\n".format(local_name))
                    f = open(local_name, "w")
                    ftp.retrbinary(
                                   "RETR {0}".format(os.path.join(rootdir, dirname, fname)),
                                   f.write
                                   )
                    f.close()
                else:
                    sys.stderr.write("{0} is present already...\n".format(local_name))
                yield local_name
   
    def download_sra(self, sra_link, gsm, outdir="./"):
        p = re.compile(r'term=(\w+)')
        m = p.search(sra_link)
        if m:
            srx = m.group(1)
            for fname in self.download_srx(srx, outdir):
                yield fname
                    
        else:
            sys.stderr.write("No SRA link found for {0}\n".format(gsm))

           
    def download(self, outdir="./"):
        outdir = os.path.join(outdir, self.gse)
        for sample in self.samples.values():
            fnames = []
            for fname in self._download_sample(sample, outdir=outdir):
                fnames.append(fname)
            sample['sra_files'] = fnames
            yield sample, fnames



    def _download_sample(self, sample, outdir="."):
        """ Download a sample of the Geo object. Returns the filename.
        :params sample The GSM identifier of the sample to download.
        :type sample string
        :params outdir The directory path to save the downloaded file in.
        :sample outdir string
        """
        for sra_link in sample['sra']:
            for fname in self.download_sra(sra_link, outdir=outdir, gsm=sample['gsm']):
                yield fname

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

        
    def check_sras(self):
        """Check the SRA files of all samples for sanity."""
        #check for vdb-validate
        try:
            sp.check_output('which vdb-validate', shell=True)
        except sp.CalledProcessError:
            print "Could not find 'vdb-validate'-tool in system path; cannot ceck sanity of downloaded files.\n"
            return []
        #check the sample files
        print 'Checking downloaded files for sanity...'
        for sample in self.samples.values():
            for sra_file in sample['sra_files']:
                status =  self._check_sra(sra_file)
                if status == False:
                   print 'File {0} seems to be not ok!\n'.format(sra_file)
    
    
    def _check_sra(self, sra):
        """Check an sra file for sanity.
        :param sra Path to sra file.
        :type  sra string
        :returns boolean to indicate sanity."""
        cmd = "vdb-validate {0}"
        p = sp.Popen(cmd.format(sra),
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
        stdout, stderr = p.communicate()
        message = stderr
        ok = []
        for line in message.splitlines():
            ok.append(line.endswith('ok') or line.endswith("consistent"))
        return not (False in ok)
        
    
    def sras2fastqs(self):
        """Convert the sra files of all samples to fastq files."""
        for sample in self.samples.values():
            for sra_file in sample['sra_files']:
                yield self._sra2fastq(sra_file, sample['name'])

    
    def _sra2fastq(self, sra, name, outdir=self.gse, keep_sra=False):
        """Convert an sra file to a fastq file. Returns a list of the fastq filenames.
        :param sra Filename of the .sra file.
        :type sra string
        :param name GSM identifier of the sample to convert.
        :type name string
        :param outdir Directory store the fastq files in.
        :type outdir string
        """
        try:
            FASTQ_DUMP = "fastq-dump"
            cmd = "{0} --split-files --gzip {1} -O {2}".format(
                                                              FASTQ_DUMP,
                                                              sra,
                                                              outdir,
                                                              )
                                                  
            p = sp.Popen(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
            stdout, stderr = p.communicate()
            if stderr:
                raise Exception(stderr)
        
            sys.stderr.write("Successfully converted {0} to fastq\n".format(sra))
        
            if not keep_sra:
                os.unlink(sra)
            
            base = os.path.splitext(os.path.basename(sra))[0]
            #os.unlink(sra)
            p = re.compile(r'(SRR.+)\.sra')
            m = p.search(sra)
            srr = m.group(1)
  
            fqs = []
            for old_fq in glob.glob(os.path.join(outdir, "*{0}*fastq.gz".format(base))):
                fastq = os.path.join(outdir, "{0}.{1}.fq.gz".format(srr, name))
                os.rename(old_fq, fastq)
                fqs.append(fastq)
    
            return fqs
    
        except Exception as e:
            sys.stderr.write("fastq-dump of {0} failed :(\n".format(sra))
            sys.stderr.write("{0}\n".format(e))
            return []



if __name__ == "__main__":
    #for k,v in Geo.search("Heeringen AND Veenstra").items():
    #    print k,v
    
    #x = Geo("GSE14025")
    x = Geo(open("tests/data/GSE14025_family.soft.gz"))
    for sample in x.samples.values():
        print sample
    

