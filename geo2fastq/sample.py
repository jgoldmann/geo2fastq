import sys
import re
import os
import glob
from ftplib import FTP
import subprocess as sp
import pdb

GEOFTP_URLBASE = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/{0}/{0}_family.soft.gz"
FTP_ROOT = "ftp-trace.ncbi.nlm.nih.gov"


class Sample:
    
    def __init__(self, gsm, data):
        self.gsm = gsm
        self.tax_id = int(data['Sample_taxid_ch1'][0])
        self.name = data['Sample_title'][0]
        self.library = data['Sample_library_strategy'][0]
        self.info = data['Sample_characteristics_ch1']
        self.sra = []
        for sra_link in [x for x in data['Sample_relation'] if x.startswith("SRA")]:
            self.sra.append(sra_link)
            
        
    def download(self, outdir):
        self.sra_files = {}
        sra_files = {}
        for sra in self.sra:
            files = self.download_sra(sra,outdir)
            self.sra_files[sra] = files
            sra_files[sra] = files
        return sra_files
        
        
    def download_sra(self, sra_link, outdir):
        p = re.compile(r'term=(\w+)')
        m = p.search(sra_link)
        if m:
            srx = m.group(1)
            return [fname for fname in self.download_srx(srx, outdir)]
        else:
            sys.stderr.write("No SRA link found for {0}\n".format(self.gsm))


    def download_srx(self, srx, outdir):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        ftp = FTP(FTP_ROOT)
        ftp.login()
        rootdir = "/sra/sra-instant/reads/ByExp/sra/SRX/{0}/{1}".format(srx[:6], srx)
        ftp.cwd(rootdir)
        dirs  = []
        ftp.retrlines('LIST', callback=lambda x: dirs.append(x.split(" ")[-1]))
        #pdb.set_trace()
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
  
    
    def check_sras(self):
        status_ok = []
        sra_files = sum(self.sra_files.values(),[]) # hack to ensure that instead of list of lists, we will get a list
        for sra_file in sra_files:
            status_ok.append(self._check_sra(sra_file))
        return not (False in status_ok)
    
    
    def _check_sra(self, sra_file):
        cmd = "vdb-validate {0}"
        p = sp.Popen(cmd.format(sra_file),
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
        stdout, stderr = p.communicate()
        message = stderr
        ok = []
        for line in message.splitlines():
            ok.append(line.endswith('ok') or line.endswith("consistent"))
        if False in ok:
            print "Warning: File {0} of samples {1} seems to be not ok!!".format(sra_file, self.gsm)
        return not (False in ok)
    
        
    def sras2fastqs(self, outdir, keep_sra=False):
        self.fastqs = []
        fastqs = []
        sra_files = sum(self.sra_files.values(),[]) # hack to ensure that instead of list of lists, we will get a list
        for sra_file in sra_files:
            fqs = self._sra2fastq(sra_file, outdir, keep_sra)
            self.fastqs.append(fqs)
            fastqs.append(fqs)
        return fastqs
        
    
    def _sra2fastq(self, sra_file, outdir, keep_sra):
        try:
            FASTQ_DUMP = "fastq-dump"
            cmd = "{0} --split-files --gzip {1} -O {2}".format(
                                                              FASTQ_DUMP,
                                                              sra_file,
                                                              outdir,
                                                              )
            if not glob.glob('{0}*.fq.gz'.format(sra_file[0:-4])):
                p = sp.Popen(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
                stdout, stderr = p.communicate()
                if stderr:
                    raise Exception(stderr)
                sys.stderr.write("Successfully converted {0} to fastq\n".format(sra_file))
                if not keep_sra:
                    os.unlink(sra_file)
                base = os.path.splitext(os.path.basename(sra_file))[0]
                p = re.compile(r'(SRR.+)\.sra')
                m = p.search(sra_file)
                srr = m.group(1)
                fqs = []
                for old_fq in glob.glob(os.path.join(outdir, "*{0}*fastq.gz".format(base))):
                    fastq = os.path.join(outdir, "{0}.{1}.fq.gz".format(srr, self.name))
                    os.rename(old_fq, fastq)
                    fqs.append(fastq)
                return fqs
            else: #if fastq exists
                print 'File {0} is converted to fastq already, skipping.'.format(sra_file)
                fqs = glob.glob(sra_file[0:-3]+'*.fq*')
                return fqs
        except Exception as e:
            sys.stderr.write("fastq-dump of {0} failed :(\n".format(sra_file))
            sys.stderr.write("{0}\n".format(e))
            return []    

        






