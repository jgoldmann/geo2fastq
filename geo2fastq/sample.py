import sys
import re
import os
import glob
from ftplib import FTP
import subprocess as sp
import shutil
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

        
    def fastqs2bam(self, config, outdir, force=False):
        # collect and preformat data
        fastqs = sum(self.fastqs, [])
        name = re.sub(r'[^a-zA-Z1-9_-]', "", self.name)
        bam = os.path.join(outdir,'{0}.{1}.bam'.format(self.gsm, name))
        genome = config['genome_build'][self.tax_id]
        if self.library in config['aligner'].keys():
            aligner = config['aligner'][self.library]
        else:
            aligner = config['aligner']['default']
        genome_dir = config['genome_dir']
        algn_cmd = config['ALIGN_CMD']
        # call mapping method
        self._map(fastqs, bam, genome, aligner, genome_dir, algn_cmd, force)
        self.bam = bam
        return bam

    
    def _map(self, fqs, bam, genome, aligner="", genome_dir="", algn_cmd="", force=False):
        """Map fastq reads to the genome and generate a bam file.
        :param fqs List of fastq filenames.
        :type  fqs list
        :param bam Path to the bam file to generate.
        :type  bam string
        :param genome Genome name to map reads to (as defined in the configfile).
        :type  genome string
        :param aligner Alignment program (as defined in the configfile).
        :type  aligner string
        :param genome_dir Directory of the genome files for the mapper (as definced in configfile).
        :type  genome_dir string
        :param algn_cmd Command structure for the alignment wrapper.
        :type  algn_cmd string
        :param force Re-map even if bamfiles exist already. Defaults to False.
        :type  force boolean
        """
        bams = [] 
        if os.path.exists(bam) and not force:
           sys.stderr.write("{0} already exists, skipping...\n".format(bam))
           return
        for fq in fqs:
            #first, map to intermediate bams
            bname = fq.replace(".fq.gz", "") 
            if os.path.exists(bname+".bam"):
                os.unlink(bname+".bam")
            #delete temporary directory of the bwa mapper
            if os.path.exists(bname):
                shutil.rmtree(bname)
            cmd = algn_cmd.format(fq, bname, genome, aligner, genome_dir)
            sys.stderr.write("Mapping {0} to {1}\n".format(fq, genome))
            p = sp.Popen(cmd, 
                             shell=True, 
                             stderr=sp.PIPE, 
                             stdout=sp.PIPE)
            stdout, stderr = p.communicate() #bwa gives a lot of routine messages on stdout, so don't check for it
            if not os.path.exists("{0}.bam".format(bname)):
                sys.stderr.write("Alignment failed\n")
                sys.stderr.write(stderr)
                sys.exit(1)
            bams.append("{0}.bam".format(bname))
        if len(bams) == 1:
            os.rename(bams[0], bam)
            os.rename(bams[0] + ".bai", bam + ".bai")
        elif len(bams) > 1:
            cmd = "{0} merge -f {1} {2} && samtools index {1}"
            cmd = cmd.format(
                             "samtools",
                             bam,
                             " ".join(bams),
                             )
            ret = sp.call(cmd, shell=True)
            if not ret:
                for bamfile in bams:
                    os.unlink(bamfile)
                    os.unlink("{0}.bai".format(bamfile))
            else:
                sys.stderr.write("Merging failed!\n")
                sys.exit(1)




