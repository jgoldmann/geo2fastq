#!/usr/bin/env python
import sys
import os
#import shutil
import subprocess as sp


#def fastq2bam(fqs, bam, genome, aligner="", genome_dir="", algn_cmd="", force=False):
#    """Map fastq reads to the genome and generate a bam file.
#    :param fqs List of fastq filenames.
#    :type  fqs list
#    :param bam Path to the bam file to generate.
#    :type  bam string
#    :param genome Genome name to map reads to (as defined in the configfile).
#    :type  genome string
#    :param aligner Alignment program (as defined in the configfile).
#    :type  aligner string
#    :param genome_dir Directory of the genome files for the mapper (as definced in configfile).
#    :type  genome_dir string
#    :param algn_cmd Command structure for the alignment wrapper.
#    :type  algn_cmd string
#    :param force Re-map even if bamfiles exist already. Defaults to False.
#    :type  force boolean
#    """
#    bams = [] 
#    
#    if os.path.exists(bam) and not force:
#        sys.stderr.write("{0} already exists, skipping...\n".format(bam))
#        return
#
#    for fq in fqs:
#        #first, map to intermediate bams
#        bname = fq.replace(".fq.gz", "") 
#        if os.path.exists(bname+".bam"):
#            os.unlink(bname+".bam")
#        
#        #delete temporary directory of the bwa mapper
#        if os.path.exists(bname):
#            shutil.rmtree(bname)
#        
#        cmd = algn_cmd.format(fq, bname, genome, aligner, genome_dir)
#        sys.stderr.write("Mapping {0} to {1}\n".format(fq, genome))
#        p = sp.Popen(cmd, 
#                             shell=True, 
#                             stderr=sp.PIPE, 
#                             stdout=sp.PIPE)
#        stdout, stderr = p.communicate() #bwa gives a lot of routine messages on stdout, so don't check for it
#        
#        if not os.path.exists("{0}.bam".format(bname)):
#            sys.stderr.write("Alignment failed\n")
#            sys.stderr.write(stderr)
#            sys.exit(1)
#        bams.append("{0}.bam".format(bname))
#    
#    if len(bams) == 1:
#        os.rename(bams[0], bam)
#        os.rename(bams[0] + ".bai", bam + ".bai")
#
#    elif len(bams) > 1:
#        cmd = "{0} merge -f {1} {2} && samtools index {1}"
#        cmd = cmd.format(
#                         "samtools",
#                         bam,
#                         " ".join(bams),
#                         )
#        ret = sp.call(cmd, shell=True)
#        if not ret:
#            for bamfile in bams:
#                os.unlink(bamfile)
#                os.unlink("{0}.bai".format(bamfile))
#        else:
#            sys.stderr.write("Merging failed!\n")
#            sys.exit(1)



def bam2bw(bam, bw, library=None, force=False):
    """ Convert bam to UCSC bigWig 
    :param bam Path to the bam file.
    :type  bam string
    :param bw Path to output bigwig file.
    :type  bw string
    :param library Library of the sample (e.g. 'RNA-seq' or 'Chip_seq')
    :type  library string
    """
    # Don't overwrite existing bigWig if not explicitly stated
    if os.path.exists(bw) and not force:
        sys.stderr.write("{0} already exists, skipping...\n".format(bw))
        return None, None
    
    if library == "RNA-Seq":
        # Do not extend, keep duplicates
        cmd = "bam2bw -i {0} -o {1} -D".format(bam, bw)
    elif library == "ChIP-Seq":
        # Extend to estimated fragment size and normalize
        cmd = "bam2bw -i {0} -o {1} -e auto -c".format(bam, bw)
    else:
        # Same as ChIP-seq for now
        cmd = "bam2bw -i {0} -o {1} -e auto -c".format(bam, bw)
    
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout,stderr = p.communicate()
    return stdout, stderr
