"""
Run mirdeep from bam file
"""
import os.path as op

import pysam
from seqcluster.libs import inputs

def _get_format(name, info):
    counts = info[int(name.replace("seq_",""))].total
    return name + "_x%" % counts

def _convert_bam_file(bam_in, ma_file):
    """
    Replace sequences with correct names
    """
    seqs = inputs.parse_ma_file_raw(ma_file)
    bam = pysam.AlignmentFile(bam_in, "rb")
    out_file = op.splitext(bam_in) + ".sam"
    with pysam.AlignmentFile(out_file, "w", template=bam) as out_handle:
        for read in bam.fetch():
            read.query_name = _get_format(read.query_name, seqs)
            out_handle.write(read)
    return out_file

def _run(bam_file, out_dir):
    cmd = ("miRDeep2.pl")

def mirdeep2(data):
    """
    Run mirdeep prediction with all samples
    """
    bam_file = op.join(dd.get_work_dir(data[0]), "align", "seqs.bam")
    out_dir = os.path.join(dd.get_work_dir(data[0]), "mirdeep")
    out_dir = os.path.abspath(safe_makedir(out_dir))
    ma_file = os.path.join(dd.get_work_dir(data[0]), "seqcluster", "prepare", "seqs.ma")
    bam_file = _convert_bam_file(bam_file, ma_file)
    _run(bam_file, out_dir)
