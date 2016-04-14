"""
Some reports for the analysis
"""
import shutil
import os
import sys
import copy
from contextlib import closing
from collections import Counter
import pandas as pd
import math

import pysam

from bcbio.utils import splitext_plus, file_exists, safe_makedir, chdir
from bcbio.log import logger
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import bam
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd

def _run_meth_extractor(bam_in, sample, workdir, config):
    """
    Run bismark_methylation_extractor command
    """
    bismark = config_utils.get_program("bismark_methylation_extractor", config)
    cores = config["algorithm"].get('cores', 1)
    memory = config["algorithm"].get('mem', 5)
    bam_in = bam.sort(bam_in, config, order="queryname")
    cmd = "{bismark}  --no_overlap --comprehensive --multicore {cores} --buffer_size {memory}G --bedGraph --counts --gzip {bam_in}"
    out_dir = os.path.join(workdir, sample)
    mbias_file = os.path.join(out_dir, os.path.basename(splitext_plus(bam_in)[0]) + '.M-bias.txt')
    if not file_exists(mbias_file):
        with tx_tmpdir() as tx_dir:
            with chdir(tx_dir):
                do.run(cmd.format(**locals()), "bismark_methylation_extractor  in %s" % bam_in)
                shutil.move(tx_dir, out_dir)
    assert os.path.exists(mbias_file), "mbias report doesn't exists:%s" % mbias_file
    return mbias_file

def _run_report(bam_in, bam_report, sample, biasm_file, workdir, config):
    """
    Run bismark2report command
    """
    bismark = config_utils.get_program("bismark2report", config)
    cmd = "{bismark} --alignment_report {bam_report} -o {tx_out} --mbias_report {biasm_file}"
    out_dir = os.path.join(workdir, sample)
    out_file = os.path.join(out_dir, sample + '.html')
    with chdir(out_dir):
        if not file_exists(out_file):
            with file_transaction(out_file) as tx_out:
                do.run(cmd.format(**locals()), "bismarkr2report  in %s" % bam_in)
    return out_dir

def _bismark_calling(data):
    workdir = safe_makedir(os.path.join(dd.get_work_dir(data), "cpg"))
    config = data["config"]
    sample = dd.get_sample_name(data)
    biasm_file = _run_meth_extractor(data["work_bam"], sample, workdir, config)
    data['bismark_report'] = _run_report(data["work_bam"], data["bam_report"], sample, biasm_file, workdir, config)
    return data

def _bsmap_calling(data):
    sample = dd.get_sample_name(data)
    workdir = safe_makedir(os.path.join(dd.get_work_dir(data), "cpg_split", sample))
    config = data["config"]
    ref = dd.get_sam_ref(data)
    work_bam = dd.get_work_bam(data)
    python = os.path.join(os.path.dirname(sys.executable), "python")
    methratio = config_utils.get_program("methratio.py", config)
    cmd = ("{python} {methratio} -n -u -p -r -m 2 --chr={chrom} --ref={ref} {work_bam} >> {out_tx}")
    chrom = data["chr_to_run"]

    out_file = os.path.join(workdir, "methyratios_%s.txt" % chrom)
    if not file_exists(out_file):
        with file_transaction(out_file) as out_tx:
            do.run(cmd.format(**locals()), "Extract methylation for: %s" % sample)
    data["cpg_file"] = out_file
    return data

def calling(data):
    if dd.get_aligner(data) == "bismark":
        data = _bismark_calling(data)
    if dd.get_aligner(data) == "bsmap":
        data = _bsmap_calling(data)
    return [[data]]

def parallel_calling(data, run_parallel):
    out = []
    for sample in data:
        work_bam = dd.get_work_bam(sample[0])
        with closing(pysam.Samfile(work_bam, "rb")) as pysam_work_bam:
             chroms = pysam_work_bam.references
             for chrom in chroms:
             # for chrom in ["chr10"]:
                 new_sample = copy.deepcopy(sample)
                 if chrom.find("_") > -1:
                     continue
                 new_sample[0]['chr_to_run'] = chrom
                 out.append(new_sample)
    out = run_parallel("cpgcalling", out)
    for sample in out:
        phenotype = dd.get_phenotype(sample[0])
        batch = dd.get_batch(sample[0])
        if phenotype == "mC":
            for sample2 in out:
                if batch in dd.get_batch(sample2[0]) and dd.get_phenotype(sample2[0]) == "hmC":
                    if sample[0]["chr_to_run"] == sample2[0]["chr_to_run"]:
                        sample[0]["control"] = sample2[0]["cpg_file"]
                        break
    run_parallel("cpg_processing", out)
    # run_parallel("cpg_processing", data)
    #_make_stats(data)

def _sync_pos(handle, tag):
    for line in handle:
        cols = line.strip().split("\t")
        if cols[3] != "CG":
            continue
        pos = int(cols[1])
        # print "read another line of 2 file: %s" % line.strip()
        if tag <= pos:
            return [pos, float(cols[4])]
    return [None, None]

def cpg_postprocessing(data):
    mC = data["cpg_file"]
    if not "control" in data:
        return [[data]]
    hmC = data["control"]
    pos_hmC = 0
    pos_mC = 0
    counts = 0
    hmC_counts = 0
    with open(mC) as mC_h:
        with open(hmC) as hmC_h:
            for line in mC_h:
                cols = line.strip().split("\t")
                if cols[3] != "CG":
                    continue
                # print "read line of 1 file %s" % line.strip()
                pos = int(cols[1])
                ratio = float(cols[4])
                if pos < pos_hmC:
                    continue
                elif pos > pos_hmC:
                    pos_hmC, ratio_hmC = _sync_pos(hmC_h, pos)
                if not pos_hmC:
                    break
                if pos == pos_hmC:
                    counts += 1
                    if ratio > 0.2 and ratio_hmC < 0.2:
                        hmC_counts += 1
                        # print "It is a hmC site %s ratio in hmC %s " % (line.strip(), ratio_hmC)
                    elif ratio > 0.2 and ratio_hmC > 0.2:
                        continue
                        # print "It is a mC site %s ratio in hmC %s " % (line.strip(), ratio_hmC)
    print "%s %s counts %s over %s" % (mC_h, hmC_h, hmC_counts, counts)

def make_stats(sample):
    dtdepth = Counter()
    dtratio = Counter()
    work_dir = dd.get_work_dir(sample)
    sample_name = dd.get_sample_name(sample)
    depth_out = os.path.join(work_dir, "cpg_split", sample_name, "depth.tsv")
    ratio_out = os.path.join(work_dir, "cpg_split", sample_name, "ratio.tsv")
    work_bam = dd.get_work_bam(sample)
    with closing(pysam.Samfile(work_bam, "rb")) as pysam_work_bam:
         chroms = pysam_work_bam.references
         for chrom in chroms:
             logger.debug("Reading %s of sample %s" % (chrom, sample_name))
             cpg_file = os.path.join(work_dir, "cpg_split", sample_name, "methyratios_%s.txt" % chrom)
             if file_exists(cpg_file):
                 with open(cpg_file) as in_handle:
                     for line in in_handle:
                         cols = line.strip().split("\t")
                         if cols[3] == "CG":
                             ratio = int(float(cols[4]) * 100)
                             dtratio[ratio] += 1
                             depth = int(math.ceil(float(cols[5]))) if float(cols[5]) < 10 else 10
                             dtdepth[depth] += 1
         pd.DataFrame(dtdepth, index=[1]).to_csv(depth_out, sep="\t")
         pd.DataFrame(dtratio, index=[1]).to_csv(ratio_out, sep="\t")

