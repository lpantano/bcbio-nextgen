import os
import pandas as pd

import numpy as np

from bcbio.utils import (file_exists, tmpfile, chdir, splitext_plus,
                         max_command_length, robust_partition_all, append_stem)
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio import broad
from bcbio.pipeline import config_utils


def checkpoint(stem):
    def check_file(f):
        def wrapper(*args, **kwargs):
            out_file = append_stem(args[0], "_summary")
            if file_exists(out_file):
                logger.debug("Skipping %s" % out_file)
                return None
            return f(*args, **kwargs)
        return wrapper
    return check_file

@checkpoint("_summary")
def _calculate_percentiles(in_file, sample):
    """
    Parse pct bases per region to summarize it in
    7 different pct of regions points with pct bases covered
    higher than a completeness cutoff (5, 10, 20, 50 ...)
    """
    out_file = append_stem(in_file, "_summary")
    out_total_file = append_stem(in_file, "_total_summary")
    dt = pd.read_csv(in_file, sep="\t", index_col=False)
    pct = dict()
    pct_bases = dict()
    size = np.array(dt["chromEnd"]) - np.array(dt["chromStart"])
    for cutoff in [h for h in list(dt) if h.startswith("percentage")]:
        a = np.array(dt[cutoff])
        for p_point in [0.01, 10, 25, 50, 75, 90, 99.9]:
            q = np.percentile(a, p_point)
            pct[(cutoff, p_point)] = q
        pct_bases[cutoff] = sum(size * a)/float(sum(size))

    with file_transaction(out_total_file) as tx_file:
        with open(tx_file, 'w') as out_handle:
            print >>out_handle, "cutoff_reads\tbases_pct\tsample"
            for k in pct_bases:
                print >>out_handle, "\t".join(map(str, [k, pct_bases[k], sample]))
    with file_transaction(out_file) as tx_file:
        with open(tx_file, 'w') as out_handle:
            print >>out_handle, "cutoff_reads\tregion_pct\tbases_pct\tsample"
            for k in pct:
                print >>out_handle, "\t".join(map(str, [k[0], k[1], pct[k], sample]))

def coverage(data):
    """
    Calculate coverage at different completeness cutoff
    for region in coverage option.
    """
    bed_file = dd.get_coverage(data)
    if not bed_file:
        return data

    work_dir = os.path.join(dd.get_work_dir(data), "report", "coverage")
    with chdir(work_dir):
        in_bam = data['work_bam']
        sample = dd.get_sample_name(data)
        logger.debug("doing coverage for %s" % sample)
        parse_file = os.path.join(sample + "_coverage.bed")
        parse_total_file = os.path.join(sample + "_cov_total.tsv")
        cores = dd.get_num_cores(data)
        if not file_exists(parse_file):
            with file_transaction(parse_file) as out_tx:
                cmd = ("sambamba depth region -F \"not unmapped\" -t {cores} -C 200 -T 1 -T 5 -T 10 -T 20 -T 40 -T 50 -T 60 -T 70 -T 80 -T 100 -L {bed_file}  {in_bam} -o {parse_file}")
                do.run(cmd.format(**locals()), "Run coverage for {sample}")
        _calculate_percentiles(parse_file, sample)
        data['coverage'] = os.path.abspath(parse_file)
        return data

def variants(data):
    if not "vrn_file" in  data:
        return data
    in_vcf = data['vrn_file']
    work_dir = os.path.join(dd.get_work_dir(data), "report", "variants")
    with chdir(work_dir):
        in_bam = data['work_bam']
        ref_file = dd.get_ref_file(data)
        assert ref_file, "Need the reference genome fasta file."
        jvm_opts = broad.get_gatk_framework_opts(data['config'])
        gatk_jar = config_utils.get_program("gatk", data['config'], "dir")
        bed_file = dd.get_variant_regions(data)
        sample = dd.get_sample_name(data)
        in_bam = data["work_bam"]
        cg_file = os.path.join(sample + "_with-gc.vcf.gz")
        parse_file = os.path.join(sample + "_gc-depth-parse.tsv")
        num_cores = dd.get_num_cores(data)
        if not file_exists(cg_file):
            with file_transaction(cg_file) as tx_out:
                cmd = ("java -jar {gatk_jar}/GenomeAnalysisTK.jar -T VariantAnnotator "
                       "-R {ref_file} "
                       "-L {bed_file} -I {in_bam} "
                       "--num_threads {num_cores} "
                       "-A GCContent --variant {in_vcf} --out {tx_out}")
                do.run(cmd.format(**locals()), " GC bias for %s" % in_vcf)

        if not file_exists(parse_file):
            with file_transaction(parse_file) as out_tx:
                with open(out_tx, 'w') as out_handle:
                    print >>out_handle, "CG\tdepth\tsample"
                cmd = ("bcftools query -f '[%GC][\\t%DP][\\t%SAMPLE]\\n' -R "
                       "{bed_file} {cg_file} >> {out_tx}")
                do.run(cmd.format(**locals()),
                       "Calculating GC content and depth for %s" % in_vcf)
                logger.debug('parsing coverage: %s' % sample)
        return data
