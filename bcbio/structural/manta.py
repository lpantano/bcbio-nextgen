"""Structural variant detection with the Manta caller from Illumina.

https://github.com/Illumina/manta
"""
import os
import sys

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.provenance import do, programs

def run(items):
    """Perform detection of structural variations with Manta.
    """
    paired = vcfutils.get_paired(items)
    work_dir = _sv_workdir(paired.tumor_data if paired else items[0])
    workflow_file = _prep_config(items, paired, work_dir)
    variant_file = _run_workflow(items, paired, workflow_file, work_dir)
    out = []
    for data in items:
        sample_file = _select_sample(data, variant_file, work_dir)
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append({"variantcaller": "manta",
                           "vrn_file": sample_file})
        out.append(data)
    return out

def _run_workflow(items, paired, workflow_file, work_dir):
    """Run manta analysis inside prepared workflow directory.
    """
    data = paired.tumor_data if paired else items[0]
    if paired:
        if paired.normal_bam:
            base_file = "somaticSV.vcf.gz"
        else:
            base_file = "tumorSV.vcf.gz"
    else:
        base_file = "diploidSV.vcf.gz"
    out_file = os.path.join(work_dir, "results", "variants", base_file)
    if not utils.file_exists(out_file):
        utils.remove_safe(os.path.join(work_dir, "workspace"))
        cmd = [sys.executable, workflow_file, "-m", "local", "-j", dd.get_num_cores(data),
               "--quiet"]
        do.run(cmd, "Run manta SV analysis")
    utils.remove_safe(os.path.join(work_dir, "workspace"))
    return out_file

def _select_sample(data, variant_file, work_dir):
    """Select current sample from original call file.
    """
    sample_name = dd.get_sample_name(data)
    out_file = os.path.join(work_dir, "%s-%s.vcf.gz" % (utils.splitext_plus(os.path.basename(variant_file))[0],
                                                        sample_name))
    if not utils.file_uptodate(out_file, variant_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = "bcftools view -s {sample_name} -O z -o {tx_out_file} {variant_file}"
            do.run(cmd.format(**locals()), "Run manta SV analysis")
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _prep_config(items, paired, work_dir):
    """Run initial configuration, generating a run directory for Manta.
    """
    assert utils.which("configManta.py"), "Could not find installed configManta.py"
    out_file = os.path.join(work_dir, "runWorkflow.py")
    if not utils.file_exists(out_file) or _out_of_date(out_file):
        cmd = [sys.executable, utils.which("configManta.py")]
        if paired:
            if paired.normal_bam:
                cmd += ["--normalBam=%s" % paired.normal_bam, "--tumorBam=%s" % paired.tumor_bam]
            else:
                cmd += ["--tumorBam=%s" % paired.tumor_bam]
        else:
            cmd += ["--bam=%s" % dd.get_align_bam(data) for data in items]
        data = paired.tumor_data if paired else items[0]
        cmd += ["--referenceFasta=%s" % dd.get_ref_file(data), "--runDir=%s" % work_dir]
        if dd.get_coverage_interval(data) not in ["genome"]:
            cmd += ["--exome"]
        do.run(cmd, "Configure manta SV analysis")
    return out_file

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "manta"))

def _out_of_date(rw_file):
    """Check if a run workflow file points to an older version of manta and needs a refresh.
    """
    with open(rw_file) as in_handle:
        for line in in_handle:
            if line.startswith("sys.path.append"):
                file_version = line.split("/lib/python")[0].split("Cellar/manta/")[-1]
                if file_version != programs.get_version_manifest("manta"):
                    return True
    return False