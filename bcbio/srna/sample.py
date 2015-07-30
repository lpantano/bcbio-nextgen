import os
import os.path as op
import shutil
from collections import Counter
from seqcluster.libs.fastq import collapse, write_output

from bcbio.utils import (splitext_plus, file_exists, append_stem, replace_directory)
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.install import _get_data_dir

# ADAPTER = "AGATCGGAAGAGCAC"


def trim_srna_sample(data):
    adapter = dd.get_adapters(data)[0]
    names = data["rgnames"]['sample']
    out_dir = os.path.join(dd.get_work_dir(data), names, "trimmed")
    in_file = data["files"][0]
    utils.safe_makedir(out_dir)
    out_file = replace_directory(append_stem(in_file, ".clean"), out_dir)
    out_noadapter_file = replace_directory(append_stem(in_file, ".fragments"), out_dir)
    out_short_file = replace_directory(append_stem(in_file, ".short"), out_dir)
    cmd = _cmd_cutadapt()
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            do.run(cmd.format(**locals()), "remove adapter")
    data["clean_fastq"] = out_file
    data["collapse"] = _collapse(data["clean_fastq"])
    data["size_stats"] = _summary(data['collapse'])
    return [[data]]

def mirbase(data):
    names = data["rgnames"]['sample']
    out_dir = os.path.join(dd.get_work_dir(data), names, "mirbase")
    utils.safe_makedir(out_dir)
    out_file = op.join(out_dir, names)
    mirbase = op.abspath(op.join(dd.get_galaxy_dir(data), '..', "genomes", "mirbase"))
    data['seqbuster'] = _miraligner(data["collapse"], out_file, dd.get_species(data), mirbase, data['config'])
    return [[data]]

def _cmd_cutadapt():
    cmd = "cutadapt --adapter={adapter} --minimum-length=8 --untrimmed-output={out_noadapter_file} -o {tx_out_file} -m 17 --overlap=8 {in_file} --too-short-output {out_short_file}"
    return cmd

def _collapse(in_file):
    seqs = collapse(in_file)
    out_file = append_stem(in_file, ".trimming").replace(".gz", "")
    write_output(out_file, seqs)
    return out_file

def _summary(in_file):
    data = Counter()
    out_file = in_file + "_size_stats"
    if file_exists(out_file):
        return out_file
    with open(in_file) as in_handle:
        for line in in_handle:
            counts = int(line.strip().split("_x")[1])
            line = in_handle.next()
            l = len(line.strip())
            in_handle.next()
            in_handle.next()
            data[l] += counts
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, 'w') as out_handle:
            for l, c in data.iteritems():
                out_handle.write("%s %s\n" % (l, c))
    return out_file

def _miraligner(fastq_file, out_file, species, db_folder, config):
    resources = config_utils.get_resources("miraligner", config)
    jvm_opts =  "-Xms750m -Xmx4g"
    if resources and resources.get("jvm_opts"):
        jvm_opts = " ".join(resources.get("jvm_opts"))

    cmd = ("miraligner {jvm_opts} -freq -sub 1 -trim 3 -add 3 -s {species} -i {fastq_file} -db {db_folder}  -o {tx_out_file}")
    if not file_exists(out_file + ".mirna"):
        with file_transaction(out_file) as tx_out_file:
            do.run(cmd.format(**locals()), "Do miRNA annotation")
            shutil.move(tx_out_file + ".mirna", out_file + ".mirna")
    return out_file + ".mirna"
