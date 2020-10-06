"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import os.path
import textwrap

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config['pacbio_runs'], dtype = str)
               .assign(pacbioRun=lambda x: x['library'] + '_' + x['run'])
               )
assert len(pacbio_runs['pacbioRun'].unique()) == len(pacbio_runs['pacbioRun'])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        process_ccs=nb_markdown('process_ccs.ipynb'),
        nt_variant_table=config['nt_variant_table_file'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        #merge_sequencing='results/summary/merge_sequencing.Rmd.md',
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:

            1. [Process PacBio CCSs]({path(input.process_ccs)}). Creates a [table]({path(input.nt_variant_table)})
               linking barcodes to the mutations in the variants.

            2. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            3. [Parse amino acid mutants and merge PacBio and Illumina sequencing data]().

            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"

rule merge_sequencing:
    input:
        config['nt_variant_table_file'],
        config['variant_counts_file']
    output:
        md='results/summary/merge_sequencing.md',
        md_files=directory('results/summary/merge_sequencing_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='merge_sequencing.Rmd',
        md='merge_sequencing.md',
        md_files='merge_sequencing_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['nt_variant_table_file'],
        config['barcode_runs']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule process_ccs:
    """Process the PacBio CCSs."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
        config['amplicons'],
        ([] if config['seqdata_source'] != 'HutchServer' else
         expand(os.path.join(config['ccs_dir'], "{pacbioRun}_report.txt"),
                pacbioRun=pacbio_runs['pacbioRun'])
         )
    output:
        config['processed_ccs_file'],
        config['nt_variant_table_file'],
        nb_markdown=nb_markdown('process_ccs.ipynb')
    params:
        nb='process_ccs.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

if config['seqdata_source'] == 'HutchServer':

    rule build_ccs:
        """Run PacBio ``ccs`` program to build CCSs from subreads."""
        input:
            subreads=lambda wildcards: (pacbio_runs
                                        .set_index('pacbioRun')
                                        .at[wildcards.pacbioRun, 'subreads']
                                        )
        output:
            ccs_report=os.path.join(config['ccs_dir'], "{pacbioRun}_report.txt"),
            ccs_fastq=os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz")
        params:
            min_ccs_length=config['min_ccs_length'],
            max_ccs_length=config['max_ccs_length'],
            min_ccs_passes=config['min_ccs_passes'],
            min_ccs_accuracy=config['min_ccs_accuracy']
        threads: config['max_cpus']
        shell:
            """
            ccs \
                --min-length {params.min_ccs_length} \
                --max-length {params.max_ccs_length} \
                --min-passes {params.min_ccs_passes} \
                --min-rq {params.min_ccs_accuracy} \
                --report-file {output.ccs_report} \
                --num-threads {threads} \
                {input.subreads} \
                {output.ccs_fastq}
            """

elif config['seqdata_source'] == 'SRA':
    raise RuntimeError('getting sequence data from SRA not yet implemented')

else:
    raise ValueError(f"invalid `seqdata_source` {config['seqdata_source']}")
