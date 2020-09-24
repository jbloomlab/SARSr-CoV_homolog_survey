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
        build_variants=nb_markdown('build_variants.ipynb'),
        codon_variant_table=config['codon_variant_table_file'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        analyze_counts=nb_markdown('analyze_counts.ipynb'),
        compute_Kd='results/summary/compute_binding_Kd.md',
        Titeseq_Kds_file=config['Titeseq_Kds_file'],
        Titeseq_Kds_homologs_file=config['Titeseq_Kds_homologs_file'],
        compute_meanF='results/summary/compute_expression_meanF.md',
        expression_sortseq_file=config['expression_sortseq_file'],
        expression_sortseq_homologs_file=config['expression_sortseq_homologs_file'],
        global_epistasis_binding=nb_markdown('global_epistasis_binding.ipynb'),
        global_epistasis_expression=nb_markdown('global_epistasis_expression.ipynb'),
        single_mut_effects='results/summary/single_mut_effects.md',
        single_mut_effects_file=config['single_mut_effects_file'],
        homolog_effects_file=config['homolog_effects_file'],
        structure_function='results/summary/structure_function.md',
        logoplots_of_muteffects='results/summary/logoplots_of_muteffects.md',
        dms_view_file_RBD=config['dms_view_file_RBD'],
        dms_view_file_spike=config['dms_view_file_spike'],
        circulating_variants='results/summary/circulating_variants.md',
        antibody_epitopes='results/summary/antibody_epitopes.md',
        sarbecovirus_diversity='results/summary/sarbecovirus_diversity.md',
        interactive_heatmap='results/summary/interactive_heatmap.md',
        interactive_heatmap_html=config['interactive_heatmap'],
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

            1. [Process PacBio CCSs]({path(input.process_ccs)}).

            2. [Build variants from CCSs]({path(input.build_variants)}).
               Creates a [codon variant table]({path(input.codon_variant_table)})
               linking barcodes to the mutations in the variants.

            3. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            4. [QC analysis of sequencing counts]({path(input.analyze_counts)}).
            
            5. [Computation of ACE2-binding *K*<sub>D</sub>]({path(input.compute_Kd)}).
               Creates files giving the ACE2-binding of each barcoded variant
               [of SARS-CoV-2 RBD]({path(input.Titeseq_Kds_file)}) and of
               [the homologs]({path(input.Titeseq_Kds_homologs_file)}).
            
            6. [Computation of expression mean fluorescence]({path(input.compute_meanF)}).
               Creates files giving the expression of each barcoded variant
               [of SARS-CoV-2 RBD]({path(input.expression_sortseq_file)}) and of
               [the homologs]({path(input.expression_sortseq_homologs_file)}).
            
            7. [Global epistasis decomposition of binding effects]({path(input.global_epistasis_binding)}).
            
            8. [Global epistasis decomposition of expression effects]({path(input.global_epistasis_expression)}).
            
            9. [Calculation of final single mutant effects on binding and expression]({path(input.single_mut_effects)}).
               Creates files giving the estimated expression and ACE2-binding of
               [single mutants to SARS-CoV-2 RBD]({path(input.single_mut_effects_file)})
               and [the homologs]({path(input.homolog_effects_file)}).
               
            10. [Structure-function analysis of mutational effects]({path(input.structure_function)}).

            11. [Logo plots of mutational effects]({path(input.logoplots_of_muteffects)}).
                Also creates input files for `dms-view` of [RBD]({path(input.dms_view_file_RBD)}) and [spike]({path(input.dms_view_file_spike)}), the visualizations of which can be seen [here](https://jbloomlab.github.io/SARS-CoV-2-RBD_DMS/structures/).

            12. [Mutational constraint within RBD antibody epitopes]({path(input.antibody_epitopes)})

            13. [RBD variation across the sarbecovirus clade]({path(input.sarbecovirus_diversity)})
            
            14. [RBD variation in circulating SARS-CoV-2 isolates]({path(input.circulating_variants)}).

            15. [Make interactive heat map]({path(input.interactive_heatmap)}).
                Creates [this heatmap](https://jbloomlab.github.io/SARS-CoV-2-RBD_DMS/).

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

rule interactive_heatmap:
    input:
        config['single_mut_effects_file']
    output:
        nb_markdown=nb_markdown('interactive_heatmap.ipynb'),
        heatmap_html=config['interactive_heatmap']
    params:
        nb='interactive_heatmap.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule sarbecovirus_diversity:
    input:
        config['single_mut_effects_file'],
        config['homolog_effects_file']
    output:
        md='results/summary/sarbecovirus_diversity.md',
        md_files = directory('results/summary/sarbecovirus_diversity_files')
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='sarbecovirus_diversity.Rmd',
        md='sarbecovirus_diversity.md',
        md_files='sarbecovirus_diversity_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule antibody_epitopes:
    input:
        config['single_mut_effects_file'],
        config['homolog_effects_file']
    output:
        md='results/summary/antibody_epitopes.md',
        md_files = directory('results/summary/antibody_epitopes_files')
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='antibody_epitopes.Rmd',
        md='antibody_epitopes.md',
        md_files='antibody_epitopes_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule circulating_variants:
    input:
        config['single_mut_effects_file']
    output:
        md='results/summary/circulating_variants.md',
        md_files = directory('results/summary/circulating_variants_files')
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='circulating_variants.Rmd',
        md='circulating_variants.md',
        md_files='circulating_variants_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule logoplots_of_muteffects:
    input:
        config['single_mut_effects_file']
    output:
        nb_markdown=nb_markdown('logoplots_of_muteffects.ipynb'),
        dms_view_file_RBD=config['dms_view_file_RBD'],
        dms_view_file_spike=config['dms_view_file_spike']
    params:
        nb='logoplots_of_muteffects.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule structure_function:
    input:
        config['single_mut_effects_file'],
        config['homolog_effects_file']
    output:
        md='results/summary/structure_function.md',
        md_files = directory('results/summary/structure_function_files')
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='structure_function.Rmd',
        md='structure_function.md',
        md_files='structure_function_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule single_mut_effects:
    input:
        config['global_epistasis_binding_file'],
        config['global_epistasis_expr_file'],
        config['Titeseq_Kds_homologs_file'],
        config['expression_sortseq_homologs_file'],
    output:
        config['single_mut_effects_file'],
        config['homolog_effects_file'],
        md='results/summary/single_mut_effects.md',
        md_files = directory('results/summary/single_mut_effects_files')
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='single_mut_effects.Rmd',
        md='single_mut_effects.md',
        md_files='single_mut_effects_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule global_epistasis_binding:
    input:
        config['Titeseq_Kds_file']
    output:
        config['global_epistasis_binding_file'],
        nb_markdown=nb_markdown('global_epistasis_binding.ipynb')
    params:
        nb='global_epistasis_binding.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule global_epistasis_expression:
    input:
        config['expression_sortseq_file']
    output:
        config['global_epistasis_expr_file'],
        nb_markdown=nb_markdown('global_epistasis_expression.ipynb')
    params:
        nb='global_epistasis_expression.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule compute_Titeseq_Kds:
    input:
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file'],
        config['Titeseq_Kds_homologs_file'],
        md='results/summary/compute_binding_Kd.md',
        md_files=directory('results/summary/compute_binding_Kd_files')
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='compute_binding_Kd.Rmd',
        md='compute_binding_Kd.md',
        md_files='compute_binding_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule compute_expression_meanFs:
    input:
        config['variant_counts_file']
    output:
        config['expression_sortseq_file'],
        config['expression_sortseq_homologs_file'],
        md='results/summary/compute_expression_meanF.md',
        md_files=directory('results/summary/compute_expression_meanF_files')
    envmodules:
        'R/3.6.1-foss-2018b'
    params:
        nb='compute_expression_meanF.Rmd',
        md='compute_expression_meanF.md',
        md_files='compute_expression_meanF_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule analyze_counts:
    """Analyze variant counts and compute functional scores."""
    input:
        config['variant_counts_file']
    output:
        nb_markdown=nb_markdown('analyze_counts.ipynb')
    params:
        nb='analyze_counts.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['codon_variant_table_file'],
        config['barcode_runs']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule build_variants:
    """Build variant table from processed CCSs."""
    input:
        config['processed_ccs_file']
    output:
        config['codon_variant_table_file'],
        nb_markdown=nb_markdown('build_variants.ipynb')
    params:
        nb='build_variants.ipynb'
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
