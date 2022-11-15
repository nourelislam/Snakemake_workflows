
rule gunzip:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT'] + ".gz"
    output:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    shell:
        """
        gunzip {input}
        """


rule gzip:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch(config['genomeDIR'] + "done__gzip_{GE}")
    shell:
        """
        gzip {input}
        """

rule all_gzip:
    input:
        expand(config['genomeDIR'] + "{GE}" + config['genomeEXT'], GE=config['genomes'])

rule assembly_stats:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        "01_assembly_stats/{GE}.stats"
    params:
        "-t"
    conda:
        "../envs/assembly_stats.yaml"
    shell:
        """
        assembly-stats {params} {input} > {output}
        """


rule quast:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        directory("01_quast/{GE}")
    log:
        "01_quast/log__quast_{GE}"
    params:
        "01_quast"
    threads:
        config['threads_quast']
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast -o {output} -t {threads} {input} > {log} 2>&1
        """


rule busco:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_busco/done__busco_{GE}")
    log:
        os.path.join(workflow.basedir, "01_busco/log__busco_{GE}")
    params:
        args = "-m genome -o busco -l metazoa",
        outdir = "01_busco"
    threads:
        config['threads_busco']
    conda:
        "../envs/busco.yaml"
    shell:
        """
        d={params.outdir}/{wildcards.GE}
        mkdir -p $d && 
        cd $d && 
        busco -i ../../{input.fas} -c {threads} {params.args} > {log} 2>&1
        """

rule EDTA:
    input:
        fasta = config['genomeDIR'] + '{GE}' + config['genomeEXT']
    output:
        directory("01_EDTA/{GE}")
        #"01_EDTA/{GE}"
    params:
        arg = "--species others --step all --sensitive 1 -anno 1 --evaluate 1"
    conda:
        "../envs/EDTA.yaml"
    threads:
        config['threads_EDTA']
    shell:
        """
        EDTA.pl --genome {input.fasta} -t {threads} {params.arg}
        """