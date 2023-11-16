## Set up directory structure based on dirs supplied in config
## Ignores non-zero exit status returned when any directories already exist
rule directory_setup:
    output:
        'progress/dirs.done'
    params:
        subdirs=SUBDIRS
    threads:1
    shell:
        '''
          mkdir {params.subdirs} -p 2> /dev/null
          touch {output}
        '''

## Build reference genome index
rule build_bwa_index:
    input:
        rules.directory_setup.output
    output:
        'progress/bwa_index_built.done'
    benchmark:
        'benchmark/build_bwa_index.tab'
    log:
        'log/build_bwa_index.log'
    conda:
        SOURCEDIR+"/../envs/bwa.yaml"
    params:
        indexseq=paths.genome.fa
    priority: 1000
    threads: max(1,NCORES)
    shell:
        '''
          echo "bwa index {params.indexseq}" > {log}
          bwa index {params.indexseq} 2>> {log}
          touch {output}

          ## export rule env details
          conda env export --no-builds > info/bwa.info
        '''

## Get genome chrom sizes for bw generation
rule genome_size:
    input:
        rules.build_bwa_index.output
    output:
        paths.genome.size
    benchmark:
        'benchmark/genome_size.tab'
    log:
        'log/genome_size.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    params:
        indexseq=paths.genome.fa
    priority: 2
    threads: 1
    shell:
        '''
          ## get genome chrom size
          echo "samtools faidx {params.indexseq}" | tee {log}
          samtools faidx {params.indexseq} 2>> {log}
          echo "cut -f1,2 {params.indexseq}.fai > {output}" | tee -a {log}
          cut -f1,2 {params.indexseq}.fai > {output} 2>> {log}
          
          ## export rule env details
          conda env export --no-builds > info/samtools.info
        '''