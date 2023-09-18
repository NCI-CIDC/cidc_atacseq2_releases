## Set up directory structure UPDATE make dirs based on config input
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
        "./envs/bwa.yaml"
    params:
        indexseq=INDEXSEQ
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