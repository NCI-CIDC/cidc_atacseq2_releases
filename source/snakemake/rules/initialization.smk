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

## Retrieve hg38 blacklist from https://github.com/Boyle-Lab/Blacklist
rule retrieve_hg38_blackist:
    output:
        'blacklist/hg38-blacklist.v2.bed'
    threads: 1
    shell:
        'wget -qO - https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz | gunzip > {output}' 

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
