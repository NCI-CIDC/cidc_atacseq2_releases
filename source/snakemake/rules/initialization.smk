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

## Retrieve hg38 blacklist from https://github.com/Boyle-Lab/Blacklist
rule retrieve_hg38_blackist:
    output:
        'blacklist/hg38-blacklist.v2.bed'
    threads: 1
    shell:
        'wget -qO - https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz | gunzip > {output}' 

## Download built bwa_index files for the specified genome
## If using different genome, need to edit rule to build using 'bwa index'
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
          echo "Downloading bwa_index files from ncbi ftp associated with genome for mapping reads to GRCh38 or hg38..." | tee {log}
          echo "wget \"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz\"" | tee -a {log}
          wget -nv "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz"; 2>> {log}
          tar -xzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz -C genome
          rm GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
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
