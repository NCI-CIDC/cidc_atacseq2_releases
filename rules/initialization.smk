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

## Download reference genome fa and supporting gtf annotation from ncbi ftp site
rule retrieve_reference_genome:
    input:
        rules.directory_setup.output
    output:
        fa=paths.genome.fa,
        gtf=paths.annot.gtf
    benchmark:
        'benchmark/retrieve_reference_genome.tab'
    log:
        'log/retrieve_reference_genome.log'
    params:
        fa_ftp=GENOME_FA_FTP,
        gtf_ftp=GENOME_GTF_FTP
    priority: 1000
    threads: 1
    shell:
        '''
          echo "downloading Genome to map reads to GRCh38 or hg38..." | tee {log}
          echo "wget -nv {params.fa_ftp} -O {output.fa}.gz" | tee -a {log}
          wget -nv {params.fa_ftp} -O {output.fa}.gz 2>> {log}
          gunzip -f {output.fa}.gz
          
          echo "downloading supporting GTF annotations..." | tee -a {log}
          echo "wget -nv {params.gtf_ftp} -O {output.gtf}.gz" | tee -a {log}
          wget -nv {params.gtf_ftp} -O {output.gtf}.gz 2>> {log}
          gunzip -f {output.gtf}.gz
        '''

## Download built bwa_index files for the specified genome
## If using different genome, need to edit rule to build using 'bwa index'
rule build_bwa_index:
    input:
        rules.directory_setup.output,
        rules.retrieve_reference_genome.output.fa
    output:
        'progress/bwa_index_built.done'
    benchmark:
        'benchmark/build_bwa_index.tab'
    log:
        'log/build_bwa_index.log'
    conda:
        SOURCEDIR+"/../envs/bwa.yaml"
    params:
        bwa_ftp=GENOME_BWA_FTP
    priority: 1000
    threads: 1
    shell:
        '''
          echo "Downloading bwa_index files from ncbi ftp associated with genome for mapping reads to GRCh38 or hg38..." | tee {log}
          echo "wget -nv {params.bwa_ftp}" | tee -a {log}
          wget -nv {params.bwa_ftp} 2>> {log}
          tar -xzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz -C genome
          rm GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
          touch {output}
          
          ## export rule env details
          conda env export --no-builds > info/bwa.info
        '''

## Get genome chrom sizes for bw generation
rule genome_size:
    input:
        genome_fa=rules.retrieve_reference_genome.output.fa
    output:
        size=paths.genome.size,
        fai=paths.genome.fai
    benchmark:
        'benchmark/genome_size.tab'
    log:
        'log/genome_size.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    threads: 1
    shell:
        '''
          ## get genome chrom size
          echo "samtools faidx {input.genome_fa}" | tee {log}
          samtools faidx {input.genome_fa} 2>> {log}
          echo "cut -f1,2 {input.genome_fa}.fai > {output.size}" | tee -a {log}
          cut -f1,2 {input.genome_fa}.fai > {output.size} 2>> {log}
          
          ## export rule env details
          conda env export --no-builds > info/samtools.info
        '''

# Filter the hg38 genome index and convert from fai to bed format 
rule create_bed:
    input:
        paths.genome.fai
    output:
        paths.genome.bed
    benchmark:
        'benchmark/create_bed.tab'
    threads: 1
    shell:
        '''
          ## Remove the entries from chrM, chrUN, _random, chrEBV in the hg38 genome index and convert fai format to bed
          grep -v -E 'chrUn|_random|chrEBV|chrM' {input} | awk -F'\t' '{{ print $1,\"0\",$2 }}' > {output}
        '''  

## Retrieve hg38 blacklist from https://github.com/Boyle-Lab/Blacklist
rule retrieve_hg38_blacklist:
    output:
        paths.genome.blacklist
    benchmark:
        'benchmark/retrieve_hg38_blacklist.tab'
    params:
        blacklist_url=GENOME_BLACKLIST_URL
    threads: 1
    shell:
        'wget -qO - {params.blacklist_url} | gunzip > {output}'

## Retrieve DHS regions list from dev GCP bucket. This might not be final location of the file.
## If file location changes, the shell directive needs to be updated.
rule retrieve_hg38_dhs:
    output:
        paths.genome.dhs
    benchmark:
        'benchmark/retrieve_hg38_dhs.tab'
    params:
        dhs_gc=GENOME_DHS_GC
    threads: 1
    shell:
        "gsutil cp {params.dhs_gc} {output}"
