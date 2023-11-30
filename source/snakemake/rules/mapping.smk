## Map reads to the reference genome using BWA and output sorted bam
rule run_bwa:
    input:
        tch=rules.build_bwa_index.output,
        fa=rules.qualityfilter.output
    output:
        paths.bam.bam
    benchmark:
        'benchmark/{sample}_run_bwa.tab'
    log:
        'log/{sample}_run_bwa.log'
    conda:
        SOURCEDIR+"/../envs/bwa.yaml"
    params:
        sample='{sample}',
        indexseq=paths.genome.fa,
        in_fa_str=expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[0] + ' ' + expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[2] if len(ENDS) == 2 else expand(paths.rqual_filter.qfilter_fastq_single, read=ENDS)[0],
        srcdir=SOURCEDIR,
        doencrypt=DOENCRYPT,
        openssl=OPENSSL,
        encrypt_pass=ENCRYPT_PASS,
        hash=HASH,
        doarchive=DOARCHIVE,
        archive=ARCHIVE,
        cloud=CLOUD
    priority: 4
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "bwa mem -t {threads} {params.indexseq} {params.in_fa_str} | samtools view -@ {threads} -Sbh | samtools sort -@ {threads} > {output}" | tee {log}
          bwa mem -t {threads} {params.indexseq} {params.in_fa_str} | samtools view -@ {threads} -Sbh | samtools sort -@ {threads} > {output} 2>> {log}

          ## encrypt and archive if needed
          python3 {params.srcdir}/python/init-encrypt-archive.py --src {params.srcdir}\
          --doencrypt {params.doencrypt} --openssl {params.openssl} --encrypt-pass {params.encrypt_pass} --hash {params.hash} \
          --doarchive {params.doarchive} --archive {params.archive} --cloud {params.cloud} \
          --output {output} \
          2>> {log}
        '''

## Perform post-alignment filtering on the sorted bam
rule filter_bam:
    input:
        bam=rules.run_bwa.output,
        blacklist=rules.retrieve_hg38_blacklist.output,
        bed=rules.create_bed.output
    output:
        dup_bam=paths.bam.dup_bam,
        metrics=paths.bam.metrics,
        temp_bam=paths.bam.temp_bam,
        filtered_bam=paths.bam.filtered_bam,
        index=paths.bam.filtered_index
    benchmark:
        'benchmark/{sample}_filer_bam.tab'
    conda:
        SOURCEDIR+"/../envs/filter_bam.yaml"
    shell:
        '''
          ## Mark duplicates and provide duplicate stats for the raw bam
          ## NOTE: The command line for Picard is going to be updated in the future. Refer to link below.
          ## https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
          picard MarkDuplicates I={input.bam} O={output.dup_bam} M={output.metrics}

          ## Index the bam with duplicates marked
          samtools index {output.dup_bam}

          ## Remove reads that are unmapped, mate unmapped (for paired-end), not primary alignment,
          ## fail platform/vendor quality checks, and PCR or optical duplicates. Also filters
          ## reads aligned to chrM, chrUN, _random, chrEBV.
          samtools view -b -F 1804 -L {input.bed} {output.dup_bam} -o {output.temp_bam}

          ## Remove alignments that are located in the hg38 blacklist
          bedtools intersect -v -abam {output.temp_bam} -b {input.blacklist} > {output.filtered_bam}
         
          ## Index the final filtered bam
          samtools index {output.filtered_bam} -o {output.index}
        '''

## Index BAM
rule index_bam:
    input:
        bam=rules.run_bwa.output
    output:
        paths.bam.index
    benchmark:
        'benchmark/{sample}_index_bam.tab'
    log:
        'log/{sample}_index_bam.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    params:
        sample='{sample}'
    priority: 5
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "samtools index -@ {threads} {input.bam}" > {log}
          samtools index -@ {threads} {input.bam} 2>> {log}
        '''

## Run FASTQC
rule fastqc:
    input:
        rules.run_bwa.output
    output:
        paths.fastqc.targz
    benchmark:
        'benchmark/{sample}_fastqc.tab'
    log:
        'log/{sample}_fastqc.log'
    conda:
        SOURCEDIR+"/../envs/fastqc.yaml"
    params:
        sample='{sample}',
        fq_base='fastqc/{sample}_fastqc',
        fq_zip='fastqc/{sample}_fastqc.zip',
        fq_html='fastqc/{sample}_fastqc.html',
        srcdir=SOURCEDIR,
        doencrypt=DOENCRYPT,
        openssl=OPENSSL,
        encrypt_pass=ENCRYPT_PASS,
        hash=HASH,
        doarchive=DOARCHIVE,
        archive=ARCHIVE,
        cloud=CLOUD
    priority: 1
    threads: 1
    shell:
        '''
          echo "fastqc {input} -q -o fastqc" > {log}
          fastqc {input} -q -o fastqc 2>> {log}

          ## unzip, remove zipped results, HTML duplicate, and tarball results
          unzip -qq {params.fq_zip} -d {params.fq_base} && tar -zcf {output} {params.fq_base} && rm -r {params.fq_zip} {params.fq_html} {params.fq_base}

          ## encrypt and archive if needed
          python3 {params.srcdir}/python/init-encrypt-archive.py --src {params.srcdir}\
          --doencrypt {params.doencrypt} --openssl {params.openssl} --encrypt-pass {params.encrypt_pass} --hash {params.hash} \
          --doarchive {params.doarchive} --archive {params.archive} --cloud {params.cloud} \
          --output {output} \
          2>> {log}

          ## export rule env details
          conda env export --no-builds > info/fastqc.info
        '''


## Run RSEQC bam_stat.py
rule bam_qc:
    input:
        bam=rules.run_bwa.output,
        idx=rules.index_bam.output
    output:
        paths.rseqc.bamqc_txt
    benchmark:
        'benchmark/{sample}_bam_qc.tab'
    log:
        'log/{sample}_bam_qc.log'
    conda:
        SOURCEDIR+"/../envs/rseqc.yaml"
    params:
        sample='{sample}',
        rseqdir=RSEQC,
        srcdir=SOURCEDIR,
        doencrypt=DOENCRYPT,
        openssl=OPENSSL,
        encrypt_pass=ENCRYPT_PASS,
        hash=HASH,
        doarchive=DOARCHIVE,
        archive=ARCHIVE,
        cloud=CLOUD
    priority: 1
    threads: 1
    shell:
        '''
          echo "{params.rseqdir}bam_stat.py -i {input.bam} > {output}" | tee {log}
          {params.rseqdir}bam_stat.py -i {input.bam} > {output} 2>> {log}

          ## encrypt and archive if needed
          python3 {params.srcdir}/python/init-encrypt-archive.py --src {params.srcdir}\
          --doencrypt {params.doencrypt} --openssl {params.openssl} --encrypt-pass {params.encrypt_pass} --hash {params.hash} \
          --doarchive {params.doarchive} --archive {params.archive} --cloud {params.cloud} \
          --output {output} \
          2>> {log}
        '''


## Run RSEQC read_gc.py
rule bam_gc:
    input:
        bam=rules.run_bwa.output,
        idx=rules.index_bam.output
    output:
        r=paths.rseqc.bamgc_r,
        txt=paths.rseqc.bamgc_txt
    benchmark:
        'benchmark/{sample}_bam_gc.tab'
    log:
        'log/{sample}_bam_gc.log'
    conda:
        SOURCEDIR+"/../envs/rseqc.yaml"
    params:
        sample='{sample}',
        rseqdir=RSEQC,
        srcdir=SOURCEDIR,
        doencrypt=DOENCRYPT,
        openssl=OPENSSL,
        encrypt_pass=ENCRYPT_PASS,
        hash=HASH,
        doarchive=DOARCHIVE,
        archive=ARCHIVE,
        cloud=CLOUD
    priority: 1
    threads: 1
    shell:
      '''
        echo "{params.rseqdir}read_GC.py -i {input.bam} -o rseqc/{params.sample}" | tee {log}
        {params.rseqdir}read_GC.py -i {input.bam} -o rseqc/{params.sample} 2>> {log}

        ## R script to get txt output info
        echo "out=as.vector(summary(gc));dta = data.frame('{params.sample}',out[1],out[2],out[3],out[4],out[5],out[6]);write.table(dta,file='{output.txt}',sep="\t",row.names=F,col.names=F,quote=F);" >> {output.r}
        sed -i "s/pdf/png/g" {output.r} 
        sed -i 's/main=""/main="{params.sample}"/g' {output.r} 
        Rscript --vanilla --quiet {output.r}

        ## encrypt and archive if needed
        python3 {params.srcdir}/python/init-encrypt-archive.py --src {params.srcdir}\
        --doencrypt {params.doencrypt} --openssl {params.openssl} --encrypt-pass {params.encrypt_pass} --hash {params.hash} \
        --doarchive {params.doarchive} --archive {params.archive} --cloud {params.cloud} \
        --output {output.txt} \
        2>> {log}

        ## export rule env details
        conda env export --no-builds > info/rseqc.info
      '''

