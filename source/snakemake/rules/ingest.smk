## Get input
rule getfile:
    input:
        rules.directory_setup.output,
        rules.build_bwa_index.output
    output:
        temp(expand(paths.input.input_fastq, read=ENDS))
    benchmark:
        'benchmark/{sample}_getfile.tab'
    log:
        'log/{sample}_getfile.log'
    conda:
        SOURCEDIR+"/../envs/getfile.yaml"
    params:
        sample='{sample}',
        srcdir=SOURCEDIR,
        ends=','.join(ENDS),
        samid=','.join(SAMID),
        fastq1=','.join(FASTQ_1),
        fastq2=','.join(FASTQ_2),
        cloud=CLOUD,
        fastqdump=FASTQDUMP,
        openssl=OPENSSL,
        decrypt_pass=DECRYPT_PASS,
        hash=HASH,
        output_joined=','.join(expand('input/{{sample}}_{pe}.fastq.gz',  pe=ENDS))
    priority: 2
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo python3 {params.srcdir}/python/init-getfile.py --src {params.srcdir} --ends {params.ends} \
          --samid {params.samid} --fastq1 {params.fastq1} --fastq2 {params.fastq2} \
          --cloud {params.cloud} --fastqdump {params.fastqdump} --openssl {params.openssl} --decrypt-pass {params.decrypt_pass} --hash {params.hash} \
          --sample {params.sample} --output {params.output_joined} \
          > {log}

          python3 {params.srcdir}/python/init-getfile.py --src {params.srcdir} --ends {params.ends} \
          --samid {params.samid} --fastq1 {params.fastq1} --fastq2 {params.fastq2} \
          --cloud {params.cloud} --fastqdump {params.fastqdump} --openssl {params.openssl} --decrypt-pass {params.decrypt_pass} --hash {params.hash} \
          --sample {params.sample} --output {params.output_joined} \
          2>> {log}

          ## export rule env details
          conda env export --no-builds > info/getfile.info
        '''

### Optionally trim adapters with cutadapt
rule trimadapters:
    input:
        fa=rules.getfile.output
    output:
        temp([x + TRIM_ADAPTERS_OUTPUT for x in expand(paths.cutadapt.cutadapt_fastq, read=ENDS)])
    benchmark:
        'benchmark/{sample}_trimadapters.tab'
    log:
        'log/{sample}_trimadapters.log'
    conda:
         SOURCEDIR+"/../envs/trimadpaters.yaml"
    params:
        sample='{sample}',
        srcdir=SOURCEDIR,
        samid=','.join(SAMID),
        trim_fp=TRIM_FP,
        trim_tp=TRIM_TP,
        fp_adapters=','.join(FP_ADAPTERS),
        tp_adapters=','.join(TP_ADAPTERS),
        cutadapt=CUTADAPT,
        input_joined=','.join(expand('input/{{sample}}_{pe}.fastq.gz',  pe=ENDS)),
        output_joined=','.join([x + TRIM_ADAPTERS_OUTPUT for x in expand('cutadapt/{{sample}}_{pe}', pe=ENDS)])
    priority: 3
    threads: max(1,min(8,NCORES))
    shell:
      '''
        echo python3 {params.srcdir}/python/init-trimadapters.py --src {params.srcdir} --trim-fp {params.trim_fp} --trim-tp {params.trim_tp}\
        --samid {params.samid} --fp-adapters {params.fp_adapters} --tp-adapters {params.tp_adapters} \
        --cutadapt {params.cutadapt} --sample {params.sample} --input {params.input_joined} --output {params.output_joined} --threads {threads}\
        > {log}

        python3 {params.srcdir}/python/init-trimadapters.py --src {params.srcdir} --trim-fp {params.trim_fp} --trim-tp {params.trim_tp}\
        --samid {params.samid} --fp-adapters {params.fp_adapters} --tp-adapters {params.tp_adapters} \
        --cutadapt {params.cutadapt} --sample {params.sample} --input {params.input_joined} --output {params.output_joined} --threads {threads} --log {log}\
        2>> {log}

        ## export rule env details
        conda env export --no-builds > info/trimadapters.info
      '''


## Optionally quality-trim reads with Trimmomatic
rule qualityfilter:
    input:
        rules.trimadapters.output if TRIM_FP or TRIM_TP else rules.getfile.output
    output:
        temp(expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])) if len(ENDS)==2 else temp(expand(paths.rqual_filter.qfilter_fastq_single), read=ENDS)
    benchmark:
        'benchmark/{sample}_qualityfilter.tab'
    log:
        'log/{sample}_qualityfilter.log'
    conda:
        SOURCEDIR+"/../envs/qualityfilter.yaml"
    params:
        sample='{sample}',
        srcdir=SOURCEDIR,
        ends=','.join(ENDS),
        trimmomatic=TRIMMOMATIC,
        qual_cutoff=QUAL_CUTOFF,
        input_joined=','.join([x + TRIM_ADAPTERS_OUTPUT for x in expand('cutadapt/{{sample}}_{pe}', pe=ENDS)] if TRIM_FP or TRIM_TP else expand('input/{{sample}}_{pe}.fastq.gz',  pe=ENDS)),
        output_joined=','.join(expand('rqual_filter/{{sample}}_{pe}{paired}_qual.fastq.gz', pe=ENDS, paired=['P','U']) if len(ENDS)==2 else expand('rqual_filter/{{sample}}_{pe}_qual.fastq.gz', pe=ENDS))
    priority: 3
    threads: max(1,min(8,NCORES))
    shell:
      '''
        echo python3 {params.srcdir}/python/init-qualityfilter.py --src {params.srcdir} --ends {params.ends} --threads {threads} \
        --trimmomatic {params.trimmomatic} --qual-cutoff {params.qual_cutoff} --input {params.input_joined} --output {params.output_joined} --log {log}\
        > {log}

        python3 {params.srcdir}/python/init-qualityfilter.py --src {params.srcdir} --ends {params.ends} --threads {threads} \
        --trimmomatic {params.trimmomatic} --qual-cutoff {params.qual_cutoff} --input {params.input_joined} --output {params.output_joined} --log {log} \
        2>> {log}

        ## export rule env details
        conda env export --no-builds > info/qualityfilter.info
      '''
