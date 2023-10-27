## Run CNV analysis with QDNAseq
rule cnv_analysis:
   input:
       bam=rules.run_bwa.output,
       idx=rules.index_bam.output
   output:
       bed=paths.cnv.bed,
       igv=paths.cnv.igv,
       csv=paths.cnv.csv
   benchmark:
       'benchmark/{sample}_cnv_analysis.tab'
   log:
       'log/{sample}_cnv_analysis.log'
   conda:
       SOURCEDIR+"/../envs/cnv_analysis.yaml"
   params:
       sample='{sample}',
       predir=PREDIR,
       annot_gtf=paths.annot.gtf,
       output_joined=','.join([paths.cnv.bed, paths.cnv.igv, paths.cnv.csv]),
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
         echo "Rscript --vanilla {params.srcdir}/r/init-qdnaseq-cnv-analysis.r {params.predir} {params.sample} {params.annot_gtf}" | tee {log}
         Rscript --vanilla {params.srcdir}/r/init-qdnaseq-cnv-analysis.r {params.predir} {params.sample} {params.annot_gtf} 2>> {log}

         ## encrypt and archive if needed
         python3 {params.srcdir}/python/init-encrypt-archive.py --src {params.srcdir}\
         --doencrypt {params.doencrypt} --openssl {params.openssl} --encrypt-pass {params.encrypt_pass} --hash {params.hash} \
         --doarchive {params.doarchive} --archive {params.archive} --cloud {params.cloud} \
         --output {params.output_joined} \
         2>> {log}

         ## export rule env details
         conda env export --no-builds > info/cnv_analysis.info
       '''

## TODO finalize - params based on SE/PE narrow/broad etc
## Run peak calling with MACS
rule call_peaks:
    input:
        bam=rules.run_bwa.output,
        idx=rules.index_bam.output
    output:
        xls=paths.peak.xls,
        peak=paths.peak.peak_narrow if PEAK_MODE=='narrow' else paths.peak.peak_broad,
        extra=paths.peak.extra_narrow if PEAK_MODE=='narrow' else paths.peak.extra_broad
    benchmark:
        'benchmark/{sample}_call_peaks.tab'
    log:
        'log/{sample}_call_peaks.log'
    conda:
        SOURCEDIR+"/../envs/call_peaks.yaml"
    params:
        sample='{sample}',
        macs=MACS,
        format='BAMPE' if len(ENDS)==2 else 'BAM',
        genome='hs',
        fdr=PEAK_FDR,
        peak_mode_str=PEAK_MODE_STR,
        output_joined=','.join([paths.peak.xls, paths.peak.peak_narrow if PEAK_MODE=='narrow' else paths.peak.peak_broad, paths.peak.extra_narrow if PEAK_MODE=='narrow' else paths.peak.extra_broad]),
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
          echo "{params.macs} callpeak -f {params.format} -g {params.genome} -t {input.bam} --outdir peak -n {params.sample} -q {params.fdr} {params.peak_mode_str}" | tee {log}
          {params.macs} callpeak -f {params.format} -g {params.genome} -t {input.bam} --outdir peak -n {params.sample} -q {params.fdr} {params.peak_mode_str} 2>> {log}

          ## encrypt and archive if needed
          python3 {params.srcdir}/python/init-encrypt-archive.py --src {params.srcdir}\
          --doencrypt {params.doencrypt} --openssl {params.openssl} --encrypt-pass {params.encrypt_pass} --hash {params.hash} \
          --doarchive {params.doarchive} --archive {params.archive} --cloud {params.cloud} \
          --output {params.output_joined} \
          2>> {log}

          ## export rule env details
          conda env export --no-builds > info/call_peaks.info
        '''

## Get peak metrics per sample with ChIPQC
rule chipqc:
    input:
        bam=rules.run_bwa.output,
        idx=rules.index_bam.output,
        peak=rules.call_peaks.output.peak,
    output:
        paths.chipqc.csv
    benchmark:
        'benchmark/{sample}_chipqc.tab'
    log:
        'log/{sample}_chipqc.log'
    conda:
        SOURCEDIR+"/../envs/chipqc.yaml"
    params:
        sample='{sample}',
        predir=PREDIR,
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
          echo "Rscript --vanilla --quiet {params.srcdir}/r/init-chipqc-sample-peak-metrics.r {params.predir} {params.sample} {input.bam} {input.peak}" | tee {log}
          Rscript --vanilla --quiet {params.srcdir}/r/init-chipqc-sample-peak-metrics.r {params.predir} {params.sample} {input.bam} {input.peak} 2>> {log}

          ## encrypt and archive if needed
          python3 {params.srcdir}/python/init-encrypt-archive.py --src {params.srcdir}\
          --doencrypt {params.doencrypt} --openssl {params.openssl} --encrypt-pass {params.encrypt_pass} --hash {params.hash} \
          --doarchive {params.doarchive} --archive {params.archive} --cloud {params.cloud} \
          --output {output} \
          2>> {log}

          ## export rule env details
          conda env export --no-builds > info/chipqc.info
        '''
    
