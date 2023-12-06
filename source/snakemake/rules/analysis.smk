## Run CNV analysis with QDNAseq
rule cnv_analysis:
   input:
       bam=rules.tn5_adjust_bam.output.adj_bam,
       idx=rules.filter_bam.output.index
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
       srcdir=SOURCEDIR,
       predir=PREDIR,
       annot_gtf=paths.annot.gtf
   priority: 1
   threads: 1
   shell:
       '''
         echo "Rscript --vanilla {params.srcdir}/r/init-qdnaseq-cnv-analysis.r {params.predir} {params.sample} {params.annot_gtf}" | tee {log}
         Rscript --vanilla {params.srcdir}/r/init-qdnaseq-cnv-analysis.r {params.predir} {params.sample} {params.annot_gtf} 2>> {log}

         ## export rule env details
         conda env export --no-builds > info/cnv_analysis.info
       '''

## Run peak calling with MACS
rule call_peaks:
    input:
        bam=rules.tn5_adjust_bam.output.adj_bam,
        idx=rules.filter_bam.output.index
    output:
        xls=paths.peak.xls,
        peak=paths.peak.peak_narrow if PEAK_MODE=='narrow' else paths.peak.peak_broad,
        extra=paths.peak.extra_narrow if PEAK_MODE=='narrow' else paths.peak.extra_broad,
        pileup_bdg=paths.peak.pileup_bdg,
        lambda_bdg=paths.peak.lambda_bdg
    benchmark:
        'benchmark/{sample}_call_peaks.tab'
    log:
        'log/{sample}_call_peaks.log'
    conda:
        SOURCEDIR+"/../envs/call_peaks.yaml"
    params:
        sample='{sample}',
        format='BAMPE' if len(ENDS)==2 else 'BAM',
        genome='hs',
        fdr=PEAK_FDR,
        peak_mode_str=PEAK_MODE_STR
    priority: 1
    threads: 1
    shell:
        '''
          ## call peaks with macs
          echo "macs2 callpeak -f {params.format} -g {params.genome} -t {input.bam} --outdir peak -n {params.sample} -q {params.fdr} {params.peak_mode_str} --keep-dup all --bdg" | tee {log}
          macs2 callpeak -f {params.format} -g {params.genome} -t {input.bam} --outdir peak -n {params.sample} -q {params.fdr} {params.peak_mode_str} --keep-dup all --bdg 2>> {log}

          ## export rule env details
          conda env export --no-builds > info/call_peaks.info
        '''

## Convert peak bdg file to bw file
rule bdg_to_bw:
    input:
        genome_size=rules.genome_size.output.size,
        peak=rules.call_peaks.output.peak,
        extra=rules.call_peaks.output.extra,
        pileup_bdg=rules.call_peaks.output.pileup_bdg
    output:
        paths.peak.bw
    benchmark:
        'benchmark/{sample}_bdg_to_bw.tab'
    log:
        'log/{sample}_bdg_to_bw.log'
    conda:
        SOURCEDIR + "/../envs/bdg_to_bw.yaml"
    params:
        sample='{sample}'
    priority: 1
    threads: 1
    shell:
        '''
          ## convert bdg to bw
          echo "bedGraphToBigWig {input.pileup_bdg} {input.genome_size} {output}" | tee {log}
          bedGraphToBigWig {input.pileup_bdg} {input.genome_size} {output} 2>> {log}

          ## export rule env details
          conda env export --no-builds > info/bdg_to_bw.info
        '''

## Get peak metrics per sample with ChIPQC
rule chipqc:
    input:
        bam=rules.sample_bam.output.sampled_bam,
        idx=rules.sample_bam.output.index,
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
        srcdir=SOURCEDIR
    priority: 1
    threads: 1
    shell:
        '''
          echo "Rscript --vanilla --quiet {params.srcdir}/r/init-chipqc-sample-peak-metrics.r {params.predir} {params.sample} {input.bam} {input.peak}" | tee {log}
          Rscript --vanilla --quiet {params.srcdir}/r/init-chipqc-sample-peak-metrics.r {params.predir} {params.sample} {input.bam} {input.peak} 2>> {log}

          ## export rule env details
          conda env export --no-builds > info/chipqc.info
        '''

## Annotate peaks and perform functional pathway enrichment
rule annot_peaks:
    input:
        peak=rules.call_peaks.output.peak
    output:
        stat=paths.peak.annot_stat,
        tab=paths.peak.annot_tab,
        gobp=paths.ptw.gobp,
        kegg=paths.ptw.kegg
    benchmark:
        'benchmark/{sample}_annot_peaks.tab'
    log:
        'log/{sample}_annot_peaks.log'
    conda:
        SOURCEDIR+"/../envs/annot_peaks.yaml"
    params:
        sample='{sample}',
        annot_gtf=paths.annot.gtf,
        srcdir=SOURCEDIR,
        predir=PREDIR
    priority: 1
    threads: 1
    shell:
        '''
          echo "Rscript --vanilla --quiet {params.srcdir}/r/init-peak-annot-ptw-enr.r {params.predir} {params.sample} {params.annot_gtf} {input.peak}" | tee {log}
          Rscript --vanilla --quiet {params.srcdir}/r/init-peak-annot-ptw-enr.r {params.predir} {params.sample} {params.annot_gtf} {input.peak} 2>> {log}

          ## export rule env details
          conda env export --no-builds > info/annot_peaks.info
        '''
