## Generate genome tracks
rule genome_track:
    input:
        peak=rules.call_peaks.output.peak,
        bw=rules.bdg_to_bw.output
    output:
        ini=paths.track.ini,
        png=paths.track.png
    benchmark:
        'benchmark/{sample}_genome_track.tab'
    log:
        'log/{sample}_genome_track.log'
    conda:
        SOURCEDIR + "/../envs/genome_track.yaml"
    params:
        sample='{sample}',
        track_region=TRACK_REGION,
        annot_gtf=paths.annot.gtf,
        output_joined=','.join([paths.track.ini, paths.track.png]),
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
          ## make track file
          echo "make_tracks_file --trackFiles {input.bw} {params.annot_gtf} -o {output.ini}" | tee {log}
          make_tracks_file --trackFiles {input.bw} {params.annot_gtf} -o {output.ini} 2>> {log}
          ## turn on gene labeling and merge transcripts into one structure per gene
          sed -i "s/^labels = false/labels = true/g" {output.ini}
          sed -i "s/^# merge_transcripts = true/merge_transcripts = true/g" {output.ini}
          ## change font size
          sed -i "s/^fontsize = 10/fontsize = 8/g" {output.ini}
          ## change annot label
          sed -i "s/^title = GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation/title = GRCh38/g" {output.ini}

          ## create track image
          echo "pyGenomeTracks --tracks {output.ini} --region {params.track_region} -o {output.png} --dpi 300" | tee -a {log}
          pyGenomeTracks --tracks {output.ini} --region {params.track_region} -o {output.png} --dpi 300 2>> {log}

          ## encrypt and archive if needed
          python3 {params.srcdir}/python/init-encrypt-archive.py --src {params.srcdir}\
          --doencrypt {params.doencrypt} --openssl {params.openssl} --encrypt-pass {params.encrypt_pass} --hash {params.hash} \
          --doarchive {params.doarchive} --archive {params.archive} --cloud {params.cloud} \
          --output {params.output_joined} \
          2>> {log}

          ## export rule env details
          conda env export --no-builds > info/genome_track.info
        '''