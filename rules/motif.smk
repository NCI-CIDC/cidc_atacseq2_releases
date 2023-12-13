# Install and configure hg38 for use with HOMER
rule configure_hg38:
    output: 
        hg38_path=PREDIR+"/motif/homer_hg38"
    conda:
       SOURCEDIR+"/../envs/homer.yaml"
    priority: 1
    threads: 1
    shell:
       """
       # Identifies path for where HOMER and the script configureHomer.pl were installed
       homer_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}')
       homer_path="${{homer_path}}/share/homer"
       
       # Installs and configures hg38 for use with HOMER
       perl ${{homer_path}}/configureHomer.pl -install hg38

       # Generates a symbolic link to the hg38 genome
       ln -s ${{homer_path}}/data/genomes/hg38/genome.fa {output}
       """

# Identify motifs with HOMER for the top 5k peaks from peaks.narrowPeak
rule find_motifs_narrowPeak:
    input:
        hg38_path=rules.configure_hg38.output,
        filtered_sorted_narrowPeak=paths.peak.filtered_sorted_narrowPeak
    output:
        directory(paths.motif.narrow_peak)
    benchmark:
        "benchmark/{sample}_find_motifs_narrowPeak.tab"
    log:
        "log/{sample}_find_motifs_narrowPeak.log"
    conda:
        SOURCEDIR+"/../envs/homer.yaml"
    params:
        preparsed_dir=PREDIR+"/motif/preparsed/"
    priority: 1
    threads: 1
    shell:
        """
        # Identify path for where HOMER was installed in order to use the findMotifsGenome.pl
        echo "homer_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}')" | tee {log}
        homer_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}') 2>> {log}

        echo "${{homer_path}}/share/homer/bin/findMotifsGenome.pl {input.filtered_sorted_narrowPeak} {input.hg38_path} {output} -size 200 -mask -preparsedDir {params.preparsed_dir}" | tee {log}
        ${{homer_path}}/share/homer/bin/findMotifsGenome.pl {input.filtered_sorted_narrowPeak} {input.hg38_path} {output} -size 200 -mask -preparsedDir {params.preparsed_dir} 2>> {log}

        # Export rule env details
        conda env export --no-builds > info/find_motifs_narrowPeak.info
        """

# Identify motifs with HOMER for the top 5k peaks from summits.bed
rule find_motifs_summits:
    input:
        hg38_path=rules.configure_hg38.output,
        filtered_sorted_summits=paths.peak.filtered_sorted_summits
    output:
        directory(paths.motif.summit)
    benchmark:
        "benchmark/{sample}_find_motifs_summits.tab"
    log:
        "log/{sample}_find_motifs_summits.log"
    conda:
        SOURCEDIR+"/../envs/homer.yaml"
    params:
        preparsed_dir=PREDIR+"/motif/preparsed/"
    priority: 1
    threads: 1
    shell:
        """
        # Identify path for where HOMER was installed in order to use the findMotifsGenome.pl
        echo "homer_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}')" | tee {log}
        homer_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}') 2>> {log}

        echo "${{homer_path}}/share/homer/bin/findMotifsGenome.pl {input.filtered_sorted_summits} {input.hg38_path} {output} -size 200 -mask -preparsedDir {params.preparsed_dir}" | tee {log}
        ${{homer_path}}/share/homer/bin/findMotifsGenome.pl {input.filtered_sorted_summits} {input.hg38_path} {output} -size 200 -mask -preparsedDir {params.preparsed_dir} 2>> {log}

        # Export rule env details
        conda env export --no-builds > info/find_motifs_summits.info
        """

# Identify motifs with HOMER for the top 5k peaks from peaks.broadPeak
rule find_motifs_broadPeak:
    input:
        hg38_path=rules.configure_hg38.output,
        filtered_sorted_broadPeak=paths.peak.filtered_sorted_broadPeak
    output:
        directory(paths.motif.broad_peak)
    benchmark:
        "benchmark/{sample}_find_motifs_broadPeak.tab"
    log:
        "log/{sample}_find_motifs_broadPeak.log"
    conda:
        SOURCEDIR+"/../envs/homer.yaml"
    params:
        preparsed_dir=PREDIR+"/motif/preparsed/"
    priority: 1
    threads: 1
    shell:
        """
        # Identify path for where HOMER was installed in order to use the findMotifsGenome.pl
        echo "homer_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}')" | tee {log}
        homer_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}') 2>> {log}

        echo "${{homer_path}}/share/homer/bin/findMotifsGenome.pl {input.filtered_sorted_broadPeak} {input.hg38_path} {output} -size 1000 -mask -preparsedDir {params.preparsed_dir}" | tee {log}
        ${{homer_path}}/share/homer/bin/findMotifsGenome.pl {input.filtered_sorted_broadPeak} {input.hg38_path} {output} -size 1000 -mask -preparsedDir {params.preparsed_dir} 2>> {log}

        # Export rule env details
        conda env export --no-builds > info/find_motifs_broadPeak.info
        """
