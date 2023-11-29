## convert GTF formatted annotations to BED [note GTF format is 1 based and BED is 0 based]
rule annot_gtf2bed:
        input:
            gtf=paths.annot.gtf
        output:
            bed=paths.annot.bed
        shell:
            "awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,"1",$10,$7,$4-1,$5}}' {input} | tr -d '";' > {output}"

## convert narrow peak to BED format
rule targets_narrowPeakToBed:
        input:
            narrow_peak_sorted=paths.format_peaks.narrowPeak
        output:
            narrow_peak_bed=paths.targets.narrowPeak_bed
        shell:
            "cut -f1,2,3,4,9 {input} > {output}"

## report all gene scores with a decay rate of 1K
rule targets_getNarrowPeaksRPScore1k:
    input:
        narrow_peak_bed=paths.targets.narrowPeak_bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_1k
    params:
        anot_bed=paths.annot.bed,
        decay=1000,
        srcdir=SOURCEDIR
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -a {params.anot_bed} -d {params.decay} -p {input} -n {output}"

## report all gene scores with a decay rate of 10K
rule targets_getNarrowPeaksRPScore10k:
    input:
        narrow_peak_bed=paths.targets.narrowPeak_bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_10k
    params:
        anot_bed=paths.annot.bed,
        decay=10000,
        srcdir=SOURCEDIR
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -a {params.anot_bed} -d {params.decay} -p {input} -n {output}"

## report all gene scores with a decay rate of 100K
rule targets_getNarrowPeaksRPScore100k:
    input:
        narrow_peak_bed=paths.targets.narrowPeak_bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_100k
    params:
        anot_bed=paths.annot.bed,
        decay=100000,
        srcdir=SOURCEDIR
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -a {params.anot_bed} -d {params.decay} -p {input} -n {output}"

## convert broad peak to BED format
rule targets_broadPeakToBed:
        input:
            broad_peak_sorted=paths.format_peaks.broadPeak
        output:
            broad_peak_bed=paths.targets.broadPeak_bed
        shell:
            "cut -f1,2,3,4,9 {input} > {output}"
     
## report all gene scores with a decay rate of 1K
rule targets_getBroadPeaksRPScore1k:
    input:
        broad_peak_bed=paths.targets.broadPeak_bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_1k
    params:
        anot_bed=paths.annot.bed,
        decay=1000,
        srcdir=SOURCEDIR
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -a {params.anot_bed} -d {params.decay} -p {input} -n {output}"

## report all gene scores with a decay rate of 10K
rule targets_getBroadPeaksRPScore10k:
    input:
        broad_peak_bed=paths.targets.broadPeak_bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_10k
    params:
        anot_bed=paths.annot.bed,
        decay=10000,
        srcdir=SOURCEDIR
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -a {params.anot_bed} -d {params.decay} -p {input} -n {output}"

## report all gene scores with a decay rate of 100K
rule targets_getBroadPeaksRPScore100k:
    input:
        broad_peak_bed=paths.targets.broadPeak_bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_100k
    params:
        anot_bed=paths.annot.bed,
        decay=100000,
        srcdir=SOURCEDIR
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -a {params.anot_bed} -d {params.decay} -p {input} -n {output}"