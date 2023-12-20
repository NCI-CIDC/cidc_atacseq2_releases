## convert GTF formatted annotations to BED [note GTF format is 1 based and BED is 0 based]
rule annot_gtf2bed:
    input:
        gtf=paths.annot.gtf
    output:
        bed=paths.annot.bed
    priority: 1
    threads: 1
    shell:
       '''
         awk 'OFS="\t" {{if ($3=="gene") {{print $1,$4-1,$5,"1",$10,$7,$4-1,$5}}}}' {input} | tr -d '";' > {output}
       '''
## convert narrow peak to BED format
rule targets_narrowPeakToBed:
    input:
        narrow_peak_sorted=paths.peak.narrowPeak
    output:
        narrow_peak_bed=paths.peak.narrowPeak_bed
    priority: 1
    threads: 1
    shell:
        "cut -f1,2,3,4,9 {input} > {output}"

## report all gene scores with a decay rate of 1K
rule targets_getNarrowPeaksRPScore1k:
    input:
        narrow_peak_bed=paths.peak.narrowPeak_bed,
        anot_bed=paths.annot.bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_1k
    params:
        decay=1000,
        srcdir=SOURCEDIR
    priority: 1
    threads: 1
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -p {input.narrow_peak_bed} -a {input.anot_bed} -d {params.decay} -n {output}"

## report all gene scores with a decay rate of 10K
rule targets_getNarrowPeaksRPScore10k:
    input:
        narrow_peak_bed=paths.peak.narrowPeak_bed,
        anot_bed=paths.annot.bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_10k
    params:
        decay=10000,
        srcdir=SOURCEDIR
    priority: 1
    threads: 1
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -p {input.narrow_peak_bed} -a {input.anot_bed} -d {params.decay} -n {output}"

## report all gene scores with a decay rate of 100K
rule targets_getNarrowPeaksRPScore100k:
    input:
        narrow_peak_bed=paths.peak.narrowPeak_bed,
        anot_bed=paths.annot.bed
    output:
        gene_targets_1k=paths.targets.narrowPeak_100k
    params:
        decay=100000,
        srcdir=SOURCEDIR
    priority: 1
    threads: 1
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -p {input.narrow_peak_bed} -a {input.anot_bed} -d {params.decay} -n {output}"

## convert broad peak to BED format
rule targets_broadPeakToBed:
    input:
        broad_peak_sorted=paths.peak.broadPeak
    output:
        broad_peak_bed=paths.peak.broadPeak_bed
    priority: 1
    threads: 1
    shell:
        "cut -f1,2,3,4,9 {input} > {output}"

## report all gene scores with a decay rate of 1K
rule targets_getBroadPeaksRPScore1k:
    input:
        broad_peak_bed=paths.peak.broadPeak_bed,
        anot_bed=paths.annot.bed
    output:
        gene_targets_1k=paths.targets.broadPeak_1k
    params:
        decay=1000,
        srcdir=SOURCEDIR
    priority: 1
    threads: 1
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -p {input.broad_peak_bed} -a {input.anot_bed} -d {params.decay} -n {output}"

## report all gene scores with a decay rate of 10K
rule targets_getBroadPeaksRPScore10k:
    input:
        broad_peak_bed=paths.peak.broadPeak_bed,
        anot_bed=paths.annot.bed
    output:
        gene_targets_1k=paths.targets.broadPeak_10k
    params:
        decay=10000,
        srcdir=SOURCEDIR
    priority: 1
    threads: 1
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -p {input.broad_peak_bed} -a {input.anot_bed} -d {params.decay} -n {output}"

## report all gene scores with a decay rate of 100K
rule targets_getBroadPeaksRPScore100k:
    input:
        broad_peak_bed=paths.peak.broadPeak_bed,
        anot_bed=paths.annot.bed
    output:
        gene_targets_1k=paths.targets.broadPeak_100k
    params:
        decay=100000,
        srcdir=SOURCEDIR
    priority: 1
    threads: 1
    shell:
        "python3 {params.srcdir}/python/targets_RegPotential_Version2.py -p {input.broad_peak_bed} -a {input.anot_bed} -d {params.decay} -n {output}"
