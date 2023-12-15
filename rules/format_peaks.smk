# Sort summits.bed from rule call_peaks by the fifth column (qvalue) and filter for top 5k summits
rule filter_5k_sorted_summits:
    input:
        summits=paths.peak.extra_narrow
    output:
        sorted_summits=paths.peak.summits,
        filtered_summits=paths.peak.filtered_sorted_summits
    priority: 1
    threads: 1
    shell:
        '''
          sort -r -n -k 5 {input.summits} > {output.sorted_summits}
          head -n 5000 {output.sorted_summits} > {output.filtered_summits}
        '''

# Sort peak.narrowPeak from rule call_peaks by the ninth column (qvalue) and filter for top 5k peaks
rule filter_5k_sorted_narrowPeak:
    input:
        narrowPeak=paths.peak.peak_narrow
    output:
        sorted_narrowPeak=paths.peak.narrowPeak,
        filtered_narrowPeak=paths.peak.filtered_sorted_narrowPeak
    priority: 1
    threads: 1
    shell:
        '''
          sort -r -n -k 9 {input.narrowPeak} > {output.sorted_narrowPeak}
          head -n 5000 {output.sorted_narrowPeak} > {output.filtered_narrowPeak}
        '''


# Sort peak.broadPeak from rule call_peaks by the ninth column (qvalue) and filter for top 5k peaks
rule filter_5k_sorted_broadPeak:
    input:
        broadPeak=paths.peak.peak_broad
    output:
        sorted_broadPeak=paths.peak.broadPeak,
        filtered_broadPeak=paths.peak.filtered_sorted_broadPeak
    priority: 1
    threads: 1
    shell:
       '''
         sort -r -n -k 9 {input.broadPeak} > {output.sorted_broadPeak}
         head -n 5000 {output.sorted_broadPeak} > {output.filtered_broadPeak}
       '''



