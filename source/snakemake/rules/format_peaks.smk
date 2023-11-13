# Sort summits.bed from rule call_peaks by the fifth column (qvalue)
rule sort_summits:
    input:
        summits=rules.call_peaks.output.extra
    output:
        sorted_summits=paths.format_peaks.summits
    shell:
        "sort -r -n -k 5 {input} > {output}"

# Format file for HOMER analysis by selecting the top 5k peaks and prepend column 1 values with "chr"
rule filter_5k_sorted_summits:
    input:
        sorted_summits=paths.format_peaks.summits
    output:
        filtered_sorted_summits=paths.format_peaks.filtered_sorted_summits,
        temp_filtered_sorted_summits=temp(paths.format_peaks.filtered_sorted_summits+".temp")
    shell:
        """
        head -n 5000 {input} > {output.temp_filtered_sorted_summits}
        awk -F'\t\' 'BEGIN {{OFS = "\t"}} {{ $1 = "chr" $1; print }}' {output.temp_filtered_sorted_summits} > {output.filtered_sorted_summits}
        """

# Sort peak.narrowPeak from rule call_peaks by the ninth column (qvalue)
rule sort_narrowPeak:
    input:
        narrowPeak=paths.peak.peak_narrow
    output:
        sorted_narrowPeak=paths.format_peaks.narrowPeak
    shell:
       "sort -r -n -k 9 {input} > {output}"

# Format file for HOMER analysis by selecting the top 5k peaks and prepend column 1 values with "chr"
rule filter_5k_sorted_narrowPeak:
    input:
        sorted_narrowPeak=paths.format_peaks.narrowPeak
    output:
        filtered_sorted_narrowPeak=paths.format_peaks.filtered_sorted_narrowPeak,
        temp_filtered_sorted_narrowPeak=temp(paths.format_peaks.filtered_sorted_narrowPeak+".temp")
    shell:
        """
        head -n 5000 {input} > {output.temp_filtered_sorted_narrowPeak}
        awk -F'\t\' 'BEGIN {{OFS = "\t"}} {{ $1 = "chr" $1; print }}' {output.temp_filtered_sorted_narrowPeak} > {output.filtered_sorted_narrowPeak}
        """

# Sort peak.broadPeak from rule call_peaks by the ninth column (qvalue)
rule sort_broadPeak:
    input:
        broadPeak=paths.peak.peak_broad
    output:
        sorted_broadPeak=paths.format_peaks.broadPeak
    shell:
       "sort -r -n -k 9 {input} > {output}"

# Format file for HOMER analysis by selecting the top 5k peaks and prepend column 1 values with "chr"
rule filter_5k_sorted_broadPeak:
    input:
        sorted_broadPeak=paths.format_peaks.broadPeak
    output:
        filtered_sorted_broadPeak=paths.format_peaks.filtered_sorted_broadPeak,
        temp_filtered_sorted_broadPeak=temp(paths.format_peaks.filtered_sorted_broadPeak+".temp")
    shell:
        """
        head -n 5000 {input} > {output.temp_filtered_sorted_broadPeak}
        awk -F'\t\' 'BEGIN {{OFS = "\t"}} {{ $1 = "chr" $1; print }}' {output.temp_filtered_sorted_broadPeak} > {output.filtered_sorted_broadPeak}
        """
