rule conservation_plotConservation:
    """generate conservation plots"""
    input:
        peak=paths.peak.filtered_sorted_narrowPeak if PEAK_MODE=='narrow' else paths.peak.filtered_sorted_broadPeak,
	bw = paths.annot.bw
    output:
        png = paths.conservation.png,
        thumb = paths.conservation.thumb,
        r = paths.conservation.r,
        score = paths.conservation.score
    params:
        srcdir = SOURCEDIR,
        width=config["width"],
#	sample = "{sample}"
        main_output_path=Path(paths.conservation.score).parent
    message: "CONSERVATION: calling conservation script"
    log: to_log(paths.conservation.score)
    benchmark: to_benchmark(paths.conservation.score)
    conda: "../envs/conservation.yaml"
    shell:
         "pwd;"
         "python3 {params.srcdir}/python/conservation_onebw_plot.py -t Conservation_at_summits -d {input.bw} "
         "-o  {params.main_output_path}/{wildcards.sample} "
	 "-l Peak_summits {input.peak} -w {params.width} > {output.score} 2>>{log}"