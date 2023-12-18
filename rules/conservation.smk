#MODULE: conservation- module to create conservation plots

# _logfile=output_path + "/logs/conservation.log"
#_numPngs is used in conservation_plot rule to see how many pngs to expect
#note: the rule plots 3 runs per png, so for example, 12 runs results in 4 pngs
#_nPerPlot = 3
#_numPngs = math.ceil(len(config['runs'].keys())/float(_nPerPlot))
#_nPngs = [n+1 for n in range(_numPngs)]

#NOTE: using the _refs from chips.snakefile
def conservation_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            # ls.append(output_path + "/peaks/%s/%s_sorted_5k_summits.bed" % (runRep,runRep))
            if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
                ls.append(output_path + "/peaks/%s/%s_sorted_5k_peaks.bed" % (runRep,runRep))
            else:
                ls.append(output_path + "/peaks/%s/%s_sorted_5k_summits.bed" % (runRep,runRep))
            ls.append(output_path + "/conserv/%s/%s_conserv.R" % (runRep,runRep))
            ls.append(output_path + "/conserv/%s/%s_conserv.png" % (runRep,runRep))
            ls.append(output_path + "/conserv/%s/%s_conserv_thumb.png" % (runRep,runRep))
    return ls


#filtered_sorted_narrowPeak
#filtered_sorted_summits

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