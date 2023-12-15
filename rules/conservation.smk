
# _logfile=output_path + "/logs/conservation.log"
#_numPngs is used in conservation_plot rule to see how many pngs to expect
#note: the rule plots 3 runs per png, so for example, 12 runs results in 4 pngs
_nPerPlot = 3
_numPngs = math.ceil(len(config['runs'].keys())/float(_nPerPlot))
_nPngs = [n+1 for n in range(_numPngs)]

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

def conservationInput(wildcards):
    run = wildcards.run
    rep = wildcards.rep
    runRep = "%s.%s" % (run,rep)
    if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
        temp = output_path + "/peaks/%s/%s_sorted_peaks.bed" % (runRep,runRep)
    else:
        temp = output_path + "/peaks/%s/%s_sorted_summits.bed" % (runRep,runRep)
    return temp

#filtered_sorted_narrowPeak
#filtered_sorted_summits


rule conservation_plotConservation:
    """generate conservation plots"""
    input:
        if macs2_broadpeaks:
           filtered_sorted_narrowPeak
        else:
           filtered_sorted_summits
    output:
        png = paths.conservation.png,
        thumb = paths.conservation.thumb,
        r = paths.conservation.r,
        score = paths.conservation.score
    params:
#        db=config['conservation'],
        db="/media/storage/hg38.phastCons100way.bw",
        width=4000,
        #run = lambda wildcards: wildcards.run,
#        run="{run}.{rep}" ,
        main_output_path=output_path
    message: "CONSERVATION: calling conservation script"
    log: output_path + "/logs/conservation/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_conservation_plotConservation.benchmark"
    conda: "../envs/conservation.yaml"
    shell:
         "source/python/conservation_onebw_plot.py -t Conservation_at_summits -d {params.db} "
         "-o  {params.main_output_path}/conserv/{params.run}/{params.run}_conserv "
         "-l Peak_summits {input} -w {params.width} > {output.score} 2>>{log}"	