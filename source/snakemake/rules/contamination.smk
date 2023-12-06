rule contamination_centrifuge_index:
    input:
        "p+h+v.tar.gz"
    output
        "centrifuge_index"
    log:
       "logs/microbiome/{sample}.centrifuge.log"
    params:
       path="set +eu;source activate %s" % config['centrifuge_root'],
    message:
       "Building Centrifuge Index"
    benchmark:
       "benchmarks/microbiome/centrifuge.index.benchmark"
    threads: _microbiome_threads     
    conda: "../envs/contamination.yml"
    shell:
        '''centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map;'''
        '''centrifuge-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 9606 -c 'reference genome' >> seqid2taxid.map;'''
        '''cat library/*/*.fna > input-sequences.fna;'''
        '''centrifuge-build -p 4 --conversion-table seqid2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 input-sequences.fna abv'''


rule contamination_centrifuge:
     input:
         fa=rules.getfile.output,
         index=rules.contamination_centrifuge_index.output
     output:
        classfication = "analysis/microbiome/{sample}/{sample}_classification.txt.gz",
        report = "analysis/microbiome/{sample}/{sample}_report.txt",
     log:
        "logs/microbiome/{sample}.centrifuge.log"
     params:
        centrifuge_index = config["centrifuge_index"],
        path="set +eu;source activate %s" % config['centrifuge_root'],
     message:
        "Running Centrifuge on {sample}"
     benchmark:
        "benchmarks/microbiome/{sample}.microbiome.benchmark"
     threads: _microbiome_threads     
     conda: "../envs/contamination.yml"
     shell:
        """{params.path}; centrifuge -x {params.centrifuge_index} -p {threads}  --host-taxids 9606 -1 {input[0]} -2 {input[1]} -S {output.classfication} --report-file {output.report} """
        """ && gzip {output.classfication} """
        """ && awk '{{print FILENAME}}' {output.report} | awk -F '/' '{{print $3}}' | paste - {output.report} | awk -F '\t' 'NR==1{{$1="sample"}}1' OFS='\t' > {output.add_sample} """
