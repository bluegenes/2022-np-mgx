import os, sys
import pandas as pd
import glob

configfile: "inputs/pacbio_reads.conf"

out_dir = config.get('output_dir', 'output.gather_pacbio_reads')
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")

basename = config.get("basename", 'np-pacbio-reads')


sample_info = pd.read_csv(config['sample_info'])
SAMPLES = sample_info["name"].to_list()
sample_info.set_index('name', inplace=True)
search_databases = config['search_databases'] # must be dictionary


# check params are in the right format, build alpha-ksize combos and param strings
alphabet_info = config['alphabet_info']
alpha_ksize_scaled=[]
all_param_str=[]
for alpha, info in alphabet_info.items():
    scaled = info["scaled"]
    ksize = info["ksize"]
    if not isinstance(scaled, list):
        scaled = [scaled]
        config["alphabet_info"][alpha]["scaled"] = scaled
    if not isinstance(ksize, list):
        ksize=[ksize]
        config["alphabet_info"][alpha]["ksize"] = ksize
    # build a parameter for the right combinations
    #alpha_ksize += expand(f"{alpha}-k{{ksize}}", ksize = ksize)
    alpha_ksize_scaled += expand(f"{alpha}-k{{ksize}}-sc{{scaled}}", ksize = ksize, scaled=scaled)

onstart:
    print("------------------------------")
    print("taxonomic classification with sourmash gather")
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas, something went wrong :C\n")


rule all:
    input:
        ancient(os.path.join(out_dir, "sigs", f"{basename}.queries.zip")),
        #expand(os.path.join(out_dir, 'gather', f"{basename}.{{aks}}.gather-pathlist.txt"), aks=alpha_ksize_scaled)


# k-mer abundance trimming
rule kmer_trim_reads:
    input: 
        #reads = ancient(outdir + '/trim/{sample}.trim.fq.gz'),
        reads = lambda w: sample_info.at[w.sample, "reads"]
    output:
        protected(os.path.join(out_dir, "abundtrim", "{sample}.abundtrim.fq.gz"))
    conda: 'conf/env/trim.yml'
    resources:
        mem_mb = int(20e9/ 1e6),
    params:
        mem = 20e9,
    shell: """
            trim-low-abund.py -C 3 -Z 18 -M {params.mem} -V \
            {input.reads} -o {output} --gzip
    """


rule sourmash_sketch:
    input: os.path.join(out_dir, "abundtrim", "{sample}.abundtrim.fq.gz")
   # input: lambda w: sample_info.at[w.sample, "reads"] # non-abundtrimmed sample
    output:
        os.path.join(out_dir, "sigs", "{sample}.sig.gz")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sketch", "{sample}.sketch.log")
    benchmark: os.path.join(benchmarks_dir, "sketch", "{sample}.sketch.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch translate {input} -p k=7,k=10,protein,abund,scaled=200 \
                                          -p k=16,k=19,dayhoff,abund,scaled=200 \
                                          --name {wildcards.sample} -o {output} 2> {log}
        """

localrules: sig_cat
rule sig_cat:
    input: expand(os.path.join(out_dir, "sigs", "{sample}.sig.gz"), sample=SAMPLES)
    output:
        zipF=os.path.join(out_dir, "sigs", f"{basename}.queries.zip"),
    log: os.path.join(logs_dir, "sig_cat", f"{basename}.sigcat.log")
    benchmark: os.path.join(benchmarks_dir, "sig_cat", f"{basename}.sigcat.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sig cat {input} -o {output} 2> {log}
        """

rule gather_sig_from_zipfile:
    input:
        query_zip=ancient(os.path.join(out_dir, "sigs", f"{basename}.queries.zip")),
        databases = lambda w: search_databases[f"{w.alphabet}-k{w.ksize}"]
    output:
        prefetch_csv=os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.prefetch.csv'),
        gather_csv=os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.txt'),
    params:
        #threshold_bp = config['sourmash_database_threshold_bp'],
        alpha_cmd = lambda w: "--" + w.alphabet 
    resources:
        mem_mb=lambda wildcards, attempt: attempt *30000,
        time=400, #240,
        partition="med2",# med2
        #partition="bml", # "low2"
    log: os.path.join(logs_dir, "gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB is {input.database}"
        echo "DB is {input.database}" > {log}

        sourmash sig grep {wildcards.sample} {input.query_zip} {alpha_cmd} \
                 --ksize {wildcards.ksize} | sourmash gather - {input.databases} \
                 -o {output.gather_csv} -k {wildcards.ksize} --scaled {wildcards.scaled} {alpha_cmd} \
                 --save-prefetch-csv {output.prefetch_csv} > {output.gather_txt} 2> {log}
        touch {output.prefetch_csv}
        touch {output.gather_txt}
        touch {output.gather_csv}
        """
                #--threshold-bp={params.threshold_bp} \ 
                 #--picklist {input.query_picklist}:name:identprefix:exclude \


rule tax_annotate:
    input:
        gather = os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.with-lineages.csv'),
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash tax annotate -g {input.gather} -t {' -t'.join(input.lineages)}
        """




localrules: annotated_gather_csvs_to_pathlist
rule annotated_gather_csvs_to_pathlist:
    input: 
        expand(os.path.join(out_dir, 'gather', '{sample}.{{aks}}.gather.with-lineages.csv'), sample=SAMPLES)
    output: 
        os.path.join(out_dir, 'gather', f"{basename}.{{aks}}.gather-pathlist.txt")
    run:
        with open(str(output), 'w') as outF:
            for inF in input:
                outF.write(str(inF)+ "\n")

