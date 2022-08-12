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
# hacky, but useful: drop k-mer trimmed read name into sample_info dict for easier finding
sample_info['abundtrim'] = f'{out_dir}/abundtrim/' + sample_info['name'] + '.abundtrim.fq.gz'
# rename reads --> raw_reads for clarity
sample_info.rename(columns={'reads': 'raw_reads'}, inplace=True)
# set name as index for easy access
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
    print("Alas!\n")


rule all:
    input:
        ancient(expand(os.path.join(out_dir, f"{basename}.{{read_type}}.queries.zip"), read_type=['raw_reads', 'abundtrim'])),
        expand(os.path.join(out_dir, '{gather_type}', f"{basename}.{{aks}}.gather-pathlist.txt"), aks=alpha_ksize_scaled, gather_type=['abundtrim-gather', 'abund-gather'])


# k-mer abundance trimming
rule kmer_trim_reads:
    input: 
        #reads = ancient(outdir + '/trim/{sample}.trim.fq.gz'),
        reads = lambda w: sample_info.at[w.sample, "raw_reads"]
    output:
        protected(os.path.join(out_dir, "abundtrim", "{sample}.abundtrim.fq.gz"))
    conda: 'conf/env/trim.yml'
    resources:
        mem_mb = int(20e9/ 1e6),
        partition = 'bmm',
        time=240,
    params:
        mem = 20e9,
    shell: """
            trim-low-abund.py -C 3 -Z 18 -M {params.mem} -V \
            {input.reads} -o {output} --gzip
    """

rule sourmash_sketch_translate:
    #input: os.path.join(out_dir, "{read_type}", "{sample}.abundtrim.fq.gz")
   # input: lambda w: sample_info.at[w.sample, "reads"] # non-abundtrimmed sample
    input: lambda w: sample_info.at[w.sample, w.read_type] # get raw reads or abundtrim reads. read_type = 'reads' or 'abundtrim'
    output:
        os.path.join(out_dir, "{read_type}", "{sample}.translate.sig.zip")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    log: os.path.join(logs_dir, "sketch_{read_type}", "{sample}.sketch_translate.log")
    benchmark: os.path.join(benchmarks_dir, "sketch_{read_type}", "{sample}.sketch_translate.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch translate {input} -p k=7,k=10,protein,scaled=200,abund \
                                          -p k=16,k=19,dayhoff,scaled=200,abund \
                                          --name {wildcards.sample} -o {output} 2> {log}
        """

rule sourmash_sketch_dna:
    #input: os.path.join(out_dir, "abundtrim", "{sample}.abundtrim.fq.gz")
   # input: lambda w: sample_info.at[w.sample, "reads"] # non-abundtrimmed sample
    input: lambda w: sample_info.at[w.sample, w.read_type] # get raw reads or abundtrim reads. read_type = 'raw_reads' or 'abundtrim'
    output:
        os.path.join(out_dir, "{read_type}", "{sample}.dna.sig.zip")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    log: os.path.join(logs_dir, "sketch_{read_type}", "{sample}.sketch_dna.log")
    benchmark: os.path.join(benchmarks_dir, "sketch_{read_type}", "{sample}.sketch_dna.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch dna {input} -p k=21,k=31,k=51,dna,scaled=1000,abund \
                                    --name {wildcards.sample} -o {output} 2> {log}
        """

localrules: sig_cat
rule sig_cat:
    input: 
        expand(os.path.join(out_dir, "{{read_type}}", "{sample}.{sketch_type}.sig.zip"), sample=SAMPLES, sketch_type = ['translate', 'dna']),
    output:
        zipF=os.path.join(out_dir, f"{basename}.{{read_type}}.queries.zip"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    log: os.path.join(logs_dir, "sig_cat_{read_type}", f"{basename}.sigcat.log")
    benchmark: os.path.join(benchmarks_dir, "sig_cat_{read_type}", f"{basename}.sigcat.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sig cat {input} -o {output} 2> {log}
        """

# consider switching this to just prefetch - think over utility of abundtrim results vs abund results
rule gather_abundtrim_sig_from_zipfile:
    input:
        query_zip=ancient(os.path.join(out_dir, f"{basename}.abundtrim.queries.zip")),
        databases = lambda w: search_databases[f"{w.alphabet}-k{w.ksize}"]
    output:
        prefetch_csv=os.path.join(out_dir, 'abundtrim-gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.prefetch.csv'),
        gather_csv=os.path.join(out_dir, 'abundtrim-gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'abundtrim-gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.txt'),
    params:
        #threshold_bp = config['sourmash_database_threshold_bp'],
        alpha_cmd = lambda w: "--" + w.alphabet 
    resources:
        mem_mb=lambda wildcards, attempt: attempt *30000,
        time=4000,
        partition="bmm",
    log: os.path.join(logs_dir, "abundtrim-gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "abundtrim-gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB is {input.databases}"
        echo "DB is {input.databases}" > {log}

        sourmash sig grep {wildcards.sample} {input.query_zip} {params.alpha_cmd} \
                 --ksize {wildcards.ksize} | sourmash gather - {input.databases} \
                 -o {output.gather_csv} -k {wildcards.ksize} --scaled {wildcards.scaled} {params.alpha_cmd} \
                 --save-prefetch-csv {output.prefetch_csv} --ignore-abundance > {output.gather_txt} 2> {log}
        touch {output.prefetch_csv}
        touch {output.gather_txt}
        touch {output.gather_csv}
        """
                #--threshold-bp={params.threshold_bp} \ 
                 #--picklist {input.query_picklist}:name:identprefix:exclude \


rule gather_raw_read_sig_using_abundtrim_prefetch:
    input:
        query_zip=ancient(os.path.join(out_dir, f"{basename}.raw_reads.queries.zip")),
        prefetch_csv=os.path.join(out_dir, 'abundtrim-gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.prefetch.csv'),
        databases = lambda w: search_databases[f"{w.alphabet}-k{w.ksize}"]
    output:
        gather_csv=os.path.join(out_dir, 'abund-gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'abund-gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.txt'),
    params:
        #threshold_bp = config['sourmash_database_threshold_bp'],
        alpha_cmd = lambda w: "--" + w.alphabet
    resources:
        mem_mb=lambda wildcards, attempt: attempt *30000,
        time=4000,
        partition="bmm",
    log: os.path.join(logs_dir, "abundtrim-gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "abundtrim-gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB is {input.databases}"
        echo "DB is {input.databases}" > {log}

        sourmash sig grep {wildcards.sample} {input.query_zip} {params.alpha_cmd} \
                 --ksize {wildcards.ksize} | sourmash gather - {input.databases} \
                 -o {output.gather_csv} -k {wildcards.ksize} --scaled {wildcards.scaled} \
                 --picklist {input.prefetch_csv}::prefetch {params.alpha_cmd} \
                  > {output.gather_txt} 2> {log}
        touch {output.gather_txt}
        touch {output.gather_csv}
        """


rule tax_annotate:
    input:
        gather = os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.with-lineages.csv'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    params:
        lineage_cmd = "-t" + "-t".join(config['database_lineage_files'])
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash tax annotate -g {input.gather} {params.lineage_cmd}
        """


localrules: annotated_gather_csvs_to_pathlist
rule annotated_gather_csvs_to_pathlist:
    input: 
        expand(os.path.join(out_dir, '{{gather_type}}', '{sample}.{{aks}}.gather.with-lineages.csv'), sample=SAMPLES)
    output: 
        os.path.join(out_dir, '{gather_type}', f"{basename}.{{aks}}.gather-pathlist.txt")
    run:
        with open(str(output), 'w') as outF:
            for inF in input:
                outF.write(str(inF)+ "\n")

