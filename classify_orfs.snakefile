import os, sys
import pandas as pd
import glob

#configfile: "inputs/spades_orfs.classify.conf"
configfile: "inputs/metaflye_orfs.classify.conf"

out_dir = config.get('output_dir', 'output.classify_orfs')
logs_dir = os.path.join(out_dir, "logs")
basename = config.get('basename', 'np-orfs')
sample_csv = config['sample_info']
sample_info = pd.read_csv(sample_csv)

SAMPLES = sample_info["name"].to_list()
search_databases = config['search_databases']

if not search_databases:
    print('Must provide one or more search databases in your config.')
    sys.exit(-1)

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
    alpha_ksize_scaled += expand(f"{alpha}-k{{ksize}}-sc{{scaled}}", ksize = ksize, scaled=scaled)
    # build param string for sketching
    if alpha in ["nucleotide", "dna", "rna"]:
        param_str = '-p dna,'
    elif alpha == "protein": 
        param_str = '-p protein,'
    elif alpha == "dayhoff": 
        param_str = '-p dayhoff,'
    ks = [ f'k={k}' for k in ksize ]
    ks = ",".join(ks)
    sketch_scaled= min(scaled)
    param_str += f"{ks},scaled={sketch_scaled},abund"
    all_param_str.append(param_str)

fromfile_param_strs = ' '.join(all_param_str)

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
        ancient(os.path.join(out_dir, f"{basename}.queries.zip")),
        expand(os.path.join(out_dir, 'gather', f"{basename}.{{aks}}.gather.lineage_summary.tsv"), aks = alpha_ksize_scaled),
        expand(os.path.join(out_dir, 'gather', '{sample}.{aks}.gather.{ext}'), sample=SAMPLES, aks = alpha_ksize_scaled, ext = ['summarized.csv', 'krona.tsv']),
    #    expand(os.path.join(out_dir, 'gather', f'{{sample}}.{sketch_type}.{alphabet}-k{{ksize}}.gather.csv'),ksize=ksize, sample=sample_names)
    #    expand(os.path.join(out_dir, 'classify', f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classifications.csv"), ksize = ksize),


rule sourmash_sketch_fromfile:
    input: config["sample_info"]
    output:
        zipF=os.path.join(out_dir, f"{basename}.queries.zip"),
    params:
        param_str = fromfile_param_strs
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sketchfromfile", f"{basename}.queries.sketch.log")
    benchmark: os.path.join(logs_dir, "sketchfromfile", f"{basename}.queries.sketch.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch fromfile {input} {params.param_str} -o {output.zipF} \
                                 --force-output-already-exists 2> {log}
        """
                                 #--already-done {output.zipF} \

rule gather_sig_from_zipfile:
    input:
        query_zip=ancient(os.path.join(out_dir, f"{basename}.queries.zip")),
        databases = lambda w: search_databases[f"{w.alphabet}-k{w.ksize}"],
    output:
        prefetch_csv=os.path.join(out_dir, 'prefetch', '{sample}.{alphabet}-k{ksize}-sc{scaled}.prefetch.csv'),
        gather_csv=os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.txt'),
    params:
        #threshold_bp = config['sourmash_database_threshold_bp'],
        alpha_cmd = lambda w: "--" + w.alphabet
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        time=4000,
        partition="med2",
    log: os.path.join(logs_dir, "gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.log")
    benchmark: os.path.join(logs_dir, "gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB is {input.databases}"
        echo "DB is {input.databases}" > {log}

        sourmash sig grep {wildcards.sample} {input.query_zip} {params.alpha_cmd} \
                 --ksize {wildcards.ksize} | sourmash gather - {input.databases} \
                 -o {output.gather_csv} -k {wildcards.ksize} --scaled {wildcards.scaled} {params.alpha_cmd} \
                 --save-prefetch-csv {output.prefetch_csv} --threshold-bp 0 > {output.gather_txt} 2> {log}
        touch {output.prefetch_csv}
        touch {output.gather_txt}
        touch {output.gather_csv}
        """

rule tax_annotate:
    input:
        gather = os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.with-lineages.csv'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, 'gather'),
    conda: "conf/env/sourmash.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax annotate -g {input.gather} -t {input.lineages} -o {params.outd}
        """

rule tax_metagenome:
    input:
        gather = os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.krona.tsv'),
        os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.summarized.csv'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'gather'),
        out_base= lambda w: f'{w.sample}.{w.alphabet}-k{w.ksize}-sc{w.scaled}.gather',
    conda: "conf/env/sourmash.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome -g {input.gather} -t {input.lineages} \
                                -o {params.out_base} --output-dir {params.outd} \
                                --output-format krona csv_summary --rank species
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


rule tax_metagenome_lineage_summary:
    input:
        gather_pathlist = os.path.join(out_dir, 'gather', f"{basename}.{{aks}}.gather-pathlist.txt"),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, 'gather', f"{basename}.{{aks}}.gather.lineage_summary.tsv")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, 'gather'),
        out_base= lambda w: f'{basename}.{w.aks}.gather',
    conda: "conf/env/sourmash.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome --from-file {input.gather_pathlist} -t {input.lineages} \
                                -o {params.out_base} --output-dir {params.outd} \
                                --output-format lineage_summary --rank species
        """


