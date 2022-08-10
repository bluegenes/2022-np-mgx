"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -s orf_classify.snakefile
"""
# classify orfs using protein sourmash gather

import pandas as pd

configfile: "inputs/spades_orfs.classify.conf"

sample_info = pd.read_csv(config['sample_info'], sep = ',')
sample_column = config.get('sample_colname', 'sample')

SAMPLES = sample_info[sample_column].to_list()
SAMPLE_FILES = sample_info["assembly"].to_list()
data_dir = config['data_dir']
out_dir = config.get('output_dir', 'output.orf_classify')
logs_dir = os.path.join(out_dir, "logs")

thumper_dirname= config.get("thumper_dir", "thumper")
#reports_dirname= config.get("reports_dir", "reports")
thumper_dir = os.path.join(out_dir, thumper_dirname)
#reports_dir = os.path.join(out_dir, reports_dirname)
basename = config.get('basename', 'spades_orfs')
#database_dir = 'databases'
database_dir = '/home/ntpierce/2021-tax/databases/'

search_databases = config.get('search_databases', ['gtdb-rs207-reps'])
input_type = config.get('input_type', 'protein')
alphabet = config.get('alphabet', 'protein')
scaled = config.get('scaled', 200)
ksizes = config.get('ksizes', [10])
sketch_type= "nucleotide" if alphabet in ['nucleotide','dna', 'rna' ] else "protein"
threshold_bp = config.get('threshold_bp', '50000')

additional_db_info = config.get('thumper_user_database_info')

sample_info.set_index(sample_column, inplace=True)


rule all:
    input: expand(os.path.join(thumper_dir, 'classify', f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classifications.csv"), ksize=ksizes),
    



rule write_protein_sample_info:
    input: expand(os.path.join(data_dir, "{assembly_file}"), assembly_file= SAMPLE_FILES)
    output:
        sample_info=os.path.join(out_dir, f"{basename}.protein.thumper.csv")
    run:
        with open(str(output), "w") as outF:
            outF.write("name,genome_filename,protein_filename\n")
            for sample_name in SAMPLES:
                fn = sample_info.at[sample_name, "assembly"] 
                full_filename = os.path.abspath(os.path.join(data_dir, fn))
                outF.write(f'{sample_name},,{full_filename}\n')
                #outF.write(f'{sample_name},{fn}\n')


rule write_thumper_config:
    input:
        sample_info=os.path.join(out_dir, f"{basename}.{alphabet}.thumper.csv")
    output:
        config= os.path.join(out_dir, f"{basename}.{alphabet}.thumper.conf")
    params:
        ksize=ksizes,
    run:
        with open(str(output), 'w') as out:
            out.write(f"strict: 1\n")
            out.write(f"force: 0\n")
            out.write(f"input_type: {input_type}\n")
            out.write(f"data_dir: {data_dir}\n")
            out.write(f"output_dir: {thumper_dir}\n")
            out.write(f"basename: {basename}\n")
            out.write(f"sample_info: {input.sample_info}\n")
            out.write(f"database_dir: {database_dir}\n")
            out.write(f"sourmash_database_threshold_bp: {threshold_bp}\n")
            out.write(f"search_databases:\n")
            for sd in search_databases:
                out.write(f"  - {sd}\n")
            if additional_db_info is not None:
                out.write(f"user_database_info: {additional_db_info}\n")
            out.write(f"alphabet: {alphabet}\n")
            out.write(f"ksize:\n")
            for k in params.ksize:
                out.write(f"  - {k}\n")
            out.write(f"search_scaled: {scaled}\n")


rule thumper_classify:
    input:
        config=rules.write_thumper_config.output.config
    output:
         expand(os.path.join(thumper_dir, 'classify', f"{basename}.{sketch_type}.{alphabet}-k{{ksize}}.classifications.csv"),ksize=ksizes)
    threads: 30
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
        time=1200,
        partition="med2",
    log: os.path.join(logs_dir, "thumper_classify", f"{basename}.{sketch_type}.{alphabet}.log")
    benchmark: os.path.join(logs_dir, "thumper_classify", f"{basename}.{sketch_type}.{alphabet}.benchmark")
    conda: "conf/env/thumper.yml"
    shell:
        """
        thumper run {input.config} --jobs {threads} genome_classify --nolock --profile slurm
        """
