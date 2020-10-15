#######################################################################
#                           Variant Calling                           #
#######################################################################


##### Target rules #####

def raw_variant_calls_input(wildcards):
    inputs = []
    for caller in config["varcall"]["callers"]:
        for aligner in config["varcall"]["aligners"]:
            for sampleset in config["varcall"]["samplesets"]:
                for ref in config["varcall"]["refs"]:
                    this_rawfiles = expand("output/variants/raw_split/{caller}~{aligner}~{ref}~{sampleset}/{region}.bcf",
                                           caller=caller, aligner=aligner, ref=ref, sampleset=sampleset, region=VARCALL_REGIONS[ref])
                    inputs.extend(this_rawfiles)
    return inputs


rule raw_variant_calls:
    input: raw_variant_calls_input

rule filtered_variants:
    input:
        expand("output/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.{ext}",
               ext=["bcf", "bcf.csi", "vcf.gz", "vcf.gz.csi"],
               caller=config["varcall"]["callers"],
               aligner=config["varcall"]["aligners"],
               ref=config["varcall"]["refs"],
               sampleset=config["varcall"]["samplesets"],
               filter=config["varcall"]["filters"]),

rule varcall:
    input:
        rules.filtered_variants.input,


##### Actual rules #####


#### EXAMPLE

rule freebayes:
    input:
        bam = "output/abra/{aligner}~{ref}~{sampleset}.bam",
        bai = "output/abra/{aligner}~{ref}~{sampleset}.bam.bai",
        sset="output/samplelists/{sampleset}.txt",
        #sset="output/samplelists/cohort.txt",
        ref=lambda wc: config['refs'][wc.ref],
    output:
        bcf="output/variants/raw_split/freebayes~{aligner}~{ref}~{sampleset}/{region}.bcf",
    log:
        "output/log/varcall/freebayes/{aligner}~{ref}~{sampleset}/{region}.log"
    benchmark:
        "output/log/varcall/freebayes/{aligner}~{ref}~{sampleset}/{region}.benchmark"
    priority: 1  # get them done earlier, normalisation is super quick
    params:
        theta=config["varcall"].get("theta_prior", 0.01),
        minmq=lambda wc: config["varcall"]["minmapq"].get(wc.aligner, 5),
        minbq=config["varcall"]["minbq"],
    shell:

###########################

pindel_types = ["D", "BP", "INV", "TD", "LI", "SI", "RP"]

rule pindel:
    input:
        ref="genome.fasta",
        # samples to call
        samples=["mapped/a.bam"],
        # bam configuration file, see http://gmt.genome.wustl.edu/packages/pindel/quick-start.html
        config="pindel_config.txt"
    output:
        expand("pindel/all_{type}", type=pindel_types)
    params:
        # prefix must be consistent with output files
        prefix="pindel/all",
        extra=""  # optional parameters (except -i, -f, -o)
    log:
        "logs/pindel.log"
    threads: 4
    wrapper:
        "0.42.0/bio/pindel/call"

rule pindel2vcf:
    input:
        ref="genome.fasta",
        pindel="pindel/all_{type}"
    output:
        "pindel/all_{type}.vcf"
    params:
        refname="hg38",  # mandatory, see pindel manual
        refdate="20170110",  # mandatory, see pindel manual
        extra=""  # extra params (except -r, -p, -R, -d, -v)
    log:
        "logs/pindel/pindel2vcf.{type}.log"
    wrapper:
        "0.42.0/bio/pindel/pindel2vcf"

rule pindel2vcf_multi_input:
    input:
        ref="genome.fasta",
        pindel=["pindel/all_D", "pindel/all_INV"]
    output:
        "pindel/all.vcf"
    params:
        refname="hg38",  # mandatory, see pindel manual
        refdate="20170110",  # mandatory, see pindel manual
        extra=""  # extra params (except -r, -p, -R, -d, -v)
    log:
        "logs/pindel/pindel2vcf.log"
    wrapper:
        "0.42.0/bio/pindel/pindel2vcf"

