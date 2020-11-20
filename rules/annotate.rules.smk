#######################################################################
#                            Annotation                               #
#######################################################################

##### Target rules #####

rule annotate:
    input:
        expand("output/annotated_variants/snpeff/{database}/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}~snpEff.vcf.gz",
               database=config["snpeff"]['database'],
               ref=config["snpeff"]["ref"],
               caller=config["snpeff"]["callers"],
               aligner=config["snpeff"]["aligners"],
               sampleset=config["snpeff"]["samplesets"],
               filter = config['snpeff']['filters']),


############### Actual Rules #################

rule snpeff:
    input:
       "output/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf.gz"
    output:
        vcf="output/annotated_variants/snpeff/{database}/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}~snpEff.vcf.gz",
        csvstats="output/annotated_variants/snpeff/{database}/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.csv",
        htmlStats="output/annotated_variants/snpeff/{database}/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.html",
    log:
        "output/log/snpeff/{database}/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}~snpeff.log"
    params:
        database=config["snpeff"]['database'],
        extra="-Xmx6g -v"
    shell:
        "( snpEff ann"
        " -csvStats {output.csvstats}"
        " -htmlStats {output.htmlStats}"
        " {params.extra}"
        " {params.database}"
        " {input}"
        " > {output.vcf}"
        " ) >'{log}' 2>&1"

############################################

#rule prepare_database:
#    input:
#        name= config['snpeff']['database'],
#    output:
#        dir(expand("genomes_and_annotations/snpeffdata/{name}", name= config['snpeff']['name']))
#    shell:
#        "mkdir genomes_and_annotations/snpeffdata/{input.name}"
