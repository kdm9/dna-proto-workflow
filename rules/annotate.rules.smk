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

## experimental for several snpeff libraries
#expand("output/annotated_variants/snpeff/{genome[database]}/#{caller}~{aligner}~{genome[ref]}~{sampleset}~filtered-{filter}~snpEff.vcf.gz",
#                genome=config["snpeff"]["genomes"].values(),
#                caller=config["snpeff"]["callers"],
#                aligner=config["snpeff"]["aligners"],
#                sampleset=config["snpeff"]["samplesets"],
#                filter = config['snpeff']['filters']),


############### Actual Rules #################

rule snpeff:
    input:
       vcf="output/variants/final/{caller}~{aligner}~{ref}~{sampleset}~filtered-{filter}.vcf.gz",
       bed="metadata/contigs_of_interest.bed",
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
        " -filterInterval {input.bed}"
        " -csvStats {output.csvstats}"
        " -htmlStats {output.htmlStats}"
        " {params.extra}"
        " {params.database}"
        " {input.vcf}"
        " > {output.vcf}"
        " ) >'{log}' 2>&1"

############################################
