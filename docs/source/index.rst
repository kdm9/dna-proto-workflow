.. dna-proto-workflow documentation master file, created by
   sphinx-quickstart on Mon Nov 16 17:07:19 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========================================
Mutant-Analysis-workflow - Documentation
========================================


.. toctree::
   :maxdepth: 4
   :caption: Contents:

.. |br| raw:: html

  <br/>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Overview PBGL's Mutant-Analysis-workflow
----------------------------------------

This collection of rules constitutes a workflow for the analysis of high-throughput sequencing data. We use it mostly for data produced with Illumina sequencing instruments. For this particular workflow we distinguish two main use cases: denovo and varcall.

Different Use Cases
^^^^^^^^^^^^^^^^^^^

In Snakemake, calling a rule will invoke all upstream rules. A rule will be run if:

* Snakemake realizes that the expected ``output`` files of a rule are not (yet) present or,
* Snakemake realizes that the expected input files or configuration parameters of a rule have changed since the last run.

For details, please consult the `Snakemake <https://snakemake.readthedocs.io/en/stable/index.html#>`_ documentation. For the differente use cases in this workflow, invoke the rule that produces the desired output:

1 - "De-Novo"
~~~~~~~~~~~~~

Choosing this option will conduct a reference-free comparision of *samples* based on the raw sequencing reads using k-mers. The workflow invokes the software tools kWIP and/or mash and outputs distance matrices and PCA plots. The workflow for this use case consists of the following steps:

+---+----------------------------+--------+------------------+
|#  |            Task            |  Rule  |    Software      |
+===+============================+========+==================+
| 1 | Prepare/clip the raw reads | readQC | AdapterRemoval   |
+---+----------------------------+--------+------------------+
| 2 | Calculate distances        | denovo | kWIP and/or mash |
+---+----------------------------+--------+------------------+
| 3 | Perform PCA and plot       | denovo | Utils/R-Scripts  |
+---+----------------------------+--------+------------------+


* Invoke this option by running the rule: ``rules.denovo.input``

Required input for ``rules.denovo.input`` are fastq files and the workflow will return distance matrices, produced by kWIP and/or mash.

2 - "Variant Calling"
~~~~~~~~~~~~~~~~~~~~~

Choosing this option will run a full re-sequencing analysis. It detects variants and genotype *samples* based on the alignments of the sequencing reads against one or several user-defined reference genome(s). Reads can be mapped with bwa and/or nextgenmapper (ngm), and variants called with freebayes and/or mpileup. If reference genome annotation is available, the effects of variants on gene integrity can also be predicted using the software `snpEff <https://pcingola.github.io/SnpEff/se_introduction/>`_.

The full workflow for this use case consists of the following steps:


+---+-----------------------------------------+---------+--------------------------+
| # |                 Task                    |  Rule   |        Software          |
+===+=========================================+=========+==========================+
| 1 | Prepare/clip the raw reads              | readQC  | AdapterRemoval           |
+---+-----------------------------------------+---------+--------------------------+
| 2 | Align the reads to the reference genome | align   | bwa and/or ngm           |
+---+-----------------------------------------+---------+--------------------------+
| 3 | Mark duplicates                         | align   | samtools fixmate |br|    |
|   |                                         |         | samtools sort    |br|    |
|   |                                         |         | samtools markdup         |
+---+-----------------------------------------+---------+--------------------------+
| 4 | Realign indels                          | align   | samtools merge  |br|     |
|   |                                         |         | abra2                    |
+---+-----------------------------------------+---------+--------------------------+
| 5 | Call variants                           | varcall | freebayes and/or mpileup |
+---+-----------------------------------------+---------+--------------------------+
| 6 | Filter variants                         | varcall | bcftools view            |
+---+-----------------------------------------+---------+--------------------------+
| 7 | Annotate variants                       | snpeff  | snpeff                   |
+---+-----------------------------------------+---------+--------------------------+

===== ======================================= ======= ========================
 #                     Task                    Rule          Software
===== ======================================= ======= ========================
1     Prepare/clip the raw reads              readQC  AdapterRemoval
2     Align the reads to the reference genome align   bwa and/or ngm
3     Mark duplicates                         align   samtools fixmate |br|
                                                      samtools sort |br|
                                                      samtools markdup 
4     Realign indels                          align   samtools merge |br|
                                                      abra2 
5     Call variants                           varcall freebayes and/or mpileup
6     Filter variants                         varcall bcftools view
7     Annotate variants                       snpeff  snpEff
===== ======================================= ======= ========================

This option can be invoked in 2 ways:

* running the ``rules.varcall.input`` will result in **several filterd vcf files**, one for each specified filter.
* running the rule ``rules.snpEff.output`` will result in **annotated vcf files for one chosen\filter** ``(config['snpeff']['filter'])`` and additional summaries; variants are filtered with the specified filter only and variant effects annotated against the provided snpEff database.

Required input files are fastq files and a genome reference (fasta). The rule ``snpeff`` in addition depends on a genome annotation for the reference genome used. For maximum flexibility and ease of trouble shooting we recommend to first run the re-sequencing analysis by invoking ``rules.varcall.input``, and upon successful completion invoke the workflow again, uncommenting ``rules.snpEff.output``.

3 - "snpEff"
~~~~~~~~~~~~

Choosing this option will attempt to annotate vcf files provided in **output/variants/final/** for a filter setting specified in the ``config.yml`` ``(config['snpeff']['filter'])`` and the chosen reference genome ``(config['snpeff']['name'])``. This workflow has only one step.

+---+-------------------+--------+----------+
| # |       Task        |  Rule  | Software |
+===+===================+========+==========+
| 1 | Annotate variants | snpeff | snpEff   |
+---+-------------------+--------+----------+

===== ================= ====== ========
 #          Task         Rule  Software
===== ================= ====== ========
1     Annotate variants snpeff snpEff
===== ================= ====== ========

Typical use case is to run ``snpEff`` after a completed run of rule ``varcall``. A snpEff run
will complete within a matter of minutes.

Hardware Requirements
^^^^^^^^^^^^^^^^^^^^^

The workflow is parallelized and snakemake will make efficient use of available resources on local machines as well as on compute clusters. It will run faster the more resources are available, but it performs fine on smaller machines. For routine applications we have used the workflow on a budged workstation, HP Z820 with 32 cores, 64 GB RAM running Ubuntu 16.04, and on a Virtual Machine in the cloud, AZURE with 16 cores and 512 GB RAM running Ubuntu 18.04.

Snakemake allows for fine-tuning resource allocation to the individual rules (number of processors and memory). We configured each rule with reasonable defaults, but they can be tailored to your particular size project and hardware. For details please consult the snakemake manual. In general though, if limited by memory, then do not parallelise too aggressively.

Software Dependencies
^^^^^^^^^^^^^^^^^^^^^

We recommend running the workflow in its own ``conda environment`` on a Linux Server. Dependencies are listed in ``envs/condaenv.yml`` and ``envs/additional.yml``. A brief explanation how to use these files to generate the conda environment is further below. For comprehensive explanation please consult the `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_ documentation. For software that was not available through conda at the time of development we make the specific binaries available in ``envs/``. Currently mainly ``abra2.jar``.

Workflow Use in a Nutshell
--------------------------

Steps
^^^^^

1. Git clone the workflow
2. Create and activate the conda environment
3. Provide reference genome(s) and annotation(s) in ``/genomes and annotation/``
4. Specify location of input files and their meta data in ``/metadata/sample2runlib.csv``
5. Provide lists of samples to analyse in ``/metadata/samplesets/``
6. Uncomment the respective workflow option for your use case in the ``/Snakefile``
7. Configure software parameters in ``/config.yml``
8. In case the ``snpeff`` rule will be called, adapt ``/snpeff.config``
9. Run snakemake

For standard applications no additional files need to be edited. The rules ``(*.smk)`` reside in
``rules/``. Most rules have explicit shell commands with transparent flag settings. Expert users
can change these for additional control. Note that, in snakemake, calling a rule will trigger run
of the upstream rules. It is therefore important to only configure the most downstream rule
``(/config.yml)``; these settings will be propagated to the upstream rules.

Workflow Use in Detail
----------------------

Creating the Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend running the workflow in its own conda environment on a Linux Server. The
files ``envs/condaenv.yml`` and ``envs/additional.yml`` can directly be used to create the
environment and install dependencies like so:

::

   $ conda create --name dna-proto
   $ conda env update --name dna-proto --file condaenv.yml
   $ conda env update --name dna-proto --file additional.yml

If env update does not work as intended or fails, then please ``conda install`` the programs individually specifying the correct channel and the version where required. Below we give examples, but please consult the `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_ (and `bioconda <https://bioconda.github.io/>`_) documentation.

Example:

:: 

   $ conda install -c bioconda samtools=1.9

In case of manual installs, it is convenient to add all required channels to the ``conda config``. To install the dependencies you will need below channels.

::

   $ conda config --show channels
   channels:
         - defaults
         - bioconda
         - conda-forge
         - r   
   $ conda config --add channels bioconda
   $ conda install samtools=1.9

The required channels are also listed in the respective ``.yml`` files. Configuring channels has the pitfall of rare ambiguities and collisions. Please consult the `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_ documentation for “managing channels”.

Reference Genome
^^^^^^^^^^^^^^^^

The workflow will look for the reference genomes and associated files in ``genomes_and_annotations/``. We provide an example directory tree in ``genomes_annotatoins/readme``. We recommend creating one subdirectory for each reference genome. Each reference genome directory must contain the necessary assembly file ``(.fa/.fna)`` and the associated files needed by the aligners. Generate the associated files in this directory from the fasta file like so:

::

   $ samtools faidx <reference-genome.fa>
   $ bwa index -a bwtsw <reference-genome.fa>
   $ ngm -r <reference-genome.fa>

Reference Genome Annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to annotate the variants (using snpEff) you will need a snpEff database (``snpEffectPredictor.bin``) for the relevant reference genome. The location of the ``snpeff`` database can be configured in ``snpeff.config``. It is currently set to ``genomes_and_annotations/snpeffdata/`` and again, we recommend subdirectories for the reference genomes in this directory. Appending ``_snpeff`` to these directory names will help avoid confusion. In order to create a ``snpeff`` database this ``reference_genome_snpeff`` directory must contain the reference genome and the annotation in files named ``sequences.fa`` and ``genes.gff``; ``protein.fa`` and other files are optional. You will also need to add the respective entry in ``snpeff.config``. Then, from the root directory of the workflow (the directory that contains the ``snpeff.config`` file) execute:

::

   $ snpEff build –gff3 genomes_and_annotations/snpeffdata/<reference-genome>_snpeff

While a ``.gtf`` file can also be used (-gtf 22), we have better experience building databases from ``.gff`` files. For a detailed explanation of the ``snpeff`` build process please consult the `snpEff <https://pcingola.github.io/SnpEff/se_introduction/>`_ documentation.

Example directory tree of ``genomes_and_annotations/`` for a cowpea reference genome downloaded from NCBI. Notice that we store the reference geomes somewhere else and use soft links. Compare also to our entries in ``snpeff.config`` under “Non-standard Databases” and replicate accordingly. The database ``snpEffectPredictor.bin`` will be generated by ``snpEff build``.

::

   genomes_and_annotations/
   ├── GCF_004118075.1_ASM411807v1 -> GCF_004118075.1_ASM411807v1/
   ├── readme
   └── snpeffdata
       └── GCF_004118075.1_ASM411807v1_snpeff
           ├── genes.gff -> GCF_004118075.1_ASM411807v1/GCF_004118075.1_ASM411807v1_genomic.gff
           ├── genes.gtf -> GCF_004118075.1_ASM411807v1/GCF_004118075.1_ASM411807v1_genomic.gtf
           ├── protein.fa -> GCF_004118075.1_ASM411807v1/GCF_004118075.1_ASM411807v1_protein.faa
           ├── sequences.fa -> GCF_004118075.1_ASM411807v1/GCF_004118075.1_ASM411807v1_genomic.fna
           └── snpEffectPredictor.bin

Providing Required Materials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Samplesets
~~~~~~~~~~

We use the term “sample” as the entity of interest. Our workflows compare “samples” with one another. Lists of sample names must be provided as text files (``.txt``) in ``/metadata/samplesets/`` and those files can be referred to as “samplesets” in ``config.yml``. There exists a glob for all samples: ‘all_samples’, which will have the effect of concatenating all ``*.txt`` files in ``/metadata/samplesets/``. The intention is to enable easy addition or removal of samples to/from an existing analysis.

Sample Definitions
~~~~~~~~~~~~~~~~~~

In praxis, oftentimes DNA gets extracted from an individual, turned into one or several sequencing libraries that are then sequenced in one or several sequencing runs. Obviously, if then those individuals are to be compared, all respective fastq files need to be assigned to the same sample. Our workflow accommodates for this, but the samples need to be defined. The sample definitions, are provided in ``/metadata/sample2runlib.csv`` and make explicit, how the fastq files constitute the samples:

The entries in columns “run” and “library” are together used as the primary key, i.e., as the unique identifier for the “sequencing run” and thus the fastq file(s). “Sequencing runs” are assigned to the same “sample” through an identical sample name in the “sample” column. fastq files can be provided either as separate forward and reverse read files (fq1, fq2) or as interleaved fastq (il_fq), with the respective other column(s) empty. Within one ``sample2runlib.csv`` file interleaved and two-file input can be mixed.

+------+---------+--------+----------------+----------------+----------------+
| run  | library | sample |      fq_1      |     fq_2       |     if_fq      |
+======+=========+========+================+================+================+
| Run1 | A-500bp | A      | <path to file> | <path to file> |                |
+------+---------+--------+----------------+----------------+----------------+
| Run1 | A-300bp | A      | <path to file> | <path to file> |                |
+------+---------+--------+----------------+----------------+----------------+
| Run2 | A-500bp | A      |                |                | <path to file> |
+------+---------+--------+----------------+----------------+----------------+
| Run2 | B-300bp | B      |                |                | <path to file> |
+------+---------+--------+----------------+----------------+----------------+
| Run1 | B-500bp | B      |                |                | <path to file> |
+------+---------+--------+----------------+----------------+----------------+
| Run1 | C-500bp | C      | <path to file> | <path to file> |                |
+------+---------+--------+----------------+----------------+----------------+

==== ======= ====== ============== ============== ==============
run  library sample      fq_1           fq_2          if_fq          
==== ======= ====== ============== ============== ==============
Run1 A-500bp A      <path to file> <path to file>
Run1 A-300bp A      <path to file> <path to file>
Run2 A-500bp A                                    <path to file>
Run2 B-300bp B                                    <path to file>
Run1 B-500bp B                                    <path to file>
Run1 C-500bp C      <path to file> <path to file>
==== ======= ====== ============== ============== ==============

Regions of Interest for Variant Calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variant calling can be restricted to particular regions of interest by providing them in ``/metadata/contigs_of_interest.bed``: A file in bed file format listing identifier, start-, and end- positions (tab delimited, no header). This is helpful for exome capture data but also to restrict the analysis to specific chromosomes. Genome assemblies often contain the main chromosomes and in addition many orphan fragments that often are not of interest. Below example will restrict variant calling to the 11 chromosomes plus the chloroplast of cowpea. The hundreds of additional contigs in the cowpea reference genome are available for read mapping but variants will not be called in them and they will then not be present in the vcf file. Note that lines starting with ‘``#``’ will be disregarded.

::

   NC_040279.1 0 42129361
   NC_040280.1 0 33908088
   NC_040281.1 0 65292630
   NC_040282.1 0 42731077
   NC_040283.1 0 48746289
   NC_040284.1 0 34463471
   NC_040285.1 0 40876636
   NC_040286.1 0 38363498
   NC_040287.1 0 43933251
   NC_040288.1 0 41327797
   NC_040289.1 0 41684185
   NC_018051.1 0 152415

Workflow Configuration
^^^^^^^^^^^^^^^^^^^^^^

config.yml
~~~~~~~~~~~

Central place for the configuration of workflow behavior and software parameters by the user is ``config.yml``. There are comments in the file that explain the configuration parameters and options.

An important configuration parameter is the location of a ``tmp/`` directory. Several rules make extensive use of the tmp/ directory to temporarily store large files. Oftentimes standard home directories on compute servers or cluster nodes are too small.

Additional Configuration
~~~~~~~~~~~~~~~~~~~~~~~~

Expert users can change the rules themselves by editing ``rules/*.rules.smk``. Use caution! We have chosen reasonable defaults and recommend modifying rules only when you know what you are doing. When allocation more cores to rules, pay attention that a) some rules are very memory intensive b) some shell commands are piped and are in fact using more than 1 core per process.

Running the Workflow
~~~~~~~~~~~~~~~~~~~~

To run the workflow, un-comment the respective rule in the ``Snakefile`` and run snakemake.

:: 
   
   $ snakemake –npr
   $ snakemake –j 6 --no-temp -kpr

For details on commandline options for snakemake please consult the snakemake manual. Un-comment only the one most downstream rule for your use case. Currently, these use case rules are

::

   # USER OPTIONS
   #     rules.denovo.input,
   #     rules.varcall.input,
   #     rules.snpeff.output,

A rule that encounters missing input files will invoke the respective upstream rule(s). E.g., if “``rules.snpeff.output``” is uncommented and snakemake is run for the first time, the entire workflow from readqc, alignment (= read alignments, duplicate removal, indel realignment), varcall (=variant calling), and ``snpeff`` (=variant annotation) will run in one go. In case no ``snpeff`` database is supplied, then rule ``rules.snpeff.output`` cannot be run. Then use ``rules.varcall.input``.

When configuring ``config.yml``, keep in mind that configuration parameters of a downstream rule take precedence because parameters will propagate upstream. I.e., you must set your alignment parameters in the ``varcall:`` section. Adjusting parameters under ``align:`` will not have the intended effect. Only in special circumstances will the “rules.align.input” be run by itself and only then will you have to adjust parameters in the ``align:`` section.

The entities that the workflows compare are “samples” as listed in ``samplesets/*.txt`` and defined in ``/metadata/sample2runlib.csv``. Sample names listed in ``samplesets/*.txt`` must correspond to the entries in the sample column in ``/metadata/sample2runlib.csv``.

Workflow Use Cases
^^^^^^^^^^^^^^^^^^

De-Novo Analysis
~~~~~~~~~~~~~~~~

De-novo analysis can be used to check the relatedness of sequencing runs and/or samples without the use of any reference genome. It will perform a comparative analysis based on kmers on the raw data and output distance matrices. We recommend performing a de-novo analysis at the very start of every project at the “sequencing run”-level prior to any merging of runs into samples. This can help to detect mix-ups and mislabels and is achieved by providing unique names for each sequencing run in the sample column of
``/metadata/sample2runlib.csv``. De-novo analysis is invoked by calling “``rules.denovo.input``” by uncomment the respective line in the Snakefile (and only this line). Pay attention to maintain the indentation.

::

   rule all:
       input:
   # USER OPTIONS
               rules.denovo.input,
   #            rules.varcall.input,
   #            rules.snpeff.output,
   # EXPERT OPTIONS
   #            rules.readqc.input,
   #            rules.align.input,
   #            rules.stats.input,

Variant Calling - Standard Re-Sequencing Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running “``rules.varcall.input``” will call variants and genotype “samples” with respect to one or several reference genomes. Varcall will compare samples listed in samplesets with one another. , as defined in the sample column. Invoke this analysis by uncommenting “``rules.varcall.input``” (and only this line). Pay attention to maintain the indentation.

::

   rule all:
       input:
   # USER OPTIONS
   #        rules.denovo.input,
           rules.varcall.input,
   #        rules.snpeff.output,
   # EXPERT OPTIONS
   #        rules.readqc.input,
   #        rules.align.input,
   #        rules.stats.input,

Variant Annotation - The Effects of Variants on Gene Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a ``snpeff`` library based on this exact reference genome is provided then the entire workflow can in principle be invoke by uncommenting “``rules.snpeff.output``”. (and only this line). Pay attention to maintain the indentation.

::

   rule all:
       input:
   # USER OPTIONS
   #       rules.denovo.input,
   #       rules.varcall.input,
          rules.snpeff.output,
   # EXPERT OPTIONS
   #       rules.readqc.input,
   #       rules.align.input,
   #       rules.stats.input,

There is the subtle difference that the rule ``snpeff`` will operate only on one single vcf file. The ``snpeff`` rule will hence propagate upstream as requirement only one reference genome and one filter setting. If variants with several different filter settings and/or against several reference genomes are desired, then it is advantageous to first run the varcall rule and subsequently the ``snpeff`` rule.

Workflow Outputs
----------------

After completion of the run, all output, including logs and stats, will be in ``output/``.

* Clipped reads (in interleaved fastq) format are in ``output/reads/``
* BAM files with In/Del-realigned alignments are in ``output/abra/``
* BCF/VCF files of the filtered variants including respective index files are in ``output/variants/final/``
* The snpEff-annotated variant file is ``output/snpeff/annotated/all.vcf.gz``

For loading into IGV, use the In/Del realigned BAM file in ``output/abra/`` and the ``*.vcf.gz`` files of the filtered variants. Note that IGV requires the ``vcf.gz.tbi`` index.

Support
-------

License
^^^^^^^

Mutant-Analysis-workflow is Copyright 2020 Norman Warthmann, and released under the GNU General Public License version 3 (or any later version).

Contributors
^^^^^^^^^^^^

This workflow was developed by Norman Warthmann, PBGL, with important contributions from Kevin D Murray (Australian National University) and Marcos Conde, PBGL. The documentation was written by Norman Warthmann with contribution from Anibal Morales, PBGL.

Reproducibility
^^^^^^^^^^^^^^^

Workflows help addressing reproducibility issues. Consider making your version of the workflow, configured for your data, available upon publication of your results.

Getting Help
^^^^^^^^^^^^

Feel free to let us know if you are using our workflow and don’t hesitate to contact us with questions: email **n.warthmann@iaea.org**
