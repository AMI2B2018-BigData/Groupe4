################################
####### mapping avec BWA #######
################################
# Tous vos fichers fastq doivent être dans raw_data/xxx.fq.gz
# 
# L'extension doit être exactement fq.gz
# 
# Il doit y avoir hg38.fa dans ./
# 
# Pour des données paired-end (2 fichiers : _R1/_R2): 
# input: 
#     hg38.fa
#     raw_data/sampleX_c3-pe_R1.fq.gz
#     raw_data/sampleX_c3-pe_R2.fq.gz
# output:
#     mapping/hg38/sampleX_c3-pe.sort.bam
# 
# processus : 
#     Index du génome de référence
#     Mapping de deux jeux de données (se_7 et c3-pe)(bwa)
#         Alignement sur l'index du génome (bwa aln)
#         Transformation des coordonnées de l'index en coordonnées 
#             génomiques (bwa sampe/samse)
#     Rangement des données par coordonnées génomique (samtools sort)
#     Indexer le fichier ordonnée (samtools index)
#
################################
################################ 

# Pour récupérer la version des outils utilisés (voir rule sort/index)
import subprocess

# Version de bwa
bwa_version = "0.7.10"

# Fichiers en sortie du pipeline ("hg38.fa.bwt" aurait pu être omis)
rule target:
    input:
        "hg38.fa.bwt",
        "mapping/hg38/se_7.sort.bam.bai",
        "mapping/hg38/c3-pe.sort.bam.bai"

# Index du génome de référence avant l'alignement des lectures (séquences)
rule bwa_index:
    input:
        "{genome}.fa" #2) cherche {wildcards.genome}.fa"
    output:
        "{genome}.fa.bwt" #1) crée la variable wildcards.genome.
    message:
        "Index du génome de référence :{wildcards.genome}"
    version: bwa_version
    threads: 1  
    shell: #3) lance la commande
        "bwa index {input}" #"{input}"="{wildcards.genome}.fa"

# Alignement des lectures sur le génome de référence
rule bwa_aln:
    input:
        bwt = "{genome}.fa.bwt",
        fasta = "{genome}.fa",
        fq = "raw_data/{samples}.fq.gz"
    output:
        temp("mapping/{genome}/{samples}.sai") #fichier temporaire
    threads: 20 #max cpu autorisé
    version: bwa_version
    log:
        "mapping/{genome}/{samples}.bwa_aln.log"
    message:
        "Alignment de {wildcards.samples} sur {wildcards.genome} "
    shell:
        "bwa aln "       # avec Snakemake
        "-t {threads} "  # vous pouvez 
        "{input.fasta} " # commenter
        "{input.fq} "    # tous les 
        ">{output} "     # arguments 
        "2>{log}"        # des commandes

# Si les données sont appariées (paired-end), autrement dit si les fragments 
# ont été séquencés dans les deux sens, cette règle sera utilisée
rule bwa_sampe_to_bam:
    input:
        fasta = "{genome}.fa",
        sai1 = "mapping/{genome}/{samples}_R1.sai",
        fq1 = "raw_data/{samples}_R1.fq.gz",
        sai2 = "mapping/{genome}/{samples}_R2.sai",
        fq2 = "raw_data/{samples}_R2.fq.gz"
    output:
        #fichier temporaire
        temp("mapping/{genome}/{samples}.bam")
    log:
        "mapping/{genome}/{samples}.samse.log"
    threads: 2
    version: bwa_version
    message:
        "PE : sai --> bam ({wildcards.genome}/{wildcards.samples})"
    shell:
        "bwa sampe "
        "{input.fasta} "
        "{input.sai1} {input.sai2} "
        "{input.fq1} {input.fq2} "
        "2>{log} | "
        "samtools view -Sbh - > {output}"

# Si on a séquencé les fragments dans un sens (Single read) cette règle sera utilisée
rule bwa_samse_to_bam:
    input:
        fasta = "{genome}.fa",
        sai = "mapping/{genome}/{samples}.sai",
        fq = "raw_data/{samples}.fq.gz"
    output:
        #fichier temporaire
        temp("mapping/{genome}/{samples}.bam")
    log:
        "mapping/{genome}/{samples}.samse.log"
    threads: 2
    message:
        "SE : sai --> bam ({wildcards.genome}/{wildcards.samples})"
    shell:
        "bwa samse {input.fasta} {input.sai} {input.fq} "
        "2>{log} | "
        "samtools view -Sbh - > {output}"

# Pour ranger les données par coordonnées génomiques
rule sort_bam:
    input:
        "{base}.bam"
    output:
        # fichier protégé en écriture
        protected("{base}.sort.bam")
    threads: 10
    log:
        "{base}.sort.log"
    message:
        "Ordonne le fichier {wildcards.base}.bam"
    version: # pour récupérer la version de l'outil avec une commande shell 
        subprocess.getoutput(
            "samtools --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    params:
        niveau_compression = "9",
        memoire_max_par_cpu = "1G"
    shell:
        "samtools sort "
        "-l {params.niveau_compression} "
        "-m {params.memoire_max_par_cpu} "
        "-@ {threads} "
        "-f {input} "
        "{output} "
        "2>{log}"

# Pour indexer le fichier trié (balise le fichier pour accéder aux données plus rapidement)
rule index_bam:
    input:
        "{base}.sort.bam"
    output:
        "{base}.sort.bam.bai"
    priority: 5 # augmente la priorité de la règle (défaut = 1)
    version: 
        subprocess.getoutput(
            "samtools --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    threads: 1
    shell:
        "samtools index {input}"