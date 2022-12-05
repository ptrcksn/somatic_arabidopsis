import os


configfile: "config.json"

chr_params = config["chr_params"]

data_dir = "data"
genome_dir = "genome"
STAR_genome_dir = "STAR"
expression_dir = "expression"
pileup_dir = "pileup"
mutations_map_dir = "mutations_map"
misc_dir = "misc"
replication_timing_dir = "replication_timing"
histone_mark_dir = "histone_mark"

genome_fa = "genome.fa"
genome_gtf = "genes.gtf"
germline_variants_vcf = "arabidopsis_2029_NoFilters.vcf"
blacklist_bed = "blacklist.bed"
tair9_genes_gff = "TAIR9_GFF3_genes.gff"
rna_edit_xls = "rna_edit.xls"
tair10_genes_gff = "TAIR10_GFF3_genes.gff"
accessions_tsv = "accessions.tsv"
weigel_raw_variants_csv = "weigel.raw_variants.csv"
late_replication_bg = "L_ratio_3.smooth.bedgraph"
early_replication_bg = "E_ratio_3.smooth.bedgraph"
results_dir = "results"


with open(data_dir+"/"+misc_dir+"/"+accessions_tsv, "r") as infile:
    accessions = [ a.strip() for a in infile.readlines() ]

chrs = [ str(c) for c in [1,2,3,4,5] ]

histone_mark_replicates = [ r[:-3] for r in os.listdir(data_dir+"/"+histone_mark_dir) if r[-3:] == ".bw" ]
histone_marks = ["ATAC", "meDIP", "H3K14ac", "H3K23ac", "H3K27ac", "H3K27me1", "H3K27me3", "H3K36ac", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K56ac", "H3K9ac", "H3K9me1", "H3K9me2", "H4K16ac"]

studies = ["epigenomes"]*6+["EU_som"]*6
refs = (["C"]*3+["T"]*3)*2
alts = ["A", "G", "T", "A", "C", "G"]*2


rule all:
    input:
        merged_model_fit_tsv = data_dir+"/"+results_dir+"/all.model_fit.tsv"

rule genome_build:
    input:
        genome_fa = data_dir+"/"+genome_dir+"/"+genome_fa,
        genome_gtf = data_dir+"/"+genome_dir+"/"+genome_gtf

    output:
        STAR_genome_dir = data_dir+"/"+STAR_genome_dir

    conda:
        "envs/STAR_env.yamlSTAR_env.yaml"

    shell:
       '''
        #mkdir -p {output.STAR_genome_dir}
        STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {output.STAR_genome_dir} --genomeFastaFiles {input.genome_fa} --sjdbGTFfile {input.genome_gtf} --sjdbOverhang 99 --genomeSAindexNbases 12
        '''

rule expression_fastq_align:
    input:
        STAR_genome_dir = data_dir+"/"+STAR_genome_dir,
        expression_fastq = data_dir+"/"+expression_dir+"/{accession}.expression_fastq.gz"

    params:
        expression_dir = data_dir+"/"+expression_dir,
        accession = data_dir+"/"+expression_dir+"/{accession}"

    output:
        expression_bam = data_dir+"/"+expression_dir+"/{accession}.Aligned.out.bam",
        STAR_splice_junction_tab = data_dir+"/"+expression_dir+"/{accession}.SJ.out.tab"

    conda:
        "envs/STAR_env.yaml"

    shell:
        '''
        #mkdir -p {params.expression_dir}
        STAR --runThreadN 16 --genomeDir {input.STAR_genome_dir} --readFilesCommand gunzip -c --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outSAMattributes MD NH --outFilterMismatchNoverLmax 0.1 --clip5pNbases 6 --quantMode GeneCounts --outFilterMismatchNmax 10 --readFilesIn {input.expression_fastq} --outFileNamePrefix {params.accession}.
        '''

rule expression_bam_sort:
    input:
        expression_bam = data_dir+"/"+expression_dir+"/{accession}.Aligned.out.bam"

    params:
        accession = data_dir+"/"+expression_dir+"/{accession}"

    output:
        expression_bam = data_dir+"/"+expression_dir+"/{accession}.SortedAligned.out.bam"

    conda:
        "envs/samtools_env.yaml"

    shell:
        '''
        samtools sort -l 0 -@ 16 -T {params.accession} -O BAM -l 9 -o {output.expression_bam} {input.expression_bam}
        '''

rule expression_bam_duplicates_remove:
    input:
        expression_bam = data_dir+"/"+expression_dir+"/"+"{accession}.SortedAligned.out.bam",

    output:
        expression_bam = data_dir+"/"+expression_dir+"/{accession}.RmdupSortedAligned.out.bam"

    shell:
        '''
        python src/removeDuplicatesBam.py -s samtools -i {input.expression_bam} -o {output.expression_bam}
        '''

rule expression_bam_index:
    input:
        expression_bam = data_dir+"/"+expression_dir+"/{accession}.RmdupSortedAligned.out.bam",

    output:
        expression_bai = data_dir+"/"+expression_dir+"/{accession}.RmdupSortedAligned.out.bam.bai"

    conda:
        "envs/samtools_env.yaml"

    shell:
        '''
        samtools index {input.expression_bam}
        '''

rule pileup_generate:
    input:
        genome_fa = data_dir+"/"+genome_dir+"/"+genome_fa,
        expression_bam = data_dir+"/"+expression_dir+"/{accession}.RmdupSortedAligned.out.bam",
        expression_bai = data_dir+"/"+expression_dir+"/{accession}.RmdupSortedAligned.out.bam.bai"

    params:
        pileup_dir = data_dir+"/"+pileup_dir,
        chr = lambda wildcards, output: output[0].split(".")[1].split("_")[0],
        start = lambda wildcards, output: output[0].split(".")[1].split("_")[1],
        end = lambda wildcards, output: output[0].split(".")[1].split("_")[2]

    output:
        pileup_tsv = data_dir+"/"+pileup_dir+"/{accession}.{chr}.pileup.tsv"

    conda:
        "envs/pysam_env.yaml"

    shell:
        '''
        #mkdir -p {params.pileup_dir}
        python src/SNVinBam_pileup.py {input.genome_fa} {input.expression_bam} {params.chr} {params.start} {params.end} "30" "N" "N" |
        awk -v FS=\"\\t\" -v OFS=\"\\t\" '{{$2=$2 \"\t\" $2; print $0}}' |
        sort -k1,1 -k2,2n > {output.snv_map}
        '''

rule pileup_merge:
    input:
        pileup_tsv = expand(data_dir+"/"+pileup_dir+"/{{accession}}.{chr}.pileup.tsv", chr=[ "_".join([ chr_params[i][j] for j in chr_params[i] ]) for i in chr_params ])

    output:
        merged_pileup_tsv = data_dir+"/"+pileup_dir+"/{accession}.all.pileup.tsv"

    shell:
        '''
        cat {input.pileup_tsv} > {output.merged_pileup_tsv}
        '''

rule pileup_vcf_generate:
    input:
        merged_pileup_tsv = data_dir+"/"+pileup_dir+"/{accession}.all.pileup.tsv",
        genome_fa = data_dir+"/"+genome_dir+"/"+genome_fa,
        expression_bam = data_dir+"/"+expression_dir+"/{accession}.RmdupSortedAligned.out.bam",
        expression_bai = data_dir+"/"+expression_dir+"/{accession}.RmdupSortedAligned.out.bam.bai"

    output:
        merged_pileup_vcf = data_dir+"/"+pileup_dir+"/{accession}.vcf"

    conda:
        "envs/bcftools_env.yaml"

    shell:
        '''
        bcftools mpileup -B --skip-indels --min-BQ 30 --max-depth 50000 --ignore-RG -T {input.merged_pileup_tsv} --fasta-ref {input.genome_fa} {input.expression_bam} > {output.merged_pileup_vcf}
        '''

rule pileup_stats_generate:
    input:
        merged_pileup_vcf = data_dir+"/"+pileup_dir+"/{accession}.vcf"

    output:
        merged_pileup_stats_tsv = data_dir+"/"+pileup_dir+"/{accession}.stats.tsv"

    conda:
        "envs/bcftools_env.yaml"

    shell:
        '''
        bcftools query -f '%CHROM\t%POS\t%POS\t[%INFO/VDB]\t[%INFO/RPB]\t[%INFO/MQB]\t[%INFO/BQB]\t[%INFO/MQSB]\n' {input.merged_pileup_vcf} | sort -k1,1 -k2,2n > {output.merged_pileup_stats}
        '''

rule pileup_augment:
    input:
        merged_pileup_tsv = data_dir+"/"+pileup_dir+"/{accession}.all.pileup.tsv",
        merged_pileup_stats_tsv = data_dir+"/"+pileup_dir+"/{accession}.stats.tsv"

    output:
        augmented_merged_pileup_tsv = data_dir+"/"+pileup_dir+"/{accession}.all.augmented_pileup.tsv"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        bedtools intersect -f 1 -wa -wb -a {input.merged_pileup_tsv} -b {input.merged_pileup_stats_tsv} | cut -f1-2,4-11,15- > {output.augmented_merged_pileup_tsv}
        '''

rule mutation_map_generate:
    input:
        augmented_merged_pileup_tsv = data_dir+"/"+pileup_dir+"/{accession}.all.augmented_pileup.tsv"

    params:
        mutations_map_dir = data_dir+"/"+mutations_map_dir

    output:
        mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.unfiltered_map.tsv"

    shell:
        '''
        #mkdir -p {params.mutations_map_dir}
        python src/printMutationMap_fromPileup.py {input.augmented_merged_pileup_tsv} TRUE 6 40 0 0.7 {output.mutations_map_tsv}
        '''

rule germline_bed_generate:
    input:
        germline_variants_vcf = data_dir+"/"+misc_dir+"/"+germline_variants_vcf

    output:
        germline_variants_bed = data_dir+"/"+misc_dir+"/arabidopsis_2029_NoFilters.bed"

    shell:
        '''
        awk -v OFS="\t" '/^[^#]/ {{ print $1,($2-1),$2,$1"."$2,".","." }}' {input.germline_variants_vcf} > {output.germline_variants_bed}
        '''

rule germline_remove:
    input:
        germline_variants_bed = data_dir+"/"+misc_dir+"/arabidopsis_2029_NoFilters.bed",
        mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.unfiltered_map.tsv"

    output:
        germline_removed_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.germline_removed.tsv"

    shell:
        '''
        awk -v OFS="\t" '{{
            if(FNR==NR){{ a[$1,$3]; next }};
            if(($1,$2) in a){{
                if($8 == ""){{
                    print $0,"germline"
                }}else if($8 == "PASS"){{
                    print $1,$2,$3,$4,$5,$6,$7,"germline_variant"
                }}else{{
                    print $1,$2,$3,$4,$5,$6,$7,$8";germline_variant"
                }}
            }}else{{
                if($8 == ""){{
                    print $0,"PASS"
                }}else{{
                    print $0
                }}
            }}
        }}' {input.germline_variants_bed} {input.mutations_map_tsv} > {output.germline_removed_tsv}
        '''

rule vaf_filter:
    input:
        germline_removed_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.germline_removed.tsv"

    output:
        vaf_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.vaf.tsv"

    conda:
        "envs/scipy_env.yaml"

    shell:
        '''
        python src/vaf_remove.py {input.germline_removed_tsv} 0.5 0.05 "vaf" > {output.vaf_filter_tsv}
        '''

rule polynucleotide_bed_generate:
    input:
        genome_fa = data_dir+"/"+genome_dir+"/"+genome_fa

    output:
        polynucleotide_bed = data_dir+"/"+misc_dir+"/polynucleotide.bed"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/polynucleotide_filter.R {input.genome_fa} {output.polynucleotide_bed}
        '''

rule polynucleotide_filter:
    input:
        vaf_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.vaf.tsv",
        polynucleotide_bed = data_dir+"/"+misc_dir+"/polynucleotide.bed"

    output:
        polynucleotide_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.polynucleotide.tsv"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        python src/filterAppend_distanceToBed.py {input.vaf_filter_tsv} {input.polynucleotide_bed} 1 "polynucleotide_region" > {output.polynucleotide_filter_tsv}
        '''

rule blacklist_filter:
    input:
        polynucleotide_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.polynucleotide.tsv",
        blacklist_bed = data_dir+"/"+misc_dir+"/"+blacklist_bed

    output:
        blacklist_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.blacklist.tsv"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        python src/filterAppend_distanceToBed.py {input.polynucleotide_filter_tsv} {input.blacklist_bed} 1 "blacklisted_region" > {output.blacklist_filter_tsv}
        '''

rule rna_edit_bed_generate:
    input:
        tair9_genes_gff = data_dir+"/"+misc_dir+"/"+tair9_genes_gff,
        rna_edit_xls = data_dir+"/"+misc_dir+"/"+rna_edit_xls

    output:
        rna_edit_bed = data_dir+"/"+misc_dir+"/rna_edit.bed"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/rna_edit_filter.R {input.tair9_genes_gff} {input.rna_edit_xls} {output.rna_edit_bed}
        '''

rule rna_edit_filter:
    input:
        blacklist_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.blacklist.tsv",
        rna_edit_bed = data_dir+"/"+misc_dir+"/rna_edit.bed"

    output:
        rna_edit_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.rna_edit.tsv"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        python src/filterAppend_distanceToBed.py {input.blacklist_filter_tsv} {input.rna_edit_bed} 1 "rna_edit_site" > {output.rna_edit_filter_tsv}
        '''

rule exon_boundary_bed_generate:
    input:
        tair10_genes_gff = data_dir+"/"+misc_dir+"/"+tair10_genes_gff

    output:
        exon_boundary_bed = data_dir+"/"+misc_dir+"/exon_boundary.bed"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/exon_boundary_filter.R {input.tair10_genes_gff} {output.exon_boundary_bed}
        '''

rule exon_boundary_filter:
    input:
        rna_edit_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.rna_edit.tsv",
        exon_boundary_bed = data_dir+"/"+misc_dir+"/exon_boundary.bed"

    output:
        exon_boundary_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.exon_boundary.tsv"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        python src/filterAppend_distanceToBed.py {input.rna_edit_filter_tsv} {input.exon_boundary_bed} 7 "exon_boundary" > {output.exon_boundary_filter_tsv}
        '''

rule STAR_splice_junction_bed_generate:
    input:
        STAR_splice_junction_tab = data_dir+"/"+expression_dir+"/{accession}.SJ.out.tab"

    output:
        STAR_splice_junction_bed = data_dir+"/"+expression_dir+"/{accession}.SJ.bed"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/STAR_splice_junction_filter.R {input.STAR_splice_junction_tab} {output.STAR_splice_junction_bed}
        '''

rule STAR_splice_junction_filter:
    input:
        exon_boundary_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.exon_boundary.tsv",
        STAR_splice_junction_bed = data_dir+"/"+expression_dir+"/{accession}.SJ.bed"

    output:
        STAR_splice_junction_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.SJ.tsv"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        python src/filterAppend_distanceToBed.py {input.exon_boundary_filter_tsv} {input.STAR_splice_junction_bed} 7 "splice_junction_STAR" > {output.STAR_splice_junction_filter_tsv}
        '''

rule sequence_error_filter:
    input:
        STAR_splice_junction_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.SJ.tsv"

    output:
        sequence_error_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.sequence_error.tsv"

    conda:
        "envs/scipy_env.yaml"

    shell:
        '''
        python src/filterAppend_sequenceError.py {input.STAR_splice_junction_filter_tsv} "0.0001" 30 "sequencing_error" > {output.sequence_error_filter_tsv}
        '''

rule clustered_mutation_filter:
    input:
        sequence_error_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.sequence_error.tsv"

    output:
        clustered_mutation_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.clustered_mutation.tsv"

    conda:
        "envs/scipy_env.yaml"

    shell:
        '''
        python src/filterAppend_clusterMutations.py {input.sequence_error_filter_tsv} 100 3 "clustered_mutation" > {output.clustered_mutation_filter_tsv}
        '''

rule bcf_filter:
    input:
        augmented_merged_pileup_tsv = data_dir+"/"+pileup_dir+"/{accession}.all.augmented_pileup.tsv",
        clustered_mutation_filter_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.clustered_mutation.tsv"

    output:
        filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/{accession}.filtered_map.tsv"

    conda:
        "envs/scipy_env.yaml"

    shell:
        '''
        nCols=`awk '{{ print NF; exit }}' {input.augmented_merged_pileup_tsv}`
        if [ "$nCols" != "" ]
        then
            if [ "$nCols" -gt "10" ]
            then
                python src/filterAppend_statFromBCF.py {input.clustered_mutation_filter_tsv} {input.augmented_merged_pileup_tsv} "bcf" 0.05 0.05 0.05 0.05 0.05 > {output.filtered_mutations_map_tsv}
            else
                cp {input.clustered_mutation_filter_tsv} {output.filtered_mutations_map_tsv}
            fi
        else
            cp {input.clustered_mutation_filter_tsv} {output.filtered_mutations_map_tsv}
        fi
        '''

rule accession_group:
    input:
        filtered_mutations_map_tsv = expand(data_dir+"/"+mutations_map_dir+"/{accession}.filtered_map.tsv", accession=accessions)

    output:
        merged_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/all.filtered_map.tsv"

    shell:
        '''
        awk -v FS="\\t" -v OFS="\\t" '{{ n_a=split(FILENAME, a, "/"); split(a[n_a], b, "."); print b[1],$0 }}' {input.filtered_mutations_map_tsv} > {output.merged_filtered_mutations_map_tsv}
        '''

rule flagged_mutations_remove:
    input:
        merged_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/all.filtered_map.tsv"

    output:
        not_flagged_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/not_flagged.filtered_map.tsv"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/passed_mutations.R {input.merged_filtered_mutations_map_tsv} {output.not_flagged_filtered_mutations_map_tsv}
        '''

rule shared_mutations_remove:
    input:
        not_flagged_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/not_flagged.filtered_map.tsv"

    output:
        unique_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/unique.not_flagged.filtered_map.tsv"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/unique_mutations.R {input.not_flagged_filtered_mutations_map_tsv} {output.unique_filtered_mutations_map_tsv}
        '''

rule outlier_derive:
    input:
        unique_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/unique.not_flagged.filtered_map.tsv"

    output:
        outliers_tsv = data_dir+"/"+misc_dir+"/outliers.tsv"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/derived_outliers.R {input.unique_filtered_mutations_map_tsv} {output.outliers_tsv}
        '''

rule outlier_remove:
    input:
        unique_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/unique.not_flagged.filtered_map.tsv",
        outliers_tsv = data_dir+"/"+misc_dir+"/outliers.tsv"

    output:
        outlier_removed_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/outlier_removed.unique.not_flagged.filtered_map.tsv"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/outlier_mutations_removed.R {input.unique_filtered_mutations_map_tsv} {input.outliers_tsv} {output.outlier_removed_filtered_mutations_map_tsv}
        '''

rule mutations_process:
    input:
        outlier_removed_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/outlier_removed.unique.not_flagged.filtered_map.tsv"

    output:
        processed_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/processed.filtered_map.tsv"

    shell:
        '''
        awk -v OFS="\t" 'BEGIN{{ print "CHROM","POS","TYPE","REF","ALT","src" }} {{ if(FNR != 1){{ print $2,$3,"SNP",$4,$5,"epigenomes" }} }}' {input.outlier_removed_filtered_mutations_map_tsv} > {output.processed_mutations_map_tsv}
        '''

rule weigel_polynucleotide_filter:
    input:
        weigel_mutations_map_csv = data_dir+"/"+mutations_map_dir+"/"+weigel_raw_variants_csv,
        polynucleotide_bed = data_dir+"/"+misc_dir+"/polynucleotide.bed"

    output:
        weigel_polynucleotide_filter_tsv = data_dir+"/"+mutations_map_dir+"/weigel.polynucleotide.tsv"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/polynucleotide_weigel_remove.R {input.weigel_mutations_map_csv} {input.polynucleotide_bed} {output.weigel_polynucleotide_filter_tsv}
        '''

rule weigel_mutations_process:
    input:
        weigel_polynucleotide_filter_tsv = data_dir+"/"+mutations_map_dir+"/weigel.polynucleotide.tsv"

    output:
        processed_weigel_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/weigel.processed.filtered_map.tsv"

    conda:
        'envs/r_env.yaml'

    shell:
        '''
        Rscript src/process_weigel_mutations.R {input.weigel_polynucleotide_filter_tsv} {output.processed_weigel_mutations_map_tsv}
        '''

rule mutations_augment:
    input:
        tair10_union_exons_gff = data_dir+"/"+misc_dir+"/TAIR10_GFF3_exons.union.noambiguous.gff",
        processed_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/processed.filtered_map.tsv",
        processed_weigel_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/weigel.processed.filtered_map.tsv"

    output:
        augmented_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/augmented.filtered_map.tsv"

    conda:
        'envs/r_env.yaml'

    shell:
        '''
        Rscript src/augment_mutations.R {input.tair10_union_exons_gff} {input.processed_mutations_map_tsv} {input.processed_weigel_mutations_map_tsv} {output.augmented_mutations_map_tsv}
        '''

rule multi_instance_bed_generate:
    input:
        merged_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/all.filtered_map.tsv"

    output:
        multi_instance_bed = data_dir+"/"+misc_dir+"/multi_instance.bed"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/multi_instance_filter.R {input.merged_filtered_mutations_map_tsv} {output.multi_instance_bed}
        '''

rule vaf_bed_generate:
    input:
        merged_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/all.filtered_map.tsv"

    params:
        accession = "{accession}"

    output:
        vaf_bed = data_dir+"/"+misc_dir+"/{accession}.vaf.bed"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/from_mutation_map_filter.R {input.merged_filtered_mutations_map_tsv} {params.accession} "vaf" {output.vaf_bed}
        '''

rule seq_error_bed_generate:
    input:
        merged_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/all.filtered_map.tsv"

    params:
        accession = "{accession}"

    output:
        seq_error_bed = data_dir+"/"+misc_dir+"/{accession}.seq_error.bed"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/from_mutation_map_filter.R {input.merged_filtered_mutations_map_tsv} {params.accession} "sequencing_error" {output.seq_error_bed}
        '''

rule bcf_bed_generate:
    input:
        merged_filtered_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/all.filtered_map.tsv"

    params:
        accession = "{accession}"

    output:
        bcf_bed = data_dir+"/"+misc_dir+"/{accession}.bcf.bed"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/from_mutation_map_filter.R {input.merged_filtered_mutations_map_tsv} {params.accession} "bcf_variant_distance_bias;bcf_read_position_bias;bcf_mapping_quality_bias;bcf_base_quality_bias;bcf_mapping_quality_vs_strand_bias" {output.bcf_bed}
        '''

rule exon_boundary_range_bed_generate:
    input:
        exon_boundary_bed = data_dir+"/"+misc_dir+"/exon_boundary.bed"

    output:
        exon_boundary_range_bed = data_dir+"/"+misc_dir+"/exon_boundary.range.bed"

    shell:
        '''
        awk -v OFS="\t" '{{ l = ($2-6); if(l < 0){{ l = 0 }}; print $1,l,($3+6),$1"."$3,".","." }}' {input.exon_boundary_bed} > {output.exon_boundary_range_bed}
        '''

rule STAR_splice_junction_range_bed_generate:
    input:
        STAR_splice_junction_bed = data_dir+"/"+expression_dir+"/{accession}.SJ.bed"

    output:
        STAR_splice_junction_range_bed = data_dir+"/"+expression_dir+"/{accession}.SJ.range.bed"

    shell:
        '''
        awk -v OFS="\t" '{{ l = ($2-6); if(l < 0){{ l = 0 }}; print $1,l,($3+6),$1"."$3,".","." }}' {input.STAR_splice_junction_bed} > {output.STAR_splice_junction_range_bed}
        '''

rule coverage_map_generate: ##### need to test
    input:
        expression_bam = data_dir+"/"+expression_dir+"/{accession}.RmdupSortedAligned.out.bam"

    output:
        coverage_map_bg = data_dir+"/"+expression_dir+"/{accession}.coverage_map.bg.gz"

    conda:
        "envs/samtools_env.yaml"

    shell:
        '''
        samtools depth -q 29 -d 0 {input.expression_bam} | \
        awk -v OFS="\t" '{{ if($3 >= 40) print $1,$2,($2+1),"NA",$3,"*" }}' | \
        sort -k 1,1 -k2,2n | \
        gzip --best -c > {ouput.depth_bed}
        '''

rule coverage_map_filter:
    input:
        coverage_map_bg = data_dir+"/"+expression_dir+"/{accession}.coverage_map.bg.gz",
        germline_bed = data_dir+"/"+misc_dir+"/"+"arabidopsis_2029_NoFilters.bed",
        blacklist_bed = data_dir+"/"+misc_dir+"/"+"blacklist.bed",
        rna_edit_bed = data_dir+"/"+misc_dir+"/"+"rna_edit.bed",
        exon_boundary_range_bed = data_dir+"/"+misc_dir+"/exon_boundary.range.bed",
        polynucleotide_bed = data_dir+"/"+misc_dir+"/"+"polynucleotide.bed",
        multi_instance_bed = data_dir+"/"+misc_dir+"/"+"multi_instance.bed",
        vaf_bed = data_dir+"/"+misc_dir+"/"+"{accession}.vaf.bed",
        STAR_splice_junction_range_bed = data_dir+"/"+expression_dir+"/{accession}.SJ.range.bed",
        seq_error_bed = data_dir+"/"+misc_dir+"/"+"{accession}.seq_error.bed",
        bcf_bed = data_dir+"/"+misc_dir+"/"+"{accession}.bcf.bed"

    params:
        chr = "{chr}"

    output:
        coverage_map_filtered_bg = data_dir+"/"+expression_dir+"/{chr}.{accession}.coverage_map.filtered.bg.gz"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        gunzip -c {input.coverage_map_bg} |
        awk '$1 == {params.chr}' |
        bedtools subtract -a - -b {input.germline_bed} |
        bedtools subtract -a - -b {input.blacklist_bed} |
        bedtools subtract -a - -b {input.rna_edit_bed} |
        bedtools subtract -a - -b {input.exon_boundary_range_bed} |
        bedtools subtract -a - -b {input.polynucleotide_bed} |
        bedtools subtract -a - -b {input.multi_instance_bed} |
        bedtools subtract -a - -b {input.vaf_bed} |
        bedtools subtract -a - -b {input.STAR_splice_junction_range_bed} |
        bedtools subtract -a - -b {input.seq_error_bed} |
        bedtools subtract -a - -b {input.bcf_bed} |
        awk -v OFS="\t" '{{ print $1,$2,$3,$5 }}' | gzip > {output.coverage_map_filtered_bg}
        '''

rule bedg_union:
    input:
        coverage_map_filtered_bg = expand(data_dir+"/"+expression_dir+"/{{chr}}.{accession}.coverage_map.filtered.bg.gz", accession=accessions)

    output:
        union_coverage_map_filtered_bg = data_dir+"/"+expression_dir+"/{chr}.union.bg.gz"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        bedtools unionbedg -i {input.coverage_map_filtered_bg} | gzip > {output.union_coverage_map_filtered_bg}
        '''

rule tair10_exons_gff_generate:
    input:
        tair10_genes_gff = data_dir+"/"+misc_dir+"/"+tair10_genes_gff

    output:
        tair10_union_exons_gff = data_dir+"/"+misc_dir+"/TAIR10_GFF3_exons.union.noambiguous.gff"

    conda:
        "envs/bedtools_env.yaml"

    shell:
        '''
        awk -v OFS="\t" '{{
            if($3 == "exon"){{
                split($9, a, "Parent=");
                split(a[2], b, ".");
                print $1,$2,$3,$4,$5,$6,$7,$8,b[1]
            }}
        }}' {input.tair10_genes_gff} |
        sort -k1,1 -k4,4n -k5,5n |
        bedtools merge -i - -c 7,9 -o distinct |
        awk -v OFS="\t" '{{ print $1,"TAIR10","exon",($2+1),$3,".",$4,".",$5 }}' |
        awk '!index($9, ",")' > {output.tair10_union_exons_gff}
        '''

rule breadth_depth_process:
    input:
        tair10_union_exons_gff = data_dir+"/"+misc_dir+"/TAIR10_GFF3_exons.union.noambiguous.gff",
        genome_fa = data_dir+"/"+genome_dir+"/"+genome_fa,
        accessions_tsv = data_dir+"/"+misc_dir+"/"+accessions_tsv,
        outliers_tsv = data_dir+"/"+misc_dir+"/outliers.tsv",
        union_coverage_map_filtered_bg = data_dir+"/"+expression_dir+"/{chr}.union.bg.gz"

    output:
        processed_breadth_depth_tsv = data_dir+"/"+expression_dir+"/{chr}.breadth_depth.tsv"

    conda:
        "envs/r_env.yaml"

    shell:
        '''
        Rscript src/breadth_depth.R {input.tair10_union_exons_gff} {input.genome_fa} {input.accessions_tsv} {input.outliers_tsv} {input.union_coverage_map_filtered_bg} {output.processed_breadth_depth_tsv}
        '''

rule breadth_depth_merge:
    input:
        processed_breadth_depth_tsv = expand(data_dir+"/"+expression_dir+"/{chr}.breadth_depth.tsv", chr=chrs)

    output:
        merged_processed_breadth_depth_tsv = data_dir+"/"+expression_dir+"/all.breadth_depth.tsv"

    shell:
        '''
        cat {input.processed_breadth_depth_tsv} | awk -v OFS="\t" 'BEGIN{{ print "id","type","chr","left","right","ref","breadth","depth","strand","transcribed_status" }} {{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }}' > {output.merged_processed_breadth_depth_tsv}
        '''

rule GC_compute:
    input:
        tair10_union_exons_gff = data_dir+"/"+misc_dir+"/TAIR10_GFF3_exons.union.noambiguous.gff",
        genome_fa = data_dir+"/"+genome_dir+"/"+genome_fa

    output:
        gc_tsv = data_dir+"/"+genome_dir+"/GC.tsv"

    conda:
        'envs/r_env.yaml'

    shell:
        '''
        Rscript src/compute_GC.R {input.tair10_union_exons_gff} {input.genome_fa} {output.gc_tsv}
        '''

rule replication_timing_process:
    input:
        tair10_union_exons_gff = data_dir+"/"+misc_dir+"/TAIR10_GFF3_exons.union.noambiguous.gff",
        late_replication_bg = data_dir+"/"+replication_timing_dir+"/"+late_replication_bg,
        early_replication_bg = data_dir+"/"+replication_timing_dir+"/"+early_replication_bg

    output:
        replication_timing_ratio_tsv = data_dir+"/"+replication_timing_dir+"/replication_timing_ratio.tsv"

    conda:
        'envs/r_env.yaml'

    shell:
        '''
        Rscript src/process_replication_timing.R {input.tair10_union_exons_gff} {input.late_replication_bg} {input.early_replication_bg} {output.replication_timing_ratio_tsv}
        '''

rule histone_bigWig_process:
    input:
        histone_mark_replicate_bw = data_dir+"/"+histone_mark_dir+"/{histone_mark_replicate}.bw"

    output:
        histone_mark_replicate_bg = data_dir+"/"+histone_mark_dir+"/{histone_mark_replicate}.bg"

    conda:
        'envs/bw_to_bg.yaml'

    shell:
        '''
        bigWigToBedGraph {input.histone_mark_replicate_bw} {output.histone_mark_replicate_bg}
        '''

rule histone_signal_process:
    input:
        tair10_union_exons_gff = data_dir+"/"+misc_dir+"/TAIR10_GFF3_exons.union.noambiguous.gff",
        histone_mark_replicate_bg = expand(data_dir+"/"+histone_mark_dir+"/{histone_mark_replicate}.bg", histone_mark_replicate=histone_mark_replicates)

    params:
        histone_mark = "{histone_mark}"

    output:
        histone_mark_signal_tsv = data_dir+"/"+histone_mark_dir+"/{histone_mark}.signal.tsv"

    conda:
        'envs/r_env.yaml'

    shell:
        '''
        Rscript src/process_histone_mark.R {input.tair10_union_exons_gff} {input.histone_mark_replicate_bg} {params.histone_mark} {output.histone_mark_signal_tsv}
        '''

rule features_augment:
    input:
        tair10_union_exons_gff = data_dir+"/"+misc_dir+"/TAIR10_GFF3_exons.union.noambiguous.gff",
        genome_fa = data_dir+"/"+genome_dir+"/"+genome_fa,
        augmented_mutations_map_tsv = data_dir+"/"+mutations_map_dir+"/augmented.filtered_map.tsv",
        merged_processed_breadth_depth_tsv = data_dir+"/"+expression_dir+"/all.breadth_depth.tsv",
        gc_tsv = data_dir+"/"+genome_dir+"/GC.tsv",
        replication_timing_ratio_tsv = data_dir+"/"+replication_timing_dir+"/replication_timing_ratio.tsv",
        histone_mark_signal_tsv = expand(data_dir+"/"+histone_mark_dir+"/{histone_mark}.signal.tsv", histone_mark=histone_marks)

    params:
        results_dir = data_dir+"/"+results_dir

    output:
        augmented_features_tsv = data_dir+"/"+results_dir+"/augmented.features.tsv"

    conda:
        'envs/r_env.yaml'

    shell:
        '''
        mkdir -p {params.results_dir}
        Rscript src/augment_features.R {input.tair10_union_exons_gff} {input.genome_fa} {input.augmented_mutations_map_tsv} {input.merged_processed_breadth_depth_tsv} {input.gc_tsv} {input.replication_timing_ratio_tsv} {input.histone_mark_signal_tsv} {output.augmented_features_tsv}
        Rscript -e 'install.packages("gamlss.lasso", repos="https://cloud.r-project.org")'
        '''

rule model_fit:
    input:
        augmented_features_tsv = data_dir+"/"+results_dir+"/augmented.features.tsv",

    params:
        study = "{study}",
        ref = "{ref}",
        alt = "{alt}"

    output:
        model_fit_tsv = data_dir+"/"+results_dir+"/{study}.{ref}.{alt}.model_fit.tsv"

    conda:
        'envs/r_env.yaml'

    shell:
        '''
        Rscript src/fit_model_gamlsslasso.R {input.augmented_features_tsv} {params.study} {params.ref} {params.alt} {output.model_fit_tsv}
        '''

rule model_fit_merge:
    input:
        model_fit_tsv = expand(data_dir+"/"+results_dir+"/{study}.{ref}.{alt}.model_fit.tsv", zip, study=studies, ref=refs, alt=alts)

    output:
        merged_model_fit_tsv = data_dir+"/"+results_dir+"/all.model_fit.tsv"

    shell:
        '''
        cat {input.model_fit_tsv} | awk -v OFS="\t" 'BEGIN{{ print "study","ref","alt","covariate","estimate" }} {{ print $1,$2,$3,$4,$5 }}' > {output.merged_model_fit_tsv}
        '''
