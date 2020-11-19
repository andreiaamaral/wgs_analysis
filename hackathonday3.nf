/*  Population genomics pipeline
 *  Usage: 
 *
 *
 *  Author: participants of Hackaton 2
 /*
 log.info """\
         POP_GENOMICS - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()
   
 */

// Setting some defaults here,
// Establishing parameters
params.reads="/home/ec2-user/environment/data/ggal/*_{1,2}.*fq"
params.ref="/home/ec2-user/environment/data/ggal/transcriptome.fa"
params.outdir="/home/ec2-user/environment/result_alignment"
params.out = "${params.outdir}/out"
//params.tmpdir = "${params.outdir}/gatk_temp"
//params.gatk = "/opt/broad/GenomeAnalysisTK.jar" 
//params.snpeff_data = "${params.outdir}/snpeff_data"
num_samples = 0



Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .tap { read_pairs_ch }
    .subscribe({ num_samples += 1 })
    
    //read_pairs_ch.view()
    


process foo {
    echo true
    input:
    path transcriptome from params.ref
    
    output:
    file("*.dict") into picard_dict_ch
    file("*.fai") into samtools_index_ch
    
    script:
    """
    echo "java -jar /data/software/picard/build/libs/picard.jar CreateSequenceDictionary R= ${transcriptome}  O=${transcriptome}.dict" > ${transcriptome}.dict

    echo "samtools faidx ${transcriptome}" > ${transcriptome}.fai
    """
 }
 
 
 //picard_dict_ch.view()
 
 //samtools_index_ch.view()


process index {
    container 'gencorefacility/variant-calling-pipeline-gatk4'
    input:
    path transcriptome from params.ref
     
    output:
    //path 'index' into index_ch
    file("*") into index_ch

    script:       
    """
    bwa index ${transcriptome}
    """
}

//index_ch.view()


process align {
    cpus 1
  
    publishDir "${params.out}/aligned_reads", mode:'copy' // tells the output of the chanel it is created automatically by nextflow
	container 'gencorefacility/variant-calling-pipeline-gatk4'
    input:
    set pair_id, file(reads) from read_pairs_ch //the chanel creates the pair ids
    file(index) from index_ch
    path(ref) from params.ref
    
     
    output:
    set val(pair_id), file("${pair_id}_aligned_reads.sam") \
	into aligned_reads_ch
	
    script:
    """
    bwa mem ${ref} ${reads} > ${pair_id}_aligned_reads.sam
    """

}

//aligned_reads_ch.view()

process sam_to_sorted_index_bam{
    publishDir "${params.out}/sort_indexed_bam", mode:'copy'
    container 'gencorefacility/variant-calling-pipeline-gatk4'
    input:
    set val(pair_id), file(aligned_reads) from aligned_reads_ch
    
    output:
    set val(pair_id),  file("${pair_id}_aligned_reads_sorted.bam"), file("${pair_id}_aligned_reads_sorted.bam.bai")  into bam_for_sorting_bam_ch
    
    script:
    """
    samtools view -S -b ${aligned_reads} > ${pair_id}_aligned_reads.bam
    samtools sort ${pair_id}_aligned_reads.bam -o ${pair_id}_aligned_reads_sorted.bam
    samtools index ${pair_id}_aligned_reads_sorted.bam
    """
   
}

//bam_for_sorting_bam_ch.view()


process markDuplicates {
    echo true
    publishDir "${params.out}/dedup_sorted", mode:'copy'
    container 'cbcrg/callings-nf@sha256:b65a7d721b9dd2da07d6bdd7f868b04039860f14fa514add975c59e68614c310'
    input:
    set val(pair_id), path(bamfile), path(baifile) from bam_for_sorting_bam_ch

    //output:
    //set val(pair_id),  path("${pair_id}_sorted_duplicates_rm.bam") into bam_for_variant_calling_ch, sorted_dedup_ch_for_metrics, bam_for_bqsr
    //set val(pair_id), path("${pair_id}_sorted_duplicates_metrics.txt") into dedup_qc_ch

    script:
    """
  
   echo "java -jar /data/software/picard/build/libs/picard.jar MarkDuplicates I=${bamfile} O=${pair_id}_sorted_duplicates_rm.bam M=${pair_id}_sorted_duplicates_metrics.txt REMOVE_SEQUENCING_DUPLICATES=true"
    
    """ 
}

process variant_calling {
    echo true
    input:
    
    output:
    
    script:
    """

 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R GCF_002310715.1_ASM231071v1_genomic.fna \
   -I E2_alnRGsorted_duplicates_removed.bam \
   -O E2_output.vcf.gz \
   -bamout E2_bamout.bam
   
9. Variant filter
   
   gatk VariantFiltration \
   -R GCF_002310715.1_ASM231071v1_genomic.fna \
   -V E2_output.vcf.gz \
   -O E2_output_filtered.vcf.gz \
   --filter-name "my_filter1" \
   --filter-expression "QUAL < 0 || DP<10.0 || MQ < 30.00 || SOR > 10.000 || QD < 2.00 || QD> 5.00|| FS > 200.000 || ReadPosRankSum < -20.000 || ReadPosRankSum > 20.000"  
   
   
   
10. After having all individual VCFs merged them in a single VCF
   
   java -jar /data/software/picard/build/libs/picard.jar MergeVcfs \
          I=E1_output_filtered.vcf.gz\
          I=E2_output_filtered.vcf.gz \
          O=E12output_variants.vcf.gz