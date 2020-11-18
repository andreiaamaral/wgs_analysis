/*  GATK4 Variant Calling Pipeline 
 *  Usage: nextflow run gencorefacility/variant-calling-pipeline-gatk4 -with-docker gencorefacility/variant-calling-pipeline-gatk4 
 *
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu >
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line
params.reads="/home/ec2-user/environment/data/ggal/*_{1,2}.*fq"
params.ref="/home/ec2-user/environment/data/ggal/transcriptome.fa"
params.outdir="/home/ec2-user/environment/result_alignment"
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/gatk_temp"
//params.snpeff_data = "${params.outdir}/snpeff_data"

// Print some stuff here
println "reads: $params.reads"
println "ref: $params.ref"
println "output: $params.out"
//println "gatk temp dir: $params.tmpdir"
//println "snpeff db: $params.snpeff_db"
//println "snpeff data: $params.snpeff_data"

// Setup the reference channel
/*
Channel
    .fromPath(params.ref)
    .ifEmpty { error "Cannot find any reads matching: ${params.ref}"  }
    .tap { ref_ch }
*/


/* Prepare the fastq read pairs for input.
 * While doing this, count number of input samples
 */
num_samples = 0
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .tap { read_pairs_ch }
    .subscribe({ num_samples += 1 })
    
    //read_pairs_ch.view()
    
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


process markDuplicatesSpark {
    publishDir "${params.out}/dedup_sorted", mode:'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    input:
    set val(pair_id), path(bamfile), path(baifile) from bam_for_sorting_bam_ch

    output:
    set val(pair_id), val(1),  path("${pair_id}_sorted_dedup.bam") into bam_for_variant_calling_ch, sorted_dedup_ch_for_metrics, bam_for_bqsr
    set val(pair_id), path("${pair_id}_dedup_metrics.txt") into dedup_qc_ch

    script:
    """
    mkdir tmp
    gatk --java-options "-Djava.io.tmpdir=tmp"  MarkDuplicatesSpark -I ${bamfile} -M ${pair_id}_dedup_metrics.txt -O ${pair_id}_sorted_dedup.bam --remove-sequencing-duplicates
    
    """ 
}