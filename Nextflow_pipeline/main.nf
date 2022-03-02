#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '1.0'

// this prints the input parameters
log.info """
AMRtcp-Pipeline ~  version ${version}
=============================================
reads                  : ${params.reads}
"""


Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: "${params.reads} }
    .set { reads_for_fastqc }

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: "${params.reads} }
    .set { reads_for_trimming }

Channel
    .value( params.illuminaclip)
    .set { illuminaclip }

Channel
    .value( params.db_16s)
    .set { db_16s }

Channel
    .value( params.db_rMLST)
    .set { db_rMLST }

Channel
    .value( params.bigsdb_rMLST)
    .set { bigsdb_rMLST }

Channel
    .value( params.abricate_dbs )
    .set { abricate_dbs }

Channel
    .value( params.ariba_db_ncbi )
    .set { ariba_db_ncbi }

Channel
    .value( params.rgi_db )
    .set { rgi_db }

Channel
    .value( params.deeparg_db )
    .set { deeparg_db }

Channel
    .value( params.resfinder_db )
    .set { resfinder_db }

Channel
    .value( params.resfinder_species)
    .set { resfinder_species }


/*
* including the modules and adding the parameter OUTPUT to pass them the folders where to publish the results
*/
include { fastqc } from "./modules/fastqc"
include { multiqc } from "./modules/multiqc"
include { trimmomaticPE } from "./modules/trimmomatic"
include { unicycler } from "./modules/unicycler"
include { bwaIndex } from "./modules/bwa_index"
include { bwaAlign } from "./modules/bwa-mem"
include { samtools } from "./modules/samtools"
include { pilon; pilon_remapping } from "./modules/pilon"
include { prokka } from "./modules/prokka"
include { quast } from "./modules/quast"
include { rMLST; call_rMLST } from "./modules/rMLST"
include { metaphlan3 } from "./modules/metaphlan"
include { make_one_contig; parse_sam_for_insertsize; coverage_pilon_corrected } from "./modules/python_functions"
include { bwaIndex as indexRemapping } from "./modules/bwa_index"
include { bwaAlign as alignRemapping } from "./modules/bwa-mem"
include { samtools as samtoolsRemapping} from "./modules/samtools"
include { typing_16S } from "./modules/typing_16S.nf"
include { abricate } from "./modules/abricate"
include { amrfinder_nuc; amrfinder_prot } from "./modules/AMRFinderPlus.nf"
include { ariba } from "./modules/ariba"
include { rgi } from "./modules/rgi"
include { sraX_basic; sraX_ext } from "./modules/sraX"
include { deeparg_LS; deeparg_SR } from "./modules/deeparg"
include { resfinder_fasta; resfinder_reads } from "./modules/resfinder"
include { summary_sample; merge_summaries } from "./modules/summary"


workflow {
  fastqc_out = fastqc(reads_for_fastqc)
  multiqc(fastqc_out.collect())
  trimm_out = trimmomaticPE(reads_for_trimming, illuminaclip)
  unicycler_out = unicycler(trimm_out.trimmed_reads)
  bwa_index_polishing = bwaIndex(unicycler_out.assembly)
  mapping = bwaAlign(trimm_out.trimmed_reads.join(bwa_index_polishing))
  bam = samtools(mapping)
  polished_assembly = pilon(bam.join(unicycler_out.assembly))
  annotation = prokka(polished_assembly.assembly)
  assembly_stats = quast(annotation.fna)
  typing_rMLST = rMLST(annotation.fna, db_rMLST)
  rmlst_out = call_rMLST(typing_rMLST, bigsdb_rMLST)
  metaphlan_out = metaphlan3(trimm_out.trimmed_reads)
  one_contig = make_one_contig(annotation.fna)
  bwa_index_remapping = indexRemapping(one_contig)
  remapping = alignRemapping(trimm_out.trimmed_reads.join(bwa_index_remapping))
  bam_remapping = samtoolsRemapping(remapping)
  insertsize = parse_sam_for_insertsize(remapping)
  remapping_polished = pilon_remapping(bam_remapping.join(one_contig))
  coverage = coverage_pilon_corrected(remapping_polished.vcf)
  typ16S = typing_16S(one_contig, db_16s)
  abricate(annotation.fna, abricate_dbs)
  amrfinder_nuc(annotation.fna)
  amrfinder_prot(annotation.faa)
  rgi(annotation.fna, rgi_db)
  sraX_basic(annotation.fna)
  sraX_ext(annotation.fna)
  deeparg_LS(annotation.fna, deeparg_db)
  deeparg_SR(trimm_out.trimmed_reads, deeparg_db)
  resfinder_fasta(annotation.fna, resfinder_db, resfinder_species)
  resfinder_reads(trimm_out.trimmed_reads, resfinder_db, resfinder_species)
  single_summary = summary_sample(trimm_out.trim_log.join(coverage).join(insertsize).join(assembly_stats.tsv).join(typ16S).join(metaphlan_out.profile).join(rmlst_out))
  summary = merge_summaries(single_summary.sample_quality.collect())
}
