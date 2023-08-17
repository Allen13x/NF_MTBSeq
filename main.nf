nextflow.enable.dsl = 2
autoMounts = true
/*
 * Define the default parameters
 */ 
	params.reads	= "$baseDir/samp"
	params.results	= "OUTPUT"
	params.SEQ	= "ILL"
	params.minbqual	= "13"
	params.RP	= "0"
	params.minphred20	= "4"
	params.mincovf	= "4"
	params.mincovr	= "4"
	params.ref="M._tuberculosis_H37Rv_2015-11-13"
	params.reff = "${baseDir}/REF/${params.ref}.fasta"
	params.bed = "$baseDir/REF/h37rv_ups_ordered.bed.gz"
	params.bedix= "$baseDir/REF/h37rv_ups_ordered.bed.gz.tbi"
	params.tgene="$baseDir/REF/target_genes.bed"
	params.pharma=false
	params.pgene="$baseDir/REF/gene_drug.csv"
	params.tdrug=10
	params.dhead="/$baseDir/REF/head"
	params.WHO="/$baseDir/REF/WHO_custom.csv"
	params.extra=true
	params.join=true
	params.sj='sample_joint'
	params.proj='def'
	params.ascii="$baseDir/REF/ascii_string"
	params.h=false


/* 
 * Import modules 
 */

include{COLLECT_READS;
	COLLECT_READS_ONT;
	MAPPING;
	MAPPING_ONT;
	REFINE;
	REFINE_ONT;
	PILE;
	PILE_ONT;
	LIST;
	VARIANTS_LOW;
	VARIANTS;
	STATS;
	STRAIN;
	MAP_STRAIN;
	JOIN;
	DEL;
	DEL_ONT;
	OUT_DEL;
	DEPTH;
	OUT_DEPTH;
	MUT_CORRECTION;
	MUT_GATHER;
	PHARMA;
	WHO;
	OUT_WHO} from "$baseDir/module.nf"
/* 
 * main pipeline logic
 */


workflow {
if (params.h){
log.info """
  --SEQ		Sequencing technology:
				  ILL: Illumina [default]
				  ONT: Oxford Nanopore technology  
  --ref        Reference Genome to use:
                  M._abscessus_CIP-104536T_2014-02-03
                  M._chimaera_DSM44623_2016-01-28
                  M._fortuitum_CT6_2016-01-08
                  M._tuberculosis_H37Rv_2015-11-13 [default]
  --join       Perform joint analysis [default: true]
  --sj         List of samples to use for joint analysis [optional]
  --proj       Name of the project for joint analysis [default: def]
  --pharma     Perform Drug resistance analysis at Custom frequencies (to use after first round of analysis) [default: false]
  --tdrug      Custom frequencies % [default: 10]
  --extra      Perform extra analysis [default: true]:
                 Deletion/Insertion detection with Delly2
                 Sequencing depth with Mosdepth
                 Pharma analysis at 10% and Comparison with WHO catalogue   
"""
}else{
if (params.pharma){

PHARMA(Channel.fromPath('Called/*corrected.tab').collect(),params.tdrug,params.pgene)


}
else{
log.info """\
================================
reads   	: $params.reads + '*_R{1,2}*.fastq.gz'
reference	: $params.ref
bed			: $params.bed
SEQ		: $params.SEQ
minbq		: $params.minbqual
REP REG		: $params.RP
minphred20	: $params.minphred20
mincovF		: $params.mincovf
mincovR		: $params.mincovr
results		: $params.results
Interesing genes: $params.bed
Target genes: $params.tgene
Drug genes: $params.pgene
Mutation Threshold: $params.tdrug
WHO Catalogue: $params.WHO
"""

if(params.SEQ == "ILL"){

reads_ch=channel.fromFilePairs(params.reads + '*_R{1,2}*.fastq.gz').map{id,file ->tuple((id - ~/_.*/),file)}
//reads_ch.view()
COLLECT_READS(reads_ch,params.SEQ,params.minbqual,params.RP,params.minphred20)
collected=COLLECT_READS.out
MAPPING(COLLECT_READS.out,params.ref)
mapped=MAPPING.out
REFINE(MAPPING.out.bam,params.ref)
refined=REFINE.out
PILE(REFINE.out.gatk,params.ref)
piled=PILE.out
}
else{
reads_ch=channel.fromPath(params.reads + '/*fastq.gz').map{file ->tuple((file.getSimpleName() - ~/_.*/),file)}
//reads_ch.view()
COLLECT_READS_ONT(reads_ch,params.SEQ,params.minbqual,params.RP,params.minphred20)
collected=COLLECT_READS_ONT.out
MAPPING_ONT(COLLECT_READS_ONT.out,params.ref)
mapped=MAPPING_ONT.out
REFINE_ONT(MAPPING_ONT.out.bam,params.ref,params.ascii)
refined=REFINE_ONT.out
PILE_ONT(REFINE_ONT.out.gatk,params.ref,params.minbqual)
piled=PILE_ONT.out}
LIST(piled.mpile,params.minbqual,params.ref)
VARIANTS_LOW(LIST.out.list,params.ref)
VARIANTS(LIST.out.list,params.mincovf,params.mincovr,params.minphred20,params.ref)
STATS(mapped.bam.join(LIST.out.list,by: 0),params.mincovf,params.mincovr,params.minphred20)
STRAIN(LIST.out.list)
map_strain=STATS.out.stats.join(STRAIN.out.strain,by:0).map{id,file1,file2 -> tuple(file1,file2)}.collect()
old_map=channel.fromPath('OUTPUT/Mapping_Classification_clean.tab')
map_strain=map_strain.concat(old_map).collect()


if (params.extra){
if (params.SEQ == "ILL"){	
DEL(mapped.bam,params.ref,params.bed,params.bedix)
deletion=DEL.out}
else{
DEL_ONT(mapped.bam,params.ref,params.bed,params.bedix)
deletion=DEL_ONT.out}
DEPTH(mapped.bam,params.tgene)
MUT_CORRECTION(VARIANTS_LOW.out.var_low)
delly=deletion.map{id,file -> file}
old_del=channel.fromPath('OUTPUT/DELETIONS.*')
delly=delly.concat(old_del).collect()
OUT_DEL(delly)
depth=DEPTH.out.map{id,file->file}
old_cov=channel.fromPath('OUTPUT/GB_cov.*')
depth=depth.concat(old_cov).collect()
OUT_DEPTH(depth)
mut=MUT_CORRECTION.out
old_mut=Channel.fromPath('Called/*corrected.tab').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}
mut=mut.concat(old_mut).unique{it[0]}.map{id,file->file}.collect()
//mut.view()
MUT_GATHER(mut)
PHARMA(mut,"10",params.pgene)
MAP_STRAIN(map_strain)
WHO(MUT_GATHER.out,params.dhead,params.WHO)
OUT_WHO(WHO.out)
}

if (params.join){
call=VARIANTS.out.var
old_call=Channel.fromPath('Called/*variants_cf4*').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}
call=call.concat(old_call).unique{it[0]}.map{id,file->file}.collect()
list=LIST.out.list
old_list=Channel.fromPath('Position_Tables/*').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}
list=list.concat(old_list).unique{it[0]}.map{id,file->file}.collect()
JOIN(call,list,params.sj,params.minbqual,params.minphred20,params.proj,params.ref)

}


}
}}