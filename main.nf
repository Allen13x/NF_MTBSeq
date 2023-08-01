nextflow.enable.dsl = 2
autoMounts = true
/*
 * Define the default parameters
 */ 
	params.reads	= "$baseDir/*_R{1,2}*.fastq.gz"
	params.results	= "OUTPUT"
	params.SEQ	= "ILL"
	params.minbqual	= "13"
	params.RP	= "0"
	params.minphred20	= "4"
	params.mincovf	= "4"
	params.mincovr	= "4"
	params.ref = "$baseDir/REF/M._tuberculosis_H37Rv_2015-11-13.fasta"
	params.bed = "/beegfs/datasets/buffer/ric.cirillo/MTB/h37rv_ups_ordered.bed.gz"
	params.bedix= "/beegfs/datasets/buffer/ric.cirillo/MTB/h37rv_ups_ordered.bed.gz.tbi"
	params.tgene="/idle/ric.cirillo/dimarco.federico/tbseq/Cov_utils/target_genes.bed"
log.info """\
C A L L I N G S  -  N F    v 2.1 
================================
reads   	: $params.reads
reference	: $params.ref
bed			: $params.bed
SEQ		: $params.SEQ
minbq		: $params.minbqual
REP REG		: $params.RP
minphred20	: $params.minphred20
mincovF		: $params.mincovf
mincovR		: $params.mincovr
results		: $params.results
"""

/* 
 * Import modules 
 */

include{COLLECT_READS;
	MAPPING;
	REFINE;
	PILE;
	LIST;
	VARIANTS_LOW;
	VARIANTS;
	STATS;
	STRAIN;
	MAP_STRAIN;
	DEL;
	OUT_DEL;
	DEPTH;
	OUT_DEPTH;
	MUT_CORRECTION} from "$baseDir/module.nf"
/* 
 * main pipeline logic
 */


workflow {
reads_ch=channel.fromFilePairs(params.reads)
COLLECT_READS(reads_ch,params.SEQ,params.minbqual,params.RP,params.minphred20)
MAPPING(COLLECT_READS.out)
REFINE(MAPPING.out.bam)
PILE(REFINE.out.gatk)
LIST(PILE.out.mpile,params.minbqual)
VARIANTS_LOW(LIST.out.list)
VARIANTS(LIST.out.list,params.mincovf,params.mincovr,params.minphred20)
DEL(MAPPING.out.bam,params.ref,params.bed,params.bedix)
DEPTH(MAPPING.out.bam,params.tgene)
MUT_CORRECTION(VARIANTS_LOW.out.var_low)
STATS(MAPPING.out.bam.join(LIST.out.list,by: 0),params.mincovf,params.mincovr,params.minphred20)
STRAIN(LIST.out.list)
map_strain=STATS.out.stats.join(STRAIN.out.strain,by:0).map{id,file1,file2 -> tuple(file1,file2)}.collect()
old_map=channel.fromPath('OUTPUT/Mapping_Classification.tab')
map_strain=map_strain.concat(old_map).collect()
MAP_STRAIN(map_strain)
delly=DEL.out.map{id,file -> file}
old_del=channel.fromPath('OUTPUT/DELETIONS.tab')
delly=delly.concat(old_del).collect()
OUT_DEL(delly)
depth=DEPTH.out.map{id,file->file}
old_cov=channel.fromPath('OUTPUT/GB_cov.csv')
depth=depth.concat(old_cov).collect()
OUT_DEPTH(depth)
}
