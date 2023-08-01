/*
* Collect the reads
*/

process COLLECT_READS {
cpus 1
executor = 'local'
tag "$replicateId"
input:
	tuple val(replicateId), path(reads) 
	val SEQ
	val minbqual
	val r
	val minphred20
output:
	tuple val(replicateId), path('*150bp_R1.fastq.gz'), path('*150bp_R2.fastq.gz')
script:
"""
for i in *_R2*.fastq.gz; do basename=`ls \$i | cut -d "_" -f 1`; mv \$i \${basename}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R2.fastq.gz; done
for i in *_R1*.fastq.gz; do basename=`ls \$i | cut -d "_" -f 1`; mv \$i \${basename}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R1.fastq.gz; done
"""

}



/*
* Bam
*/

process MAPPING {
//errorStrategy 'ignore'
cpus 8
memory "20GB"
//conda "--copy bioconda::mtbseq"
tag "$replicateId"
publishDir "Bam", mode:'copy', pattern: "*bam*"
input:
	tuple val(replicateId), path(reads1), path(reads2)
output:
	tuple val(replicateId), path("Bam"), emit: BAM
	tuple val(replicateId), path("*bam*"), emit: bam
	val 'done', emit:done
script:
"""
USER=a perl /opt/conda/bin/MTBseq --step TBbwa --threads 8 || echo "processed \$?"
ln -s Bam/*bam* .
#touch ${replicateId}_a.bam
#touch ${replicateId}_a.bam.bai
#mkdir Bam

"""
}





process REFINE {
cpus 8
memory "20GB"
//conda "--copy bioconda::mtbseq"
tag "$replicateId"
publishDir "GATK_Bam", mode:'copy', pattern: "*gatk*"
input:
        tuple val(replicateId), path(bam1)
output:
        tuple val(replicateId), path("*gatk*"), emit: gatk
        tuple val(replicateId), path("GATK_Bam"), emit: GATK_Bam
	val 'done', emit:done
script:
"""
mkdir GATK_Bam
mkdir Bam
mv *bam* Bam/
USER=a perl /opt/conda/bin/MTBseq --step TBrefine --threads 8 || echo "processed \$?"
ln -s GATK_Bam/* .
#touch ${replicateId}_a.gatk.bam ${replicateId}_a.gatk.bai

"""
}




process PILE {
memory "20GB"
cpus 8
//conda "--copy bioconda::mtbseq"
tag "$replicateId"
publishDir "Mpileup", mode:'copy', pattern: "*mpileup*"
input:
        tuple val(replicateId), path(bam)
output:
        tuple val(replicateId), path("*mpileup*"), emit: mpile
        tuple val(replicateId), path("Mpileup"), emit: MPILE
	val 'done', emit:done
script:
"""
mkdir Mpileup
mkdir GATK_Bam
mv *gatk* GATK_Bam
USER=a perl /opt/conda/bin/MTBseq --step TBpile --threads 8 || echo "processed \$?"
ln -s Mpileup/* .
"""
}


process LIST {
memory "20GB"
cpus 8
//conda "--copy bioconda::mtbseq"
tag "$replicateId"
//publishDir "Position_Tables", mode:'copy', pattern: "*position_table*"
input:
        tuple val(replicateId), path(pile)
		val(minbq)
output:
        tuple val(replicateId), path("*position_table*"), emit: list
        tuple val(replicateId), path("Position_Tables"), emit: LIST
script:
"""
mkdir Mpileup
mv *mpileup* Mpileup/
mkdir Position_Tables
USER=a perl /opt/conda/bin/MTBseq --step TBlist --threads 8 --minbqual $minbq || echo "processed \$?"
ln -s Position_Tables/* .
"""
}

process VARIANTS_LOW {
memory "20GB"
cpus 1
//conda "--copy bioconda::mtbseq"
tag "$replicateId"
publishDir "Called", mode:'copy', pattern: "*tab"
input:
        tuple val(replicateId), path(list)
output:
        tuple val(replicateId), path("*tab"), emit: var_low
        tuple val(replicateId), path("Called"), emit: VAR_LOW
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Called
USER=a perl /opt/conda/bin/MTBseq --step TBvariants  --mincovf 1 --mincovr 1 --lowfreq_vars --minfreq 5 --minphred20 1 || echo "processed \$?"
ln -s Called/* .
"""
}


process VARIANTS {
memory "20GB"
cpus 1
//conda "--copy bioconda::mtbseq"
tag "$replicateId"
publishDir "Called", mode:'copy', pattern: "*tab"
input:
        tuple val(replicateId), path(list)
	val(mincovf)
	val(mincovr)
	val(minphred)
output:
        tuple val(replicateId), path("*tab"), emit: var
        tuple val(replicateId), path("Called"), emit: VAR
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Called
USER=a perl /opt/conda/bin/MTBseq --step TBvariants  --mincovf $mincovf --mincovr $mincovr --minphred20 $minphred || echo "processed \$?"
ln -s Called/* .
"""
}



process STATS {
memory "20GB"
cpus 1
//conda "--copy bioconda::mtbseq"
tag "$replicateId"
//publishDir "Statistics", mode:'copy', pattern: "*tab"
input:
        tuple val(replicateId), path(bam), path(list)
		val(mincovf)
        val(mincovr)
        val(minphred)
output:
        tuple val(replicateId), path("*tab"), emit: stats
        tuple val(replicateId), path("Statistics"), emit: STATS
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Bam
mv *bam* Bam/
mkdir Statistics
USER=a perl /opt/conda/bin/MTBseq --step TBstats  --mincovf $mincovf --mincovr $mincovr --minphred20 $minphred || echo "processed \$?"
mv Statistics/Mapping_and_Variant_Statistics.tab Statistics/${replicateId}_Mapping_and_Variant_Statistics.tab
ln -s Statistics/* .
"""
}


process STRAIN {
memory "20GB"
cpus 1
tag "$replicateId"
//publishDir "Classification", mode:'copy', pattern: "*tab"
input:
        tuple val(replicateId), path(list)
output:
        tuple val(replicateId), path("*.tab"), emit: strain
		tuple val(replicateId), path("Classification"), emit: STRAIN
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Classification
USER=a perl /opt/conda/bin/MTBseq --step TBstrains || echo "processed \$?"
mv Classification/Strain_Classification.tab Classification/${replicateId}_Strain_Classification.tab
ln -s Classification/* .
"""
}




process MAP_STRAIN{
memory '5GB'
cpus 1
publishDir "OUTPUT", mode: 'copy', pattern: "Mapping_Classification.tab"
input:
	path(stats)
output:
	path("Mapping_Classification.tab")
script:
"""
for i in \$(ls -1 *Mapping_and_Variant_Statistics* | cut -f1 -d '_' | sort -u)
do
paste \${i}_Mapp* \${i}_Strain* >> Mapping_Classification.tab
done

sort -u -r -o Mapping_Classification.tab Mapping_Classification.tab
"""
}


process DEL {
tag "$replicateId"
cpus 1
memory '20GB'
input:
	tuple val(replicateId), path(bam)
	path(ref)
	path(bed)
	path(bedix)
output:
	tuple val(replicateId), path("*.dels")
script:
"""
delly call -o ${replicateId}.bcf -g ${ref} ${replicateId}*.bam

bcftools annotate -a h37rv_ups_ordered.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') ${replicateId}.bcf | grep -P "<DEL>|<INS>" | grep "PASS" | grep -i -v "lowQual" |cut -f2,4,5,8,10 | awk '{ if (\$4 ~ "GENE") {print \$1,gensub(";.*","","g",gensub(".*;END=","","g",\$4)),gensub(">|<","","g",\$3),\$2,gensub(".*:([^:]*):([^:]*\$)","\\\\1;\\\\2","g",gensub(":0:0\$","","g",\$5)),gensub("(^[PI]).*","\\\\1","g",\$4),gensub(".*GENE=","","g",\$4)} else {print \$1,gensub(";.*","","g",gensub(".*;END=","","g",\$4)),gensub(">|<","","g",\$3),\$2,gensub(".*:([^:]*):([^:]*\$)","\\\\1;\\\\2","g",gensub(":0:0\$","","g",\$5)),gensub("(^[PI]).*","\\\\1","g",\$4)}}' OFS=';' | awk -F ';' '{if (\$2-\$1 < 100000) {print \$1,\$2, \$2-\$1,\$3,\$4,\$5,\$6,\$6/(\$6+\$5),\$7,\$8}}' OFS=';' > ${replicateId}.dels
"""
}

process OUT_DEL {
cpus 1
memory '20GB'
publishDir 'OUTPUT', mode: "copy", pattern: "DELETIONS.tab"
input:
	path(del)
output:
	path("DELETIONS.tab")
script:
"""

awk '{print gensub("\\.dels","","g",FILENAME),\$0}' OFS='\\t' *.dels | tr ';' '\\t' >> DELETIONS.tab

"""
}



process DEPTH {
cpus 8
memory "40GB"
tag "$replicateId"
input:
	tuple val(replicateId), path(bam)
	path(tgene)
output:
	tuple val(replicateId), path("gbcov_*")
script:
"""
mosdepth -b target_genes.bed -n -T ${task.cpus} ${replicateId} ${replicateId}*.bam

cat <(zcat ${replicateId}.thresholds.bed.gz | tail -n +2 | awk '{print \$1,\$2,\$3,\$4,\$5/(\$3-\$2)}' OFS='\\t' | sed 's/GB/GB_perc/g') <(zcat ${replicateId}.regions.bed.gz) | cut -f4,5 | awk -F '\t' '{for (i=1; i<=NF; i++) {a[NR,i] = \$i}}; NF>p { p = NF }; END {for (j=1; j<=p; j++) {str=a[1,j]; for (i=2; i<=NR; i++){str=str";"a[i,j];};print str}}' > gbcov_${replicateId}

rm ${replicateId}*

"""

}



process OUT_DEPTH {
cpus 1
memory '20GB'
publishDir 'OUTPUT', mode: "copy", pattern: "GB_cov.csv"
input:
	path(cov)
output:
	path("GB_cov.csv")
script:
"""
for i in gbcov_*;do awk '{if (NR == 1) {print "000000",\$0} else {print gensub(".*_","","g",FILENAME),\$0} }' OFS=';' \${i}; done | sort -u -k1 -t ';' | sed 's/000000/Samples/g' >> GB_cov.csv

sort -u -o GB_cov.csv GB_cov.csv
sort -k2 -r -t ';' -o GB_cov.csv GB_cov.csv

"""
}


process MUT_CORRECTION {
cpus 1
memory "20GB"
tag "$replicateId"
publishDir "Called", mode:'copy', pattern: "*corrected.tab"
input:
	tuple val(replicateId), path(var)
output:
	tuple val(replicateId), path('*corrected.tab')

script:
"""
#!/usr/bin/env Rscript

library(tidyverse)
library(Biostrings)



l=list.files()



l<-l[str_detect(l,'variants_cf1_cr1_fr5_ph1_outmode001.tab\$')]

for (i in l){
 print(i)
  s=str_extract(i,'[^_]*_[^_]*')


  read_delim(i,delim='\\t',show_col_types = FALSE) %>%
        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character)) %>%
        filter(Subst!=' ',Type=='SNP') %>%
    extract(Subst,c('refA','pA','varA','codR','cV1','cV2','cV3'),'([A-z]*)([0-9]*)([A-z]*) \\\\((.*)/(.)(.)(.)\\\\)') %>%
    group_by(Gene,pA) %>%
    filter(n()>1,sd(Freq)<15) %>%
    mutate(across(starts_with('cV'),function(x){ifelse(sum(LETTERS%in%x)>0,x[x%in%LETTERS],x)})) %>%
    unite('codV',cV1:cV3,sep='') %>%
    {if(dim(.)[1]>0) mutate(.,varA=plyr::mapvalues(as.character(translate(DNAString(codV[1]))),from = names(AMINO_ACID_CODE),to=AMINO_ACID_CODE,warn_missing = F)) else
            mutate(.,varA='A')}%>%
    arrange(`#Pos`) %>%
    mutate(Ref=ifelse(n()==2 &max(diff(`#Pos`))==2,codR,(str_c(Ref,collapse=''))),
           Allel=ifelse(n()==2 &max(diff(`#Pos`))==2,codV,(str_c(Allel,collapse='')))#,
           #`#Pos`=str_c(`#Pos`,collapse=',')
    ) %>%
    mutate(across(c('Ref','Allel'),function(x){ifelse(str_detect(Gene,'c'),as.character(reverseComplement(DNAString(x[1]))),x[1])})) %>%
    unite('Subst',refA,pA,varA,sep='') %>% mutate(Subst=paste(Subst,' (',codR,'/',codV,')',sep='')) %>% select(-codR,-codV) %>%
    mutate(across(where(is.logical),function(x){x=as.character(x)})) %>%
    rbind(read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                       across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
        filter(Subst!=' ',Type=='SNP') %>%
            extract(Subst,c('refA','pA','varA','codR','cV1','cV2','cV3'),'([A-z]*)([0-9]*)([A-z]*) \\\\((.*)/(.)(.)(.)\\\\)') %>%
            group_by(Gene,pA) %>%
            filter(n()<2|sd(Freq)>=15) %>%  unite('codV',cV1:cV3,sep='') %>%
            unite('Subst',refA,pA,varA,sep='') %>% mutate(Subst=paste(Subst,' (',codR,'/',codV,')',sep='')) %>%
            #mutate(`#Pos`=as.character(`#Pos`)) %>%
            select(-codR,-codV),
          read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
            #mutate(`#Pos`=as.character(`#Pos`)) %>%
            filter(Subst==' ',Type=='SNP'),
          read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
        filter(Type!='SNP') %>% arrange(Type,`#Pos`) %>%
            group_by(Type) %>%
            mutate(p=cumsum(c(1,diff(`#Pos`)>1))) %>%
            group_by(Type,Subst,Gene,GeneName,Product,ResistanceSNP,PhyloSNP,InterestingRegion,p) %>%
            summarise(`#Pos`=min(`#Pos`),Insindex=0,Ref='_',Allel=str_to_upper(Type),across(CovFor:Cov,mean)) %>% arrange(Type,p) %>% ungroup() %>%
            select(-p)
    ) ->a
    #%>% write_delim(paste(i,'corrected.tab',sep='_'),delim='\\t')

  if (file.exists(paste('../delly/',s,'.dels',sep=''))) {
    a %>% bind_rows(read_delim(paste('../delly/',s,'.dels',sep=''),show_col_types = FALSE,delim='\\t',
                         col_names = c('Start','End','Length','Type','Ref','RefR','VarR','Freq','Precision','Gene')) %>%
                {if(dim(.)[1]>0) mutate(.,Insindex=0,Ref='_',Allel=Type,Type=str_to_title(Type),Subst=" ",GeneName='-', Product=" ",Freq=Freq*100,Qual20=RefR+VarR) %>%
                select(`#Pos`=Start,Insindex,Ref,Type,Allel,Subst,Gene,GeneName,Product,Freq,Qual20)})->a}

  a %>%
write_delim(paste(i,'corrected.tab',sep='_'),delim='\\t')
}

"""
}
