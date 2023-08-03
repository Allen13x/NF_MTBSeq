/*
* Collect the reads
*/

process COLLECT_READS {
cpus 1
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
for i in *_R2*.fastq.gz; do basename=`ls \$i | cut -d "_" -f 1`;mv \$i \${basename}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R2.fastq.gz; done
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
publishDir "Position_Tables", mode:'copy', pattern: "*position_table*"
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



process JOIN {
memory "60GB"
cpus 1
publishDir "Joint", mode: "copy", pattern: "Joint/*"
publishDir "Amend", mode: "copy", pattern: "Amend/*"
publishDir "Groups", mode: "copy", pattern: "Groups/*"
input:
        path(cc)
        path(ll)
        val(sample_joint)
	val(minbqual)
	val(minphred20)
        val(proj)
output:
        path("Joint/*"), emit: joint
        path("Amend/*"), emit: amend
        path("Groups/*"), emit: groups
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Called
mv *cf4* Called/
mkdir Joint
mkdir Amend
mkdir Groups
ls -1 Called/*_variants_cf4* | cut -f2 -d'/' | cut -f1,2 -d '_' | tr '_' '\\t' | sort -r | sort -u -k1,1 > sample_joint
USER=a perl /opt/conda/bin/MTBseq --step TBjoin --continue --samples ${sample_joint} --distance 5 --project ${proj} --minbqual ${minbqual} --minphred20 ${minphred20} || echo "processed \$?"
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

sort -u -r -k2 Mapping_Classification.tab | cut --complement -f1,25 |uniq >Mapping_Classification.tab.tmp
mv Mapping_Classification.tab.tmp Mapping_Classification.tab

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


sort -k1,1 -k2,2 -n -u -o DELETIONS.tab DELETIONS.tab


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


process MUT_GATHER{
cpus 1
memory "10GB"
input:
        path(tabs)
output:
        path('all_cf1N.tab')
script:
"""
#!/usr/bin/env Rscript
library(tidyverse)
files <- list.files(path = ".", pattern = "corrected.tab")
lapply(files, function(x){
    read_delim(x,delim='\\t') %>%
    mutate(File=str_remove_all(x,'_.*'))%>%
    relocate(File)
})->l
bind_rows(l)%>%
write_delim('all_cf1N.tab',delim='\\t',col_names = F)
"""
}



process PHARMA {
cpus 1
memory "10GB"
publishDir "Called", mode: 'copy', pattern: 'pharma_gene*'
publishDir "OUTPUT", mode: 'copy', pattern: '*format*'
input:
        path (tabs)
        val (thresholds)
        path (genes)

output:
        path("pharma_gene*"), emit:pgene
        path("*format*"), emit:format
script:
"""
#!/usr/bin/env Rscript
library(tidyverse)

t=${thresholds}

# Read in the data
#read all the file that ends with "corrected.tab"
files <- list.files(path = ".", pattern = "corrected.tab")

lapply(files, function(x){
    read_delim(x,delim='\\t') %>%
    filter(Freq>=t) %>%
    mutate(ID=str_remove_all(x,'_.*'))
})->l

genes<-read_delim('gene_drug.csv',col_names=c('Gene_name','Start','Stop','Gene','S'))

#from genes get a vector containing all the numbers in the intervavals between Start and Stop

lapply(genes\$Start,function(x){
    seq(x,genes\$Stop[genes\$Start==x][1])
})->l2

unlist(l2)->gene_intervals



bind_rows(l) %>%
    filter(`#Pos` %in% gene_intervals) %>%
    relocate(ID)->p


genes %>% 
  group_by(Gene_name) %>% 
  mutate(i=str_c(seq(Start,Stop),collapse=',')) %>% 
  separate_rows(i,sep=',') %>% 
  transmute(Name=Gene_name,
            Range=ifelse(Start>Stop,
                         paste('r[',Stop,'-',Start,']',sep=''),
                         paste('r[',Start,'-',Stop,']',sep='')),
            Locus=Gene,`#Pos`=as.numeric(i),
            Start=ifelse(S=='rev',max(Start,Stop),min(Start,Stop))) %>% 
  right_join(p) %>% 
  select(Name,Range,Locus,`#Pos`,ID,Start,Freq,Ref,Allel,Type,Subst) %>%
  mutate(
    Type=ifelse(Type=='SNP',Type,paste('__',Type,sep='')),
    Mut=case_when(str_detect(Name,'_ups')~paste0(Ref,'-',abs(`#Pos`-Start),Allel,'-',`#Pos`),
                       (Subst!=' ')~str_remove(Subst,' .*'),
                  (Type=='SNP')~paste0(`#Pos`,'_',Ref,'>',Allel,'_',Type),
                  (Type!='SNP')~paste0(`#Pos`,'_',str_to_upper(Type))),
    Mut=ifelse(Freq<75,paste0(Mut,':',Freq),Mut)) %>% 
  ungroup() %>% 
  select(Name,Range,LocusName=Locus,ID,Mut) %>% 
  group_by(Name,Range,LocusName,ID) %>% 
  summarise(Mut=str_c(Mut,collapse=' ')) %>%ungroup() %>%  
  select(Name,ID,Mut) %>% 
  spread(Name,Mut) %>% rename(Name=ID)->tail


genes %>% 
  transmute(Name=Gene_name,
            Range=ifelse(Start>Stop,
                         paste('r[',Stop,'-',Start,']',sep=''),
                         paste('r[',Start,'-',Stop,']',sep='')),
            LocusName=Gene)->head


as.data.frame(t(head[,-1]))->h
colnames(h)<-head\$Name
h %>% rownames_to_column('Name') %>% as_tibble()->head



bind_rows(list(head,tail)) %>% 
  write_delim(paste0('pharma_gene_format_',t,'.tab'),delim='\\t',na = '')

p%>%
write_delim('pharma_gene.tab',delim='\\t',col_names=F)


"""
}



process  WHO {
cpus 1
memory "10GB"
input:
        path(tab)
        path(head)
        path(WHO)
output:
        path('WHO_raw.csv')
script:
"""
#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from functools import reduce


#Check wheter a subs is Syn
def syn(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res = temp.match(s).groups()
    if res[0] == res[2]:
        print(s)

#Check wheter a subs is NoSyn
def nsyn(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res = temp.match(s).groups()
    if res[0] != res[2]:
        print(s)
        
def nsyn2(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res1 = temp.match(s).group(1)
    res2 = temp.match(s).group(2)
    res3 = temp.match(s).group(3)
    if res1 != res3:
        print(s)
        
# Check whether there is STOP codon       
def nsyn3(s):
    import re
    temp = re.compile("([a-zA-Z_]+)([0-9]+)([a-zA-Z_]+)")
    #check whether s IS NOT white space
    if s and not s.isspace():
        res1 = temp.match(s).group(1)
        res2 = temp.match(s).group(2)
        res3 = temp.match(s).group(3)
        if res1 != res3:
            return(True)
        else:
            return(False)
    else:
        return(False)

# remove content between () in s
def cleanvar(s):
    import re
    return(re.sub(r"\\([^()]*\\)", "", s))  

# convert tab to space
def tab2space(s):
    import re
    return(s.replace('\\t', ' '))

def removespace(s):
    import re
    return(s.replace(' ', ''))

def remove_(s):
    S=s.split('_')
    return(S[0])

def remove1_(s):
    S=s.split('_')
    return(S[1])

def first(s):
    return(s[0])

def last(s):
    return(s[-1])

def one2tree(s):
    from Bio.PDB.Polypeptide import one_to_three
    s=s.rstrip()
    if s.isupper() and not s.isnumeric() :
        try:
            return(one_to_three(s))
        except:
            return(s)
    else:
        return(s)
    
def codon(s):
    import re
    temp = re.findall(r'\\d+', s)
    return(temp[0])


myfile = open("head", "r")
head1 = myfile.read()
head_list1=head1.split('\\t')
head_list1[16]=head_list1[16].strip()
head_list1.insert(0,"File")
all=pd.read_csv('all_cf1N.tab',sep='\\t',header=None,names=head_list1)


#clean substitution and add a new column SubstClean
all['SubstClean']=all.Subst.apply(cleanvar)
all['SubstClean2'] = all.SubstClean.apply(removespace)
#all['ID']=all.File.apply(remove_)
all['ID']=all.File
#all['NSYN'] = all.SubstClean2.apply(nsyn3)
all.rename({'#Pos': 'Genome position'}, axis=1,inplace=True)
all['genome_index']=all['Genome position'].astype(int)
all.genome_index.astype(int)
whoG=pd.read_csv('${WHO}',sep='\\t')

whoG.genome_index=whoG.genome_index.astype(str).str.split(',')
whoG = whoG.explode('genome_index').reset_index(drop=True)
whoG.genome_index=whoG.genome_index.astype(int)
whoG['Allel']=whoG['final_annotation.AlternativeNucleotide'].str.upper()
whoG['Ref']=whoG['final_annotation.ReferenceNucleotide'].str.upper()
whoG['common_name']=whoG['final_annotation.Gene'] + '_' + whoG['final_annotation.TentativeHGVSNucleotidicAnnotation']

all_whoG=pd.merge(all,whoG,on=['genome_index','Allel','Ref'])

all_whoG.sort_values('File',inplace=True)


all_whoG.filter(regex=r'(File|_Conf_Grade)')
drugs_conf=all_whoG.filter(regex=r'(_Conf_Grade)').columns
# create array of drugs
drugs=drugs_conf.str.split('_').str[0]

cutoff=10
A = {}
Asubst = {}
Asubst1 = {}
Asubst5 = {}
Ali=[]
Alisubst=[]
Alisubst5=[]
Alisubst1=[]
for i in drugs_conf:
    t = i.split('_')[0]
    A[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File',i]]
    Asubst[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    Asubst1[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim') | (all_whoG[i]=='3) Uncertain significance')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    Asubst5[t] = all_whoG[(all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    Ali.append(A[t])
    Alisubst.append(Asubst[t])
    Alisubst5.append(Asubst5[t])
    Alisubst1.append(Asubst1[t])


A5=pd.concat(Alisubst5)

A5.fillna(0,inplace=True)


#A5a=A5.drop('Freq',axis=1).drop_duplicates()
A5a=A5.drop_duplicates()
#A5a['variant']=A5a['variant'] + '_FR' + A5a['Freq'].astype(str)
A5a['variant']=np.where(A5a['Freq'] > 75, A5a['variant'],A5a['variant'] + ':' + A5a['Freq'].astype(str))
A5b=A5a.groupby(['File','RIF_Conf_Grade','INH_Conf_Grade', 'EMB_Conf_Grade', 'PZA_Conf_Grade', 'LEV_Conf_Grade',       'MXF_Conf_Grade', 'BDQ_Conf_Grade', 'LZD_Conf_Grade', 'CFZ_Conf_Grade',       'DLM_Conf_Grade', 'AMI_Conf_Grade', 'STM_Conf_Grade', 'ETH_Conf_Grade',       'KAN_Conf_Grade', 'CAP_Conf_Grade'])['variant'].apply(', '.join).reset_index()


#A5b=A5.groupby(['File','RIF_Conf_Grade','INH_Conf_Grade', 'EMB_Conf_Grade', 'PZA_Conf_Grade', 'LEV_Conf_Grade',       'MXF_Conf_Grade', 'BDQ_Conf_Grade', 'LZD_Conf_Grade', 'CFZ_Conf_Grade',       'DLM_Conf_Grade', 'AMI_Conf_Grade', 'STM_Conf_Grade', 'ETH_Conf_Grade',       'KAN_Conf_Grade', 'CAP_Conf_Grade'])['variant'].apply(', '.join).reset_index()


A5b.insert(1, "ID", A5b.File.str.split('_').str[0], True)
#T5b['ID']=T5b.File.str.split('_').str[0]
#A5b['new'] = A5b['ID'].str.extract('(\\d+)').astype(int)
#https://stackoverflow.com/questions/66134896/python-pandas-sort-an-alphanumeric-dataframe
#A5b = A5b.sort_values(by=['new'], ascending=True).drop('new', axis=1)
A5b = A5b.sort_values(by=['ID'], ascending=True)#.drop('new', axis=1)

A5b.to_csv('WHO_raw.csv',index=False, sep=';')

"""
}

process OUT_WHO {
cpus 1
memory "20G"
publishDir "OUTPUT", mode: 'copy', pattern: 'res_WHO.csv'
input:
        path(WHO)
output:
        path('res_WHO.csv')
script:
"""
#!/usr/bin/env Rscript
library(tidyverse)

read_delim('${WHO}',delim=';') %>%
  select(ID) %>% distinct() %>%
  left_join(
    (read_delim('${WHO}',delim=';') %>%
       select(ID, ends_with('Grade'),variant) %>%
       gather(drug,grade,ends_with('Grade')) %>%
       filter(grade!='0') %>%
       mutate(drug=str_remove_all(drug,'_.*')) %>%
       group_by(ID,drug,grade) %>%
       summarise(v=str_c(variant,collapse=', ')))) %>%
  mutate(grade=str_replace_all(grade, ' ','_')) %>%
  group_by(ID) %>% 
  complete(drug,grade) %>% 
  ungroup() %>% 
  unite('drug',drug,grade,sep='_') %>%
  spread(drug,v) %>%
  mutate(Interpretation_126=case_when((!if_all(matches('^RIF_1|^RIF_2|^RIF_6'),is.na)&
                                       !if_all(matches('^INH_1|^INH_2|^INH_6'),is.na)&
                                       !if_all(matches('^LEV_1|^LEV_2|^LEV_6|^MXF_1|^MXF_2|^MXF_6'),is.na)&
                                       !if_all(matches('^BDQ_1|^BDQ_2|^BDQ_6'),is.na))~'XDR',
                                    (!if_all(matches('^RIF_1|^RIF_2|^RIF_6'),is.na)&
                                       !if_all(matches('^INH_1|^INH_2|^INH_6'),is.na)&
                                       !if_all(matches('^LEV_1|^LEV_2|^LEV_6|^MXF_1|^MXF_2|^MXF_6'),is.na))~'Pre-XDR',
                                    (!if_all(matches('^RIF_1|^RIF_2|^RIF_6'),is.na)&
                                       !if_all(matches('^INH_1|^INH_2|^INH_6'),is.na))~'MDR',
                                    (!if_all(matches('^RIF_1|^RIF_2|^RIF_6'),is.na))~'RR')) %>% 
  select(ID,starts_with('Interpretation'),matches('\\\\)')) %>% 
  write_delim('res_WHO.csv',delim=';',na='')


"""

}