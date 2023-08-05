# NF_TBSEQ

Nextflow implementation of MTBseq pipeline with additional features.

## Requirements
  - Singularity
  - Nextflow

## Quick Usage
Put the .fastqs in a folder named "samp" and then use:

`nextflow run https://github.com/Allen13x/NF_TBSEQ -with-singularity library://allen13x/mtbseq/nf_mtbseq:1.0.0 -resume --join false --extra false -r main --reads $(pwd)/samp/`

This will perform the basic single-sample-level characterizion from MTBseq pipeline

## Parameters

```
  -h          Show this help
  -ref        Reference Genome to use:
                  M._abscessus_CIP-104536T_2014-02-03
                  M._chimaera_DSM44623_2016-01-28
                  M._fortuitum_CT6_2016-01-08
                  M._tuberculosis_H37Rv_2015-11-13 [default]
  -join       Perform joint analysis [default: true]
  -sj         List of samples to use for joint analysis [optional]
  -proj       Name of the project for joint analysis [default: def]
  -pharma     Perform Drug resistance analysis at Custom frequencies (to use after first round of analysis) [default: false]
  -tdrug      Custom frequencies % [default: 10]
  -extra      Perform extra analysis [default: true]:
                 Deletion/Insertion detection with Delly2
                 Sequencing depth with Mosdepth
                 Pharma analysis at 10% and Comparison with WHO catalogue            
```
