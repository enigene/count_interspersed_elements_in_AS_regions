#!/usr/bin/env nextflow

params.input_cenSat_bed="http://public.gi.ucsc.edu/~hloucks/CenSat/CHM13/CHM13-noHOR_alphaSummary_detailed_mixed.sorted.bed"
params.input_RM_bigBed="http://t2t.gi.ucsc.edu/chm13/dev/t2t-chm13-v2.0/rmsk/rmsk.bigBed"
params.input_ASat_SF_bigBed="http://t2t.gi.ucsc.edu/chm13/dev/t2t-chm13-v2.0/as_annotation/ASat_SF.bigBed"
params.output="./output"
cenSat_bed=file("${params.input_cenSat_bed}").getBaseName()
RM_bigBed=file("${params.input_RM_bigBed}").getBaseName()
ASat_SF_bigBed=file("${params.input_ASat_SF_bigBed}").getBaseName()

process download_input_bed {
  input:
    path input_cenSat_bed
    path input_RM_bigBed
    path input_ASat_SF_bigBed
  output:
    tuple path('cenSat.bed'),
          path('RM.bigBed'),
          path('ASat_SF.bigBed')
  shell:
  '''
  wget -O cenSat.bed !{params.input_cenSat_bed}
  wget -O RM.bigBed !{params.input_RM_bigBed}
  wget -O ASat_SF.bigBed !{params.input_ASat_SF_bigBed}
  '''
}

process bigBed_to_bed {
  input:
    tuple path('cenSat.bed'), path('RM.bigBed'), path('ASat_SF.bigBed')
  output:
    tuple path('cenSat.bed'), path('RM.bed'), path('ASat_SF.bed')
  shell:
  '''
#  module add kent-tools
  bigBedToBed RM.bigBed RM.bed
  bigBedToBed ASat_SF.bigBed ASat_SF.bed
  '''
}

process filter_bed {
  input:
    tuple path('cenSat.bed'), path('RM.bed'), path('ASat_SF.bed')
  output:
    tuple path('cenSat.SF4.wo_mixed.bed'),
          path('cenSat.SF5.wo_mixed.bed'),
          path('cenSat.SF6.wo_mixed.bed'),
          path('ASat_SF.only_Ga.bed'),
          path('ASat_SF.only_R1R2.bed'),
          path('ASat_SF.only_Ha.bed'),
          path('RM.wo_ALR_Alpha.bed'),
          path('cenSat.only_mons.bed')
  shell:
  '''
  fgrep 'SF4' cenSat.bed | fgrep -v 'mix' > cenSat.SF4.wo_mixed.bed
  fgrep 'SF5' cenSat.bed | fgrep -v 'mix' > cenSat.SF5.wo_mixed.bed
  fgrep 'SF6' cenSat.bed | fgrep -v 'mix' > cenSat.SF6.wo_mixed.bed
  fgrep -v 'ALR_Alpha' RM.bed > RM.wo_ALR_Alpha.bed
  fgrep 'Ga' ASat_SF.bed > ASat_SF.only_Ga.bed
  grep -E 'R1|R2' ASat_SF.bed > ASat_SF.only_R1R2.bed
  fgrep 'Ha' ASat_SF.bed > ASat_SF.only_Ha.bed
  fgrep 'mon' cenSat.bed > cenSat.only_mons.bed
  '''
}

process intersect_bed {
  input:
    tuple path('cenSat.SF4.wo_mixed.bed'),
          path('cenSat.SF5.wo_mixed.bed'),
          path('cenSat.SF6.wo_mixed.bed'),
          path('ASat_SF.only_Ga.bed'),
          path('ASat_SF.only_R1R2.bed'),
          path('ASat_SF.only_Ha.bed'),
          path('RM.wo_ALR_Alpha.bed'),
          path('cenSat.only_mons.bed')
  output:
    tuple path('int_RM_cenSat.SF4.bed'),
          path('int_RM_cenSat.SF5.bed'),
          path('int_RM_cenSat.SF6.bed'),
          path('ASat_SF.only_Ga_in_mons.bed'),
          path('ASat_SF.only_R1R2_in_mons.bed'),
          path('ASat_SF.only_Ha_in_mons.bed')
  shell:
  '''
#  module add BEDTools
  bedtools intersect -a RM.wo_ALR_Alpha.bed -b cenSat.SF4.wo_mixed.bed \
  -f 0.1 -wa -wb | grep -E 'L1PA3|L1PA2|L1HS' | sort -k1,1 -k2,2n > int_RM_cenSat.SF4.bed
  bedtools intersect -a RM.wo_ALR_Alpha.bed -b cenSat.SF5.wo_mixed.bed \
  -f 0.1 -wa -wb | grep -E 'L1PA2|L1HS' | sort -k1,1 -k2,2n > int_RM_cenSat.SF5.bed
  bedtools intersect -a RM.wo_ALR_Alpha.bed -b cenSat.SF6.wo_mixed.bed \
  -f 0.1 -wa -wb | grep -E 'L1PA4|L1PA3|L1PA2|L1HS' | sort -k1,1 -k2,2n > int_RM_cenSat.SF6.bed

  bedtools intersect -a ASat_SF.only_Ga.bed -b cenSat.only_mons.bed \
  -f 0.1 -wa | sort -k1,1 -k2,2n > ASat_SF.only_Ga_in_mons.bed
  bedtools intersect -a ASat_SF.only_R1R2.bed -b cenSat.only_mons.bed \
  -f 0.1 -wa | sort -k1,1 -k2,2n > ASat_SF.only_R1R2_in_mons.bed
  bedtools intersect -a ASat_SF.only_Ha.bed -b cenSat.only_mons.bed \
  -f 0.1 -wa | sort -k1,1 -k2,2n > ASat_SF.only_Ha_in_mons.bed
  '''
}

process search_closest {
  publishDir "$params.output", mode: 'link', pattern: 'closest*.bed'
  input:
    tuple path('int_RM_cenSat.SF4.bed'),
          path('int_RM_cenSat.SF5.bed'),
          path('int_RM_cenSat.SF6.bed'),
          path('ASat_SF.only_Ga_in_mons.bed'),
          path('ASat_SF.only_R1R2_in_mons.bed'),
          path('ASat_SF.only_Ha_in_mons.bed')
  output:
    tuple path('closest_L1_to_Ga.bed'),
          path('closest_L1_to_R1R2.bed'),
          path('closest_L1_to_Ha.bed'),
          path('ASat_SF.only_Ga_in_mons.bed'),
          path('ASat_SF.only_R1R2_in_mons.bed'),
          path('ASat_SF.only_Ha_in_mons.bed')
  shell:
  '''
#  module add BEDTools
  bedtools closest -d -a ASat_SF.only_Ga_in_mons.bed -b int_RM_cenSat.SF4.bed -t all | \
  awk '\$10!="."&&+\$NF<500' | \
  awk -v FS="\\t" '{for(i=10;i<=22;i++){if(i<22){printf("%s\\t",\$i)}else{printf("%s\\n",\$i)}}}' | \
  sort -u -k1,1V -k2,2n | \
  awk -v OFS="\\t" '{print \$1,\$2,\$3,\$4,\$6,\$3-\$2}' > closest_L1_to_Ga.bed

  bedtools closest -d -a ASat_SF.only_R1R2_in_mons.bed -b int_RM_cenSat.SF5.bed -t all | \
  awk '\$10!="."&&+\$NF<500' | \
  awk -v FS="\\t" '{for(i=10;i<=22;i++){if(i<22){printf("%s\\t",\$i)}else{printf("%s\\n",\$i)}}}' | \
  sort -u -k1,1V -k2,2n | \
  awk -v OFS="\\t" '{print \$1,\$2,\$3,\$4,\$6,\$3-\$2}' > closest_L1_to_R1R2.bed

  bedtools closest -d -a ASat_SF.only_Ha_in_mons.bed -b int_RM_cenSat.SF6.bed -t all | \
  awk '\$10!="."&&+\$NF<500' | \
  awk -v FS="\\t" '{for(i=10;i<=22;i++){if(i<22){printf("%s\\t",\$i)}else{printf("%s\\n",\$i)}}}' | \
  sort -u -k1,1V -k2,2n | \
  awk -v OFS="\\t" '{print \$1,\$2,\$3,\$4,\$6,\$3-\$2}' > closest_L1_to_Ha.bed
  '''
}

process output_stats {
  publishDir "$params.output", mode: 'link', pattern: '*.txt'
  input:
    tuple path('closest_L1_to_Ga.bed'),
          path('closest_L1_to_R1R2.bed'),
          path('closest_L1_to_Ha.bed'),
          path('ASat_SF.only_Ga_in_mons.bed'),
          path('ASat_SF.only_R1R2_in_mons.bed'),
          path('ASat_SF.only_Ha_in_mons.bed')
  output:
    file("*txt")
  shell:
  '''
  awk -v OFS="\t" '{a[\$4]++;n++}END{for(i in a)print i, a[i]; print "Total", n}' \
  closest_L1_to_Ga.bed > closest_L1_to_Ga.txt
  awk -v OFS="\t" '{a[\$4]++;n++}END{for(i in a)print i, a[i]; print "Total", n}' \
  closest_L1_to_R1R2.bed > closest_L1_to_R1R2.txt
  awk -v OFS="\t" '{a[\$4]++;n++}END{for(i in a)print i, a[i]; print "Total", n}' \
  closest_L1_to_Ha.bed > closest_L1_to_Ha.txt

  awk '{l+=\$3-\$2;n++;a[\$1]++}END{print \
  "The length of regions: " l " bp, \
  total regions: " n ", on " length(a) " chromosomes"}' \
  ASat_SF.only_Ga_in_mons.bed > ASat_SF.only_Ga_in_mons.txt

  awk '{l+=\$3-\$2;n++;a[\$1]++}END{print \
  "The length of regions: " l " bp, \
  total regions: " n ", on " length(a) " chromosomes"}' \
  ASat_SF.only_R1R2_in_mons.bed > ASat_SF.only_R1R2_in_mons.txt

  awk '{l+=\$3-\$2;n++;a[\$1]++}END{print \
  "The length of regions: " l " bp, \
  total regions: " n ", on " length(a) " chromosomes"}' \
  ASat_SF.only_Ha_in_mons.bed > ASat_SF.only_Ha_in_mons.txt
  '''
}

workflow {
  download_input_bed(params.input_cenSat_bed,
                     params.input_RM_bigBed,
                     params.input_ASat_SF_bigBed)
  bigBed_to_bed(download_input_bed.out)
  filter_bed(bigBed_to_bed.out)
  intersect_bed(filter_bed.out)
  search_closest(intersect_bed.out)
  output_stats(search_closest.out)
}
