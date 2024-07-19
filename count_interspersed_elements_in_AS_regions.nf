#!/usr/bin/env nextflow

params.input_cenSat_bed="http://public.gi.ucsc.edu/~hloucks/CenSat/CHM13/CHM13-noHOR_alphaSummary_detailed_mixed.sorted.bed"
params.input_RM_bigBed="http://t2t.gi.ucsc.edu/chm13/dev/t2t-chm13-v2.0/rmsk/rmsk.bigBed"
params.output="./output"
cenSat_bed=file("${params.input_cenSat_bed}").getBaseName()
RM_bigBed=file("${params.input_RM_bigBed}").getBaseName()

process download_input_bed {
  input:
    path input_cenSat_bed
    path input_RM_bigBed
  output:
    tuple path('cenSat.bed'),
          path('RM.bigBed')
  shell:
  '''
  wget -O cenSat.bed !{params.input_cenSat_bed}
  wget -O RM.bigBed !{params.input_RM_bigBed}
  '''
}

process bigBed_to_bed {
  input:
    tuple path('cenSat.bed'), path('RM.bigBed')
  output:
    tuple path('cenSat.bed'), path('RM.bed')
  shell:
  '''
#  module add kent-tools
  bigBedToBed RM.bigBed RM.bed
  '''
}

process filter_bed {
//  publishDir "$params.output", mode: 'link', pattern: '*.bed'
  input:
    tuple path('cenSat.bed'), path('RM.bed')
  output:
    tuple path('cenSat.SF4.bed'),
          path('cenSat.SF5.bed'),
          path('cenSat.SF6.bed'),
          path('RM.wo_ALR_Alpha.bed')
  shell:
  '''
  fgrep 'SF4' cenSat.bed > cenSat.SF4.bed
  fgrep 'SF5' cenSat.bed > cenSat.SF5.bed
  fgrep 'SF6' cenSat.bed > cenSat.SF6.bed
  fgrep -v 'ALR_Alpha' RM.bed > RM.wo_ALR_Alpha.bed
  '''
}

process intersect_bed {
  publishDir "$params.output", mode: 'link', pattern: '*.bed'
  input:
    tuple path('cenSat.SF4.bed'),
          path('cenSat.SF5.bed'),
          path('cenSat.SF6.bed'),
          path('RM.wo_ALR_Alpha.bed')
  output:
    tuple path('int_RM_cenSat.SF4.bed'),
          path('int_RM_cenSat.SF5.bed'),
          path('int_RM_cenSat.SF6.bed')
  shell:
  '''
#  module add BEDTools
  grep -E 'L1PA3|L1PA2|L1HS' RM.wo_ALR_Alpha.bed | \
  bedtools intersect -a - -b cenSat.SF4.bed -f 0.1 -wa -wb | \
  sort -k1,1 -k2,2n > int_RM_cenSat.SF4.bed
  grep -E 'L1PA2|L1HS' RM.wo_ALR_Alpha.bed | \
  bedtools intersect -a - -b cenSat.SF5.bed -f 0.1 -wa -wb | \
  sort -k1,1 -k2,2n > int_RM_cenSat.SF5.bed
  grep -E 'L1PA4|L1PA3|L1PA2|L1HS' RM.wo_ALR_Alpha.bed | \
  bedtools intersect -a - -b cenSat.SF6.bed -f 0.1 -wa -wb | \
  sort -k1,1 -k2,2n > int_RM_cenSat.SF6.bed
  '''
}

process output_stats {
  publishDir "$params.output", mode: 'link', pattern: '*.txt'
  input:
    tuple path('int_RM_cenSat.SF4.bed'),
          path('int_RM_cenSat.SF5.bed'),
          path('int_RM_cenSat.SF6.bed')
  output:
    file("*txt")
  shell:
  '''
  awk -v OFS="\t" '{a[\$4]++;n++}END{for(i in a)print i, a[i]; print "Total", n}' \
  int_RM_cenSat.SF4.bed > L1_in_SF4.txt
  awk -v OFS="\t" '{a[\$4]++;n++}END{for(i in a)print i, a[i]; print "Total", n}' \
  int_RM_cenSat.SF5.bed > L1_in_SF5.txt
  awk -v OFS="\t" '{a[\$4]++;n++}END{for(i in a)print i, a[i]; print "Total", n}' \
  int_RM_cenSat.SF6.bed > L1_in_SF6.txt
  '''
}

workflow {
  download_input_bed(params.input_cenSat_bed,
                     params.input_RM_bigBed)
  bigBed_to_bed(download_input_bed.out)
  filter_bed(bigBed_to_bed.out)
  intersect_bed(filter_bed.out)
//  search_closest(intersect_bed.out)
  output_stats(intersect_bed.out)
}
