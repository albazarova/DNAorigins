# DNAorigins
Supplementary code for paper "Bayesian inference of origin firing time distributions, origin interference and licensing probabilities from NGS data".

The repository consists of two .cpp source files essential for running an algorithm: algorithm_with_NA.cpp, algorithm_plain.cpp  

algorithm_with_NA.cpp is designed to handle the data with unsequencable region, algorithm_plain.cpp - without it.

The files are to be compiled with C++11

Processed data files are provided. Where not stated otherwise data are obtained from paired-end NGS Okazaki Fragment data

chrX_f/r.txt

X - chromosome number, f corresponds to the data on the forward strand, r - to the data on the reverse strand

left/right indicates whether the data corresponds to the left or right triplet as described in the paper

For chromosome 10 numbers XYZ correspond to the origins included in the triplet; wt/rat1 corresponds to the datasets analysed from single-end sequenced data. Wt corresponds to wt, rat1 - to temperature sensitive mutant;

chr2 ARS207,207.5,207.8
chr7_left ARS717, 718, 719
chr7_right ARS718, 719, 720
chr8_left ARS813, 815, 816
chr8_right ARS815, 816, 818
chr10_123 ARS1010, 1011, 1013
chr10_234 ARS1011, 1013, 1014
chr10_345 ARS1013, 1014, 1015
chr10_456 ARS1014,1015,1018 missing values between 1098 and 1336 for paired-end sequenced data, 1108 and 1336 for wt/rat1
chr10_567 ARS1015, 1018, 1019 missing values 588 and 826 for paired-end sequenced data, 605 and 833 for wt/rat1
chr10_678 ARS1018, 1019, 1021



