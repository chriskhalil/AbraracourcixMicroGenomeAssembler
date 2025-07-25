## AbraracourcixMicroGenomeAssembler
  A high-performance bioinformatics tool for approximate pattern matching and genome assembly, developed as part of the AUB Bioinformatics Graduate course.
   Uses BWT for seeding/pattern matching and a greedy score based aligner.

## SYNOPSIS
   **random_assembler** [-ms <match_score>] [-mp <mismatch_penalty>] [-gp <gap_penalty>] [-ep <extended_penalty>] [-vb] [-mascot] -irf <input_reference> -ird <input_reads> -m <allowed_mismatches> -t <num_threads>

## DESCRIPTION
   **Random Assembler** is an ultrafast and memory-efficient tool for assebmling reads to long reference sequences.
   **bwt_storm** is a faster and more optimized approach based on our random_assembler.
    
    Note: bwt_storm has not been test as much as the base random_assembler as it was done after the competition as training and fun project.
## OPTIONS
   ### Mandatory Arguments:
   -irf <input_reference>
          Specify the input reference file.

   -ird <input_reads>
          Specify the input reads file.

   -m <allowed_mismatches>
          Specify the allowed number of mismatches.

   -t <num_threads>
          Specify the number of threads to be used.

   ### Optional Arguments:
   -ms <match_score>
          Specify the match score. (Default: 1)

   -mp <mismatch_penalty>
          Specify the mismatch penalty. (Default: 4)

   -gp <gap_penalty>
          Specify the gap penalty. (Default: 6)

   -ep <extended_penalty>
          Specify the extended penalty. (Default: 1)

   -vb
          Enable verbose mode.

   -mascot
          Print mascot.

## AUTHORS
   This program, **abraracourcix**, was created by Christophe Khalil (cak29@mail.aub.edu) and Ali Chehab (amc33@mail.aub.edu) on 20/11/2023.

## LICENSE
   This program is licensed under the MIT License. If you use this program, please cite the authors.

## REPORTING BUGS
   Report bugs to cak29@mail.aub.edu and amc33@mail.aub.edu.

## SEE ALSO
   More information can be found in the report.pdf.

## EXAMPLES
   ```
   random_assebmler -irf reference.txt -ird reads.txt -m 2 -t 8 -vb
   ```
---

**compile flags**:
-std=c++17 -fopenmp -O2

Assumptions:
- inputs are in text format and contain no symbols
- reads do not come with a quality score (no quality score provided)
- single stranded input
