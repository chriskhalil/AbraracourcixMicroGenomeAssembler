# Abraracourcix - A small genome assembler
This is our implementation of a reference based genome assembler.
Uses BWT for seeding/pattern matching and a greedy score based aligner.

Assumptions:
- inputs are in text format and contain no symbols
- reads do not come with a quality score
- single stranded input
