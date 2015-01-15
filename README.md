# TRACE
Trace Primer Design Software
TRACE Primer-Linker Design Documentation v 1.0
Last Updated: January 2014

The TRACE software is a standalone MATLAB application that can design orthogonal primers and linkers. This software can be downloaded here.



## Step 1: Designing Primers from Search Space:

Sequences that are upstream and downstream of the sites of interest can be input into the table labelled genome sites in the Step 1 table (5'-3').  This search space should be from the same strand on the genome and should flank the site:

Upstream Search Space (5'-3') - Site - Downstream Search Space (5'-3')

Any number of sites (1-10) can be chosen.  Fill sites in sequential order.  Primers can also be loaded from a .csv file with column 1 as site name, column 2 as upstream sequence and column 3 as downstream sequence.  Rows correspond to different sites (row 1 -> site 1).  After inputting search space click Calculate Primers once.



#### User Defined Variables:

Primer Melting Temperature: The melting temperature primers are designed for. Default = 60.

Homodimer Assoc (max):  The maximum association of a primer with its reverse complement. Reduce this to improve primer selection stringency but reduce search flexibility. Default = 7.

Nucleotides from 3' end to check: Homology is calculate based on 3' homology.  The more nucleotides calculated for, the more stringent design space.  Minimum primer size is set to 19 nt, so please keep this value 10 or below. Default = 12.



## Step 2: Designing Linkers onto Primers.

Linkers are randomly chosen and tested for compatibility with the primer set designed from Step 1.  The software searches for N-1 linkers.  The program dynamically updates the number of random linker sequences tested. Once primers are designed and linker parameters are chosen click "Calculate Linkers" once.



#### User Defined Variables:

Linker Length: The length of the linker search space. Linker length, construct density and melt temperature are related. Default = 29 nt.

Linker Tm: The linker melt temperature.  Higher melt temperatures require longer linkers but have more stable linker-linker interactions. Default = 70.



Maximum Association: The maximum association primer-linker sets can have.  If maximum association is set lower than the homodimer association or heterodimer association then no solution will be found. Default = 8.



#### Final Sequences with Linkers:

Primer-Linkers can be exported to a .csv file with the "Export" command.



## TroubleShooting:

No acceptable primers found for site N: Increase search space or loosen primer design constraints.
