This directory is focused on testing for shared divergence time among taxa that span the rivers of Mainland Southeast Asia. To do this, we use the phylogenetics software tool ecoevolity (https://phyletica.org/ecoevolity/index.html).

In short, we took the alignment files that indicated riverine divergence and used those to test for shared divergence times.

Ecoevolity analyses alignment files in nexus format. Files should be formatted with the population label at the end of the individual name, and population labels cannot be shared across alignment files. Because our fasta alignments all contain labels such as "E" (east) and "W" (west) for population labels, we simply change these to "E1" and "W1" for one alignment file, "E2" and "W2" for other alignment files, and so on.

We then created the config.yml file, which contains the parameters for running ecoevolity. To run the analysis, we performed the following command:

```
ecoevolity --relax-triallelic-sites --relax-missing-sites config.yml
```

We ran two chains to assess convergence and boost posterior sample sizes. 
