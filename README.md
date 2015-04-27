# Population genetics and paternity analysis for moss populations


These scripts were used to determine multilocus microsatellite genotypes, assign haploid paternity to diploid offspring, and estimate mating pattern statistics in moss populations.

# `clonecheck.py`

##Description

In moss populations, there can be considerable asexual "clonal" reproduction. By collecting many individual "stems" (ramets), the same genetic individual (genet) may be sampled multiple times. This script will return a list of unique multilocus genotypes, as well as a file assigning each haploid individual to one or more clones.

This script features a method for dealing with missing data and the potential for one or more mismatched loci in assigning samples to clones. First, the script generates unique multilocus genotypes from all individuals with no missing data. All other individuals are assigned to these genotypes unless they have unique alleles, in which case a new multilocus genotype is created. 

If the number of "mismatched loci" is set to more than zero, then an individual can be assigned to a clone even if it does not match at every locus. 

An individual may be assigned to more than one clone if there is missing data or the number of mismatched loci is set to one or greater.

##Input

The file must be a comma-delimited file with no header. 
On each row, the first field should be an individual ID, followed by the multilocus genotype:

`MJ874,158,193,191,254,194,0,249,129,276,180,290,195,198`

Indicate missing data with a `0` 

For an example input file, see: `examples/tenellum_paternity.csv`

##Output

The file stem (before `.csv`) will be used to write the output files.

###`filestem_clones.csv`
	
Contains only the unique multilocus genotypes, one per line, in comma-delimited format.

###`filestem_assign.csv`

Contains the original genotype information for each individual, plus additional columns assigning each individual to one or more clones. 

Individuals that can be unambiguously sorted into just one clone will have only one additional column.

##Example Usage

```
python clonecheck.py examples/tenellum_paternity.csv 0
```

This will determine the clones with no allowance for mismatches between individuals and clones.

# `paternity.py`

##Description
In most mosses (excluding polyploid species), the paternity of sporophytes can be determined by simple subtraction. Sporophytes remain attached to their haploid maternal gametophyte, and thus for most species can be sampled simultaneously. Any alleles that appear in the sporophyte that are not present in the maternal gametophyte must be of paternal origin. 

This script determines the haploid paternity of diploid offspring using the subtraction method. For heterozygotes, the paternal allele is the one not present in the maternal individual. For homozygotes, the paternal alelle is assumed to be the same as teh maternal allele. This method should work for any codominant genetic marker (such as microsatellites or SSRs).

##Input
You must provide the multilocus genotypes of the haploid maternal parent and the diploid offspring.

The file is a comma-separated-value (CSV) format.
The first line of the document should be comma-separated list of loci names, beginning with "Loci":

	```Loci,17a,65a,78a,1a,7a,68a,10a,18a,19a,9a,20a,14a,56a```

The multilocus genotypes begin on the next line. The first field is the individual ID. 
The second field is the ID of the maternal parent.
Then there should be two columns per locus. For maternal genotypes, simply list each allele twice.

###Example Maternal Genotype:

```MJ872,,158,158,193,193,191,191,254,254,190,190,219,219,249,249,129,129,273,273,180,180,290,290,195,195,198,198```


Note that the second field is blank. This identifies this individual as a maternal haploid.

###Example Offspring Genotype:

```MJ873,MJ872,158,158,193,193,191,191,254,254,182,190,219,219,249,249,129,129,273,276,180,180,290,293,195,195,198,198```

The second field identifies that this individual (MJ873) has individual MJ872 as a maternal haploid parent.


##Output

Two output files are generated for each input:

`inputfilestem_het.csv `

outputs the percent heterozygosity for each diploid offspring.

`inputfilestem_paternity.csv`

contains the haploid multilocus genotype of all maternal and inferred paternal individuals. This file may be used directly in the `clonecheck.py` script.

The script will print an error if a maternal genotype contains heterozygous loci.
The script will also list all offspring for which the maternal allele cannot be found.

##Example Usage

```
python paternity.py examples/tenellum.csv 13
```

The number indicates the number of loci in the input file.

# `matingpatterns.R`

##Description
This script calculates population genetics statistics in natural populations of moss sporophytes. If multiple sporophytes are sampled that share a maternal haploid parent, they will not be independent draws from the population. This script addresses this bias by generating "pseudopopulations," choosing one sporophyte from each "brood" many times. This generates a range of values for population statistics such as allelic diversity and mating statistics like the inbreeding coefficient FIS.

Imagine a population with three maternal gametophytes (A, B, and C), each of which bears four sporophytes (numbered 1-4). A subset population would contain one sporophyte randomly drawn from each mother, such as: [A1, B2, C3]. We calculate diversity and inbreeding statistics on this subset population. We then repeat this procedure by randomly drawing another sporophyte (with replacement) from each mother and recalculating each statistic on this new subset, such as: [A3, B2, C4]. 

The script will also calculate the linear relationship between percent heterozygosity and sporophyte size (a proxy for fitness) in each pseudopopulation. 

The following statistics are calculated for the pseudopopulations:

* Allelic Diversity (Simpson's diversity index or Effective Number of Alleles, Ne)
* Information Index (Shannon's information index, I)
* Inbreeding coefficient (FIS) = 1 / ExpectedHeterozygosity - ObservedHeterozygosity
* Number of Fathers per brood (K)
* Effective number of fathers per brood (Ke)
* Paternity skew (K/Ke)
* Inbreeding Depression (linear correlation between sporophyte size and % heterozygosity)

The script returns the mean confidence interval (2.5% to 97.5%) for FIS, Ne, I, K, Ke, and K/Ke, as well as the slope, regression coefficient, and p-value for inbreeding depression.
 
##Input
The file must be a tab-delimited text file containing: 

* Sporophytes (their IDs)
* Maternal Shoot ID (grouping)
* Microsatellite loci (pairs of columns)
* Size (Fitness Proxy)
* Mom Clone ID
* Dad Clone ID
* Percent Heterozygosity

An example input file can be found here: `examples/tenellum_2ns.txt`

Indicate the full path to the input file on the first line of the script:

```filename= "path/to/your/file.txt"```

You can also set the number of replicate populations below that.

##Output
A summary of the pseudopopulations is stored in the variable `all_stats`, and is also copied to the clipboard. See the end of the script for the order of the statistics and an explanation of each.



##Example usage
This script is not used from the command line, you must open it in R and change the paramaters noted above, then run the script.


Relevant Publications
========

**clonecheck.py** has been used in the following publications:

* **Johnson, M. G.**, and A. J. Shaw. 2015. Genetic diversity, sexual condition, and microhabitat preference determine mating patterns in *Sphagnum* (Sphagnaceae) peat‐mosses. *Biological Journal of the Linnean Society* 115:96–113.
* Mikulášková, E., M. Hájek, A. Veleba, **M. G. Johnson**, T. Hájek, and J. A. Shaw. 2014. Local adaptations in bryophytes revisited: the genetic structure of the calcium‐tolerant peatmoss *Sphagnum warnstorfii* along geographic and pH gradients. *Ecology and Evolution*.
* Ricca, M., P. Szövényi, E. M. Temsch, **M. G. Johnson**, and A. J. Shaw. 2011. Interploidal hybridization and mating patterns in the *Sphagnum subsecundum* complex. *Molecular Ecology* 20:3202–3218.


**paternity.py** has been used in the following publications:

* **Johnson, M. G.**, and A. J. Shaw. 2015. Genetic diversity, sexual condition, and microhabitat preference determine mating patterns in *Sphagnum* (Sphagnaceae) peat‐mosses. *Biological Journal of the Linnean Society* 115:96–113.

**matingpatterns.R** has been used in the following publications:

* **Johnson, M. G.**, and A. J. Shaw. 2015. Genetic diversity, sexual condition, and microhabitat preference determine mating patterns in *Sphagnum* (Sphagnaceae) peat‐mosses. *Biological Journal of the Linnean Society* 115:96–113.
