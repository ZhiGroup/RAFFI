
## RAFFI (v.0.1):

RAFFI: accurate and fast familial relationship inference in large scale biobank studies using RaPID

RAFFI is an efficient and robust tool to infer relatedness in large phased genotype data. The efficienct of RAFFI is based on the observation that for very large cohorts, it is much faster to identify all pairwise IBD segments long enough than to identify all pairwise genotype similarity high enough. The kinship coefficients are computed using detected Identical by Descent (IBD) segments from RaPID. Unlike existing IBD-based metbods, RAFFI adopts a data-driven approach that adjusts the estimation based on phasing and genotyping quality. 

![kinship thresholds](https://github.com/ZhiGroup/RAFFI/blob/master/kinship_threshold.png)
|:--:| 
| *RAFFI adjusts the kinship coefficient thresholds to infer the degrees of relatedness. Kinship coefficients computed by the total sum of IBDs using RaPID in a dataset with phasing and genotyping errors with the expected kinship coefficient thresholds (left), and adjusted kinship coefficient thresholds for different degrees of relatedness accounting for phasing/genotyping errors (right)*.|

![simulation_results](https://github.com/ZhiGroup/RAFFI/blob/master/king_raffi_simulation_res-1.png)
|:--:| 
| *Precision, Recall and F1 values (harmonic mean of precision and recall) for calling categories of close relatedness for RAFFI and KING using simulated datasets with different genotyping/phasing errors, maker densities and admixture populations*.|

### Usage:
A binary version of the code is provided in Debug folder (RAFFI_v.0.1). The libraries have been linked statically. If you compile the source code without -static tag (dynamic linking), then you will need to add the path to the compiled boost library to $LD_LIBRARY_PATH:
<br>
`export LD_LIBRARY_PATH=<boost_installation_path>/boost/lib/:$export LD_LIBRARY_PATH`
<br>
RaPID binary (v.1.7) and a Python script to estimate the parameters for RaPID are also provided in the bin directory.

The IBD segments from RaPID are given to the program (using -O tag). If the tag is not used, then the program estimates the parameters for RaPID (by executing the python script) and runs RaPID. The input phased data should be provided in a folder containing compressed (.gz) VCF files. The genetic mapping files should also be provided in a folder containing the genetic location for each site.

Genetic Mapping File Format (tab-delimited):
<br>
`<site_number> <genetic_location>`
<br>
Each line contains the site index and genetic location of a site, in the same order as the VCF input file. Please note that the genetic mapping file should be monotonically increasing, and the genetic location should be provided for each site.


<pre>
Required parameters:
-i [input folder]
        Path to the folder containing gzipped VCF files.
-v [vcf prefix]
        Prefix of gzipped VCF files. e.g. hap1_chr for hap1_chr1.vcf.gz
-g [genetic mapping file path]
        Directory containing genetic maps that were given to RaPID as inputs.
        Map for Chromosome i must be named chr{i}.rMap.
Optional parameters:
-O [output directory of RaPID results]
        Using this tag, RaPID will not be run and the provided outputs will be used directly for relatedness inference.
-o [output directory]
        Output will be written to {output directory}/predictions.txt.
        Default is current direcotry.
-d [max degree]
        Maximum target degree (4 is largest supported degree).
        Default is 4.
-t [number of threads]
        Optimal number of threads is the number of chromosomes.
        Default is 22.
-p [Python version]
        Python path
        Default is python3.6
</pre>

### Output file:
The inferred relationships are stored in a text file along with the computed kinship coefficients and the probability of zero IBD (IBD0), IBD1, and IBD2 between any two pair of related individuals:

ID1   | ID2   |  KINSHIP |IBD0  |  IBD1 |   IBD2   | TYPE
----- | ------|--------|--------|--------|--------|----
7596  |  7642 | 0.1004 | 0.7155 | 0.1674 | 0.1171 | FS
7596  |  9552 | 0.0100 | 0.9598 | 0.0402 | 0.0000 | 4th
7596  |  8643 | 0.0119 | 0.9525 | 0.0475 | 0.0000 | 4th
7596  |  8114 | 0.0103 | 0.9589 | 0.0411 | 0.0000 | 4th


### Installation:
To compile the source code, you will need to install the boost library and modify the boost library path in the Make file. C++14 support is also required to compile the code.

### Citation:
Naseri et al. PLoS Genetics, 2021 accepted.
