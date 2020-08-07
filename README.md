
## RAFFI (v.0.1):
### Usage:
A binary version of the code is provided in the Debug folder. The libraries have been linked statically. If you compile the source code without -static tag (dynamic linking), then you will need to add the path to the compiled boost library to $LD_LIBRARY_PATH$: export LD_LIBRARY_PATH=<boost_installation_path>/boost/lib/ 
RaPID binary (v.1.7) and a Python script to estimate the parameters for RaPID are also provided in the bin directory.

The IBD segments from RaPID are given to the program (using -O tag). If the tag is not used, then the program estimates the parameters  for RaPID (by executing the python script) and run RaPID. The input phased data should be provided in a folder containing compressed (.gz) VCF files. The genetic mapping files should also be provided in a folder containing the genetic location for each site. 

Genetic Mapping File Format (tab-delimited):
`<site_number> <genetic_location>`
Each line contains the site index and genetic location of a site, the same order as the VCF input file. Please note that the genetic mapping file should be monotically increasing and the genetic locaion should be provided for each site. We have provided two Python scripts to filter the genetic mapping file and also a python script to interpolate the genetic locaions. 

### Compilation:
In order to compile the source code you will need to install boost library and modify the boost library path in the Make file. C++14 support is also required to compile the code.
