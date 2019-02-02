### The CLI works fairly simply. 

To dump a sequence:

* nodefile - Absolute Path to the nodefile you wish to load
* edgefile - Aboslute path to the edgefile you wish to load
* chromosome - Number of the chromosome you are currently evaluating
* sequenceName - The name of the strain you wish to retrieve the path for

`python CCGG_CLI.py nodefile edgefile chromosome [--sequence strain path filedest]`

Example:
In the following example I extract the sequence for strain A from chromosome 14 following the 0 path (because A is not heterozygous so there is only one path) and dump it to the file named "sequenceAChr14.txt"

`python CCGG_CLI.py NodeFile_Chr14.fa EdgeFile_Chr14.fa 14 --sequence A 0 sequenceAChr14.txt`