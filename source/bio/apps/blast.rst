#######
 Blast
#######

**************
 Introduction
**************

BLAST (Basic Local Alignment Search Tool) finds regions of similarity
between biological sequences. The program compares nucleotide or protein
sequences to sequence databases and calculates the statistical
significance.

***********
 Databases
***********

Local copies of blast databases can be found in the directory
`/cluster/tufts/biocontainers/datasets/blast/latest/`.

In blast modulefiles, we defined the environment variable `BLASTDB` as
`/cluster/tufts/biocontainers/datasets/blast/latest`, so that you can
use the blast databases such as `nr` and `nt` without specifying the
full path.

**********
 Versions
**********

-  2.15.0
-  2.16.0

**********
 Commands
**********

-  amino-acid-composition
-  between-two-genes
-  blastdb_aliastool
-  blastdbcheck
-  blastdbcmd
-  blast_formatter
-  blastn
-  blastp
-  blastx
-  cleanup-blastdb-volumes.py
-  deltablast
-  dustmasker
-  eaddress
-  eblast
-  get_species_taxids.sh
-  legacy_blast.pl
-  makeblastdb
-  makembindex
-  makeprofiledb
-  psiblast
-  rpsblast
-  rpstblastn
-  run-ncbi-converter
-  segmasker
-  tblastn
-  tblastx
-  update_blastdb.pl
-  windowmasker

*******************************
 Memory requirements for blast
*******************************

Ensure that you allocate sufficient memory for the BLAST database being
used. Insufficient memory will force the database to be repeatedly read
from disk, leading to intensive I/O operations that significantly slow
down your job. This can also cause severe strain on the filesystem,
especially when running multiple jobs concurrently. Such issues may
disrupt other users and prompt intervention from the Research Technology
(RT) staff.

Calculate the memory required by totalling the size of the .nsq (for
nucleotide databases) or .psq (for protein databases) files. e.g.:

.. code::

   $ du -sh --total /fdb/blastdb/nt*.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.000.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.001.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.002.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.003.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.004.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.005.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.006.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.007.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.008.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.009.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.010.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.011.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.012.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.013.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.014.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.015.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.016.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.017.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.018.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.019.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.020.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.021.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.022.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.023.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.024.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.025.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.026.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.027.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.028.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.029.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.030.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.031.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.032.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.033.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.034.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.035.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.036.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.037.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.038.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.039.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.040.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.041.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.042.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.043.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.044.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.045.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.046.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.047.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.048.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.049.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.050.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.051.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.052.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.053.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.054.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.055.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.056.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.057.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.058.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.059.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.060.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.061.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.062.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.063.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.064.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.065.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.066.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.067.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.068.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.069.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.070.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.071.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.072.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.073.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.074.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.075.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.076.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.077.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.078.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.079.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.080.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.081.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.082.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.083.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.084.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.085.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.086.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.087.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.088.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.089.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.090.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.091.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.092.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.093.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.094.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.095.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.096.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.097.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.098.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.099.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.100.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.101.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.102.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.103.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.104.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.105.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.106.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.107.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.108.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.109.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.110.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.111.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.112.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.113.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.114.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.115.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.116.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.117.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.118.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.119.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.120.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.121.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.122.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.123.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.124.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.125.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.126.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.127.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.128.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.129.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.130.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.131.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.132.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.133.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.134.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.135.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.136.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.137.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.138.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.139.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.140.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.141.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.142.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.143.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.144.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.145.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.146.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.147.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.148.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.149.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.150.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.151.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.152.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.153.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.154.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.155.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.156.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.157.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.158.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.159.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.160.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.161.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.162.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.163.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.164.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.165.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.166.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.167.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.168.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.169.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.170.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.171.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.172.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.173.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.174.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.175.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.176.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.177.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.178.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.179.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.180.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.181.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.182.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.183.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.184.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.185.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.186.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.187.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.188.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.189.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.190.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.191.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.192.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.193.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.194.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.195.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.196.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.197.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.198.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.199.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.200.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.201.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.202.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.203.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.204.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.205.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.206.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.207.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.208.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.209.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.210.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.211.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.212.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.213.nsq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nt.214.nsq
   760M        /cluster/tufts/biocontainers/datasets/blast/latest/nt.215.nsq
   601G        total

   $ du -sh --total /cluster/tufts/biocontainers/datasets/blast/latest/nr*.psq
   3.4G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.000.psq
   3.5G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.001.psq
   3.2G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.002.psq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.003.psq
   3.2G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.004.psq
   3.1G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.005.psq
   3.0G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.006.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.007.psq
   47M     /cluster/tufts/biocontainers/datasets/blast/latest/nr.008.psq
   3.2G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.009.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.010.psq
   649M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.011.psq
   3.2G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.012.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.013.psq
   556M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.014.psq
   3.3G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.015.psq
   3.4G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.016.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.017.psq
   756M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.018.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.019.psq
   126M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.020.psq
   3.6G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.021.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.022.psq
   517M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.023.psq
   3.1G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.024.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.025.psq
   320M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.026.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.027.psq
   970M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.028.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.029.psq
   563M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.030.psq
   2.6G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.031.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.032.psq
   308M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.033.psq
   3.6G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.034.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.035.psq
   484M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.036.psq
   3.5G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.037.psq
   3.0G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.038.psq
   3.0G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.039.psq
   3.0G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.040.psq
   3.3G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.041.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.042.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.043.psq
   103M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.044.psq
   2.7G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.045.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.046.psq
   120M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.047.psq
   3.7G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.048.psq
   2.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.049.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.050.psq
   62M     /cluster/tufts/biocontainers/datasets/blast/latest/nr.051.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.052.psq
   113M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.053.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.054.psq
   212M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.055.psq
   2.7G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.056.psq
   3.4G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.057.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.058.psq
   384M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.059.psq
   3.7G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.060.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.061.psq
   96M     /cluster/tufts/biocontainers/datasets/blast/latest/nr.062.psq
   3.7G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.063.psq
   3.5G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.064.psq
   3.6G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.065.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.066.psq
   140M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.067.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.068.psq
   140M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.069.psq
   3.7G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.070.psq
   3.1G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.071.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.072.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.073.psq
   51M     /cluster/tufts/biocontainers/datasets/blast/latest/nr.074.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.075.psq
   740M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.076.psq
   3.3G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.077.psq
   3.2G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.078.psq
   3.4G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.079.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.080.psq
   385M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.081.psq
   3.7G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.082.psq
   3.4G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.083.psq
   3.1G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.084.psq
   3.3G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.085.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.086.psq
   269M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.087.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.088.psq
   13M     /cluster/tufts/biocontainers/datasets/blast/latest/nr.089.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.090.psq
   106M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.091.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.092.psq
   390M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.093.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.094.psq
   513M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.095.psq
   2.9G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.096.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.097.psq
   895M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.098.psq
   3.2G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.099.psq
   3.3G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.100.psq
   3.3G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.101.psq
   2.5G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.102.psq
   2.9G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.103.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.104.psq
   136M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.105.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.106.psq
   3.2G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.107.psq
   3.5G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.108.psq
   3.5G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.109.psq
   3.6G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.110.psq
   3.1G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.111.psq
   3.7G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.112.psq
   3.1G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.113.psq
   3.8G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.114.psq
   702M        /cluster/tufts/biocontainers/datasets/blast/latest/nr.115.psq
   1.3G        /cluster/tufts/biocontainers/datasets/blast/latest/nr.116.psq
   302G        total

As you can see, if you want to run **blastn** with the **nt** database,
you need to allocate at least **601GB** of memory, or better yet,
**610GB** for safety.

If you want to run **blastp** with the **nr** database, you need to
allocate at least **302GB** of memory, or better yet, **310GB** for
safety.

*************
 Example job
*************

If you want to run a blast job with the **nr** or **nr** database, it is
recommended to use the **largemem** partition.

Adjust slurm options based on job requirements (`slurm cheat sheet
<https://slurm.schedmd.com/pdfs/summary.pdf>`_):

.. code::

   #!/bin/bash
   #SBATCH -p largemem  # batch, gpu, preempt, mpi or your group's own partition
   #SBATCH -t 1:00:00  # Runtime limit (D-HH:MM:SS)
   #SBATCH -N 1   # Number of nodes
   #SBATCH -n 1   # Number of tasks per node
   #SBATCH -c 24  # Number of CPU cores per task
   #SBATCH --mem=400G     # Memory required per node
   #SBATCH --job-name=blast       # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load blast/XXXX ### Latest version is recommended.

   blastp -query query.fasta -db nr -out blastp_results.txt -evalue 1e-5 -outfmt 6 -num_threads 24
