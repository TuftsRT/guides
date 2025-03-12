#######
 Homer
#######

**************
 Introduction
**************

HOMER is a suite of tools for Motif Discovery and next-gen sequencing
analysis.

**********
 Versions
**********

-  4.11

**********
 Commands
**********

-  addDataHeader.pl
-  addData.pl
-  addGeneAnnotation.pl
-  addInternalData.pl
-  addOligos.pl
-  adjustPeakFile.pl
-  adjustRedunGroupFile.pl
-  analyzeChIP-Seq.pl
-  analyzeHiC
-  analyzeRepeats.pl
-  analyzeRNA.pl
-  annotateInteractions.pl
-  annotatePeaks.pl
-  annotateRelativePosition.pl
-  annotateTranscripts.pl
-  assignGeneWeights.pl
-  assignGenomeAnnotation
-  assignTSStoGene.pl
-  batchAnnotatePeaksHistogram.pl
-  batchFindMotifsGenome.pl
-  batchFindMotifs.pl
-  batchMakeHiCMatrix.pl
-  batchMakeMultiWigHub.pl
-  batchMakeTagDirectory.pl
-  batchParallel.pl
-  bed2DtoUCSCbed.pl
-  bed2pos.pl
-  bed2tag.pl
-  blat2gtf.pl
-  bridgeResult2Cytoscape.pl
-  changeNewLine.pl
-  checkPeakFile.pl
-  checkTagBias.pl
-  chopify.pl
-  chopUpBackground.pl
-  chopUpPeakFile.pl
-  cleanUpPeakFile.pl
-  cleanUpSequences.pl
-  cluster2bedgraph.pl
-  cluster2bed.pl
-  combineGO.pl
-  combineHubs.pl
-  compareMotifs.pl
-  condenseBedGraph.pl
-  configureHomer.pl
-  cons2fasta.pl
-  conservationAverage.pl
-  conservationPerLocus.pl
-  convertCoordinates.pl
-  convertIDs.pl
-  convertOrganismID.pl
-  createIGVhtml.pl
-  duplicateCol.pl
-  eland2tags.pl
-  fasta2tab.pl
-  fastq2fasta.pl
-  filterListBy.pl
-  filterTADsAndCPs.pl
-  filterTADsAndLoops.pl
-  filterTagDirectory.pl
-  findcsRNATSS.pl
-  findGO.pl
-  findGOtxt.pl
-  findHiCCompartments.pl
-  findHiCDomains.pl
-  findHiCInteractionsByChr.pl
-  findKnownMotifs.pl
-  findMotifsGenome.pl
-  findMotifs.pl
-  findPeaks
-  findRedundantBLAT.pl
-  findTADsAndLoopsFromRelMatrix
-  findTADsAndLoops.pl
-  findTopMotifs.pl
-  flipPC1toMatch.pl
-  freq2group.pl
-  genericConvertIDs.pl
-  genomeOntology
-  GenomeOntology.pl
-  getChrLengths.pl
-  getConservedRegions.pl
-  getDifferentialBedGraph.pl
-  getDifferentialPeaks
-  getDifferentialPeaksReplicates.pl
-  getDiffExpression.pl
-  getDistalPeaks.pl
-  getFocalPeaks.pl
-  getGenesInCategory.pl
-  getGenomeTilingPeaks
-  getGWASoverlap.pl
-  getHiCcorrDiff.pl
-  getHomerQCstats.pl
-  getLikelyAdapters.pl
-  getMappableRegions
-  getMappingStats.pl
-  getPartOfPromoter.pl
-  getPeakTags
-  getPos.pl
-  getRandomReads.pl
-  getSiteConservation.pl
-  getTopPeaks.pl
-  gff2pos.pl
-  go2cytoscape.pl
-  groupSequences.pl
-  homer
-  homer2
-  HomerConfig.pm
-  HomerSVGLogo.pm
-  homerTools
-  joinFiles.pl
-  loadGenome.pl
-  loadPromoters.pl
-  makeBigBedMotifTrack.pl
-  makeBigWig.pl
-  makeBinaryFile.pl
-  makeHiCWashUfile.pl
-  makeMetaGeneProfile.pl
-  makeMultiWigHub.pl
-  makeTagDirectory
-  makeUCSCfile
-  map-fastq.pl
-  merge2Dbed.pl
-  mergeData.pl
-  mergePeaks
-  motif2Jaspar.pl
-  motif2Logo.pl
-  old
-  parseGTF.pl
-  pos2bed.pl
-  preparseGenome.pl
-  prepForR.pl
-  profile2seq.pl
-  qseq2fastq.pl
-  randomizeGroupFile.pl
-  randomizeMotifs.pl
-  randRemoveBackground.pl
-  removeAccVersion.pl
-  removeBadSeq.pl
-  removeOutOfBoundsReads.pl
-  removePoorSeq.pl
-  removeRedundantPeaks.pl
-  renamePeaks.pl
-  resizePosFile.pl
-  revoppMotif.pl
-  rotateHiCmatrix.pl
-  runHiCpca.pl
-  sam2spliceJunc.pl
-  scanMotifGenomeWide.pl
-  scrambleFasta.pl
-  SIMA.pl
-  Statistics.pm

*************
 Example job
*************

Adjust slurm options based on job requirements (`slurm cheat sheet
<https://slurm.schedmd.com/pdfs/summary.pdf>`_):

.. code::

   #!/bin/bash
   #SBATCH -p partitionName  # batch, gpu, preempt, mpi or your group's own partition
   #SBATCH -t 1:00:00  # Runtime limit (D-HH:MM:SS)
   #SBATCH -N 1   # Number of nodes
   #SBATCH -n 1   # Number of tasks per node
   #SBATCH -c 4   # Number of CPU cores per task
   #SBATCH --mem=8G       # Memory required per node
   #SBATCH --job-name=homer       # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load homer/XXXX ### Latest version is recommended.
