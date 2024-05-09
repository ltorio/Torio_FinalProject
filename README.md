# PHYLOGENY OF UCE HABRONATTUS DATA

# Background

Ultra-Conserved Elements (UCE) are regions of the genome that are conserved across many species. Probe sets are designed to target UCEs of a specific clade for sequencing. This UCE data or more specifically the flanking regions around the UCE locus are used for phylogenomic analysis. 
My UCE data was provided to me by Rodrigo Monjaraz-Ruedas in Dr. Hedin’s Lab. The data was collected using a probe set for the retrolateral tibial apophysis clade of spiders developed by Zhang et al that targets exonic single copy orthologs.

  Zhang, J., Li, Z., Lai, J., Zhang, Z. and Zhang, F. (2023), A novel probe set for the phylogenomics and evolution of RTA spiders. Cladistics, 39: 116-128. https://doi.org/10.1111/cla.12523

The UCE data was collected from 67 spiders including 3 species of Habronattus jumping spiders as well as one potential novel Habronattus species. Each of the 3343 unaligned FASTA files corresponds to a UCE locus. Each file contains up to 67 sequences that correspond to the given UCE locus. Additionally there is a configuration file Hamicus_taxon-set.conf, containing the names of all the spiders that were sampled, that was used by Rodrigo some preliminary data manipulation and that will be using again.
These overall steps are based on Borowiec’s paper that analyzed UCE data on army ants, including the custom python script uce_to_protein.py

  Marek L Borowiec, Convergent Evolution of the Army Ant Syndrome and Congruence in Big-Data Phylogenetics, Systematic Biology, Volume 68, Issue 4, July 2019, Pages 642–656, https://doi.org/10.1093/sysbio/syy088

Also included is a reference protein sequence for Trichonephila antipodiana (inside the Hamicus_unaligned directory).

  Fan Z; Yuan T; Liu P; Wang LY; Jin JF; Zhang F; Zhang ZS (2021): A chromosome‐level genome of the spider Trichonephila antipodiana GigaScience Database. https://doi.org/10.5524/100868

The goal of this project is to develop a phylogeny to see how the novel species is related to the three known Habronattus species. 

# environments and packages

You will need to create two separate conda environments. The testEnvPhylogny environment will contain PHYLUCE 1.6.8 , ASTER 1.16 , and ETE Toolkit 3.1.2 . The phylogenyRaxml environment will contain RAxML 8.2.13 , Python3 or newer, SQLite 3.45, Biopython 1.83, BLAST 2.12

These packages were installed thru either the bioconda or conda forge channels. These channels may be enabled with the code in the 1st code block below. 

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
```
conda create -n testEnvPhylogny
conda activate testEnvPhylogny
conda install phyluce
conda install ASTER
conda install ete3
conda deactivate
conda create -n phylogenyRaxml
conda activate phylogenyRaxml
conda install raxml
conda install python
conda install sqlite
conda install blast ############### conda install blast=2.12.
```

# UCE to Protein

First all the sequence names need to be cleaned to remove UCE locus information. This is done to properly format the sequences for the custom python script.

```
sed -i 's/uce-[0-9]*_//g' *.unaligned.fasta
sed -i 's/ |uce-[0-9]*//g' *.unaligned.fasta
```

Then we need to create a blast database from the reference sequence.
```
cd Hamicus_unaligned/
../uce_to_protein.py blastdb -i Trichonephila_antipodiana.pep.fasta
```
Next we blast each UCE locus against the database. This returns .xml files that contain blast results for each UCE locus and a configuration file that matches a FASTA file to its corresponding .xml file.
```
../uce_to_protein.py queryblast -i *unaligned.fasta
```

From the previous blast results stored in .xml file, we retrieve the best hits. A single SQLite database is created from each .xml file storing the locus name, taxon name, name of annotated protein matched in the reference, nucleotide query trimmed from introns, trimmed and untrimmed protein queries, and protein subject that was matches the reference set.
```
../uce_to_protein.py parse -i fasta_to_xml.conf -o my_uce_hits.sqlite
```

Finally we convert the SQLite database information about the best hits into a more typical FASTA format, which requires a taxon configuration file, Hamicus_taxon-set.conf, that simply lists all the sample names.
```
../uce_to_protein.py queryprot -d my_uce_hits.sqlite -c Hamicus_taxon-set.conf -g all
```

Three different kinds of sequences are produced from these FASTA files: 1) trimmed protein sequence found, 2) its corresponding nucleotide sequence, and 3) the nucleotide sequence with 3rd codon position removed. Future steps will proceed with the nucleotide sequences, so we will move the nucleotide sequence to their own directory.
```
mkdir ../nucleotideHits
cp all-coding-nuc-uce-* ../nucleotideHits/
cd ../nucleotideHits/
```

# Alignment and Cleaning

Before continuing we will filter out loci with fewer than fewer than 70% of the data (fewer taxa 47 taxa). We use a phyluce command phyluce_assembly_get_fasta_lengths that gets taxa count data for each loci, among other information.

```
conda deactivate
conda activate testEnvPhylogny

for i in *.fasta; do phyluce_assembly_get_fasta_lengths --input $i --csv >> nucInfo.csv; done # trying to write output data to a csv file for easy of use
mkdir lowTaxaLoci
awk 'BEGIN { FS = "," } ; { if ($2<=47) {print $1, "has low taxa"; system("mv " $1 " lowTaxaLoci/" ) }}' nucInfo.csv 
```

Now we align sequences at the same locus using MAFFT which conda should have installed with other packages.
```
mkdir ../alignedHits
for i in *.fasta; do mafft $i > ../alignedHits/$i; done
ls ../alignedHits/ | wc -l
ls *fasta | wc -l # double check the numbers are the same for each file count

cd ../alignedHits
rename 's/.unaligned/.aligned/' *.unaligned.fasta
ls | head

```

Gap regions are removed using Gblocks 0.91b which should have also been installed with previous packages.

```
for i in *.fasta; do Gblocks $i b1=0.5 b2=0.5 b3=12 b4=7; done
```
# Phylogeny

Gblocks also produced gb.reduced files which are not properly formatted for RAxML and are deleted. Now the files are ready to be use to create a phylogeny for each of the loci using RAxML (50 bootstraps were used due to computational limitations). 
```
conda deactivate
conda activate phylogenyRaxml
mkdir ../../trees50/
rm *gb.reduced

for i in *.fasta-gb; do raxmlHPC -m GTRGAMMA -f a -# 50 -p 12345 -x 12345 -s $i -n raxml-$i; done

mkdir ../../trees5
mv RAxML_* ../../trees50/
conda deactivate
```

A summary tree is made from the best trees produced by RAxML using ASTER.
```
conda activate testEnvPhylogny
ls RAxML_bestTree.raxml-all-coding-nuc-uce-* > asterInputTest.txt

for i in RAxML_bestTree.raxml*; do cat $i >> asteralInput.txt; done

astral -t 6 -o asteral_output.txt -i asteralInput.txt 2>asteral.log

cp asteral_output.txt ..
cd ..
```

Finally the summary tree is visualized. Once the below command is run you must make sure that branch info, branch support and show leaf names are enabled and that G3077_H_sp_nov leaf selected to produce the same final tree.
```
ete3 view -t asteral_output.txt 
```
All code presented here is also available in the collatedCode.txt file.
