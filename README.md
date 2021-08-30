# RfxCas13d guide RNA scoring and selection

## Installation and Dependencies

To Download & Install:

```bash
git clone https://github.com/chewlabSB2/cas13gRNAtor.git
cd cas13gRNAtor 
python setup.py install
```
Only available on Linux OS. setup.py will automatically download the required R packages and python requirements

## Cas13gRNAtor Simple Usage
With reference to the test dataset (EV71 Virus)

```bash
mafft --add NonMSAEV71.fasta --keeplength --reorder --nomemsave --thread $thread referenceEV71.fasta 1> MSAEV71.fasta 2> error.log
cas13gRNAtor -r referenceEV71.fasta -m MSAEV71.fasta -bc NonMSAEV71.fasta -t 8 -p cas13gRNAtor --mismatch 5 -bo Homo_sapiens.GRCh38.cdna.all.fa
```
cas13gRNAtor includes [Wessel's gRNA scoring](https://gitlab.com/sanjanalab/cas13/-/tree/master/Cas13designGuidePredictor) and [bowtie](http://bowtie-bio.sourceforge.net/manual.shtml) to filter for offtargets with up to 3 mismatches. Only 3rd and 4th quartile from Wessel's gRNA scoring and prediction are included in the final output. Mafft's output is used to generate a consensus sequence to calculate the Shannon Entropy and conservation scores for each guide RNA. cas13gRNAtor prints the top 100 guide RNA based on a weighted score. 

## Cas13gRNAvalidate Simple Usage

```bash
cas13gRNAvalidate -i investigatedEV71.txt -m MSAEV71.fasta -bc NonMSAEV71.fasta -t 8 -p cas13gRNAvalidate-test --mismatch 5 -r referenceEV71.fasta -bo Homo_sapiens.GRCh38.cdna.all.fa 
```

cas13gRNAvalidate evaluates the conservation and entropy scores of the gRNAs previously designed for further inspection at a later time as more viral genomes (such as SARS-CoV-2) are being sequenced. cas13gRNAvalidate will print the conservation and entropy score for each gRNA. The script can find offtargets as well with the -bo flag. 

## Expected Output

1. bestgRNAs.csv: Includes the top (default: 100) candidate gRNAs based on the weighted scoring schemer
2. combined_results.csv: Includes all candidate gRNAs sequence, offtarget Information and entropy/ conservation/ guide (from Wessels)/ weighted Scores
3. alternatives.csv: Include candidate gRNAs sequence with multiple target on reference
4. unfound.csv: Include candidate gRNAs sequence which are not found in the consensus sequence
5. entropy_score.csv and conservation_score.csv: Includes the individual entropy and conservation score of every nucleotide in all candidate gRNAs as supplementary data
6. position_score.csv: Includes the nucleotide, individual entropy and conservation score of every nucleotide of the consensus sequence generated from the multiple sequence alignment as supplementary data 
7. plot the conservation/ entropy and gRNA occurrence of the entire consensus sequence
8. plot the top (default: 100) candidate gRNAs and highlight the intolerant region for cas13 targetting activities in a sub-directory called best_gRNA

Add the flag --keep_tmp to retain the temporary files 

## Possible Issues

If you face any issue with libgsl.so.0, refer to [this page](https://stackoverflow.com/questions/22222666/error-while-loading-shared-libraries-libgsl-so-0-cannot-open-shared-object-fil). 

All other dependencies are included in the misc & data sub-directory