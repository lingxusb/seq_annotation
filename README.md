# seq_annotation
This aims to annotates the start/end sites of transcripts and start sites of peptides from short-read sequencing data.

Currently I used a rule-based system to identify Transcription start/end sites from RNAseq data, which includes sites
* near gene coding region in gff file
```python
exp1,sites1 = refine_TSS(BHI_R1f,pos_TIS1)
```
* between gene coding region
```python
exp2,sites2 = forphan(BHI_R1f,pos_TIS1)
```
* inside gene coding region
```python
exp3,sites3 = finternal(BHI_R1f,pos_TIS1)
```
