# seq_annotation
This aims to annotates the start/end sites of transcripts and start sites of peptides from short-read sequencing data.

### programming
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
### number of different kinds of promoters
In a literiture about E. coli the proportion of internal promoter is 30%, while promoters upstream of annotated genes contribute 17% [1]. It's a high estimation for internal promoter, though.

**reference**
1. [Wade JT. Where to begin? Mapping transcription start sites genome-wide in Escherichia coli. J Bacteriol. 2015;197(1):4-6. Epub 2014/10/22. doi: 10.1128/JB.02410-14] (http://jb.asm.org/content/197/1/4.full)
