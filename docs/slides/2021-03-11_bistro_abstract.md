
# OptiFit: a fast method for fitting amplicon sequences to existing OTUs for comparing microbial communities

Amplicon sequencing has become a mainstay of microbial ecology and
host-associated microbiome research.
Researchers can affordably generate millions of sequences to characterize the
composition of hundreds of samples from culture-independent microbial
communities. 
In a typical analysis pipeline, 16S rRNA gene sequences are assigned to
Operational Taxonomic Units (OTUs) to facilitate comparison of taxonomic
composition between communities.
Assigning amplicon sequences to OTUs is an important step in characterizing the
composition of microbial communities across large datasets.
OptiClust, a _de novo_ OTU clustering method in the mothur program, has been
shown to produce higher quality OTU assignments than other methods and at
comparable or faster speeds.
A notable difference between _de novo_ clustering and database-dependent methods
is that OTU assignments clustered with _de novo_ methods are not stable when new
sequences are added to a dataset.
However, in some cases one may wish to incorporate new samples into a previously
clustered dataset without performing clustering again on all sequences, such as
using already-trained machine learning models to make predictions
on new data or comparing OTUs across studies.
To provide an efficient and robust method to fit amplicon sequence data to
existing OTUs, we developed the OptiFit algorithm as a new component of the
mothur program.
We tested the OptiFit algorithm on four datasets isolated from soil, marine,
mouse gut, and human gut microbiome samples and compared OTU quality and 
performance against existing _de novo_ and reference-based algorithms.
OptiFit overcomes the limitations of _de novo_ clustering in OTU stability 
without sacrificing OTU quality, opening up new possibilities for comparing
microbial communities.
