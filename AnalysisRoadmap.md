# Analysis Roadmap

## 0. Relevant repositories to get data & see how modeling was done
* De novo (peerj and msystems papers): [(GitHub repo)](https://github.com/SchlossLab/Schloss_Cluster_PeerJ_2015)
    - Includes code for open/closed referencing using VSEARCH
* OptiClust: [(GitHub repo)](https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017)
    - Includes human, murine, marine, and soil datasets
* Baxter: [(GitHub repo)](https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015)
    - Describes RF modeling and CRC
* [MiSeq SOP](https://mothur.org/wiki/MiSeq_SOP)

## 1. Independent database
* Questions:
    - How well does fitting into a reference perform vs. de novo
    - How does fitting into a reference using OptiFit perform vs. VSEARCH-based open/closed clustering
    - How important is the composition of the reference database?
    - How sensitive is OptiFit to random variation?
* Use SILVA database to generate full length database (like gg_97):
    - [x] Download full length silva v132 (https://www.mothur.org/wiki/Silva_reference_files)
    - [x] Remove archaea and eukarya and select for bacteria (i.e. get.lineage(taxon=Bacteria))
    - [x] Filter sequences for full length sequences (i.e. summary.seqs, screen.seqs, filter.seqs)
    - [x] Cluster the silva reference (OptiClust)
    - [x] Find representative sequence from each OTU (i.e. get.oturep w/ method=abundance)
* Use RDP database to generate full length database (like gg_97):
    - [x] Download rdp (https://www.mothur.org/wiki/RDP_reference_files)
    - [x] Remove archaea and eukarya and select for bacteria
    - [x] Align to silva SEED reference alignment (i.e. align.seqs)
    - [x] Filter sequences for full length sequences
    - [x] Cluster the rdp reference (OptiClust)
    - [x] Find representative sequence from each OTU
* Use greengenes database to generate full length database (like gg_97):
    - [x] Download greengenes (https://mothur.org/wiki/Greengenes-formatted_databases)
    - [x] Remove archaea and eukarya and select for bacteria
    - [x] Realign to silva SEED reference alignment (original gg ref alignment is bad)
    - [x] Filter sequences for full length sequences
    - [x] Cluster the greengenes reference (OptiClust)
    - [x] Find representative sequence from each OTU
* Create reference databases and OTUs for V4 regions
    - [x] Trim SILVA to V4 region of 16S rRNA gene (i.w. pcr.seqs - ask Pat for start/stop coordinates)
    - [x] Trim RDP to V4 region of 16S rRNA gene
    - [x] Trim greengenes to V4 region of 16S rRNA gene
    - [x] Use OptiClust to cluster reference sequences
    - [x] Find representative sequence from each OTU
* For each dataset (human, murine, marine, and soil) and region (V4, full length), with 10 replicates, perform & compare MCCs:
    - [x] OptiClust: De novo clustering on samples
    - [x] OptiFit
        - [x] Open reference clustering
        - [x] Closed reference clustering
    - [ ] VSEARCH / QIIME2 against representative sequences from each OTU
        - [ ] Open
        - [ ] Closed
    - [ ] VSEARCH / QIIME2 against core_gg_97 collection
        - [ ] Open
        - [ ] Closed   
    - [ ] Find % of sequences that map to references for all closed methods.

## 2. Importance of database
* Questions:
    - Is there an optimal way to select sequences to be in the database?
    - Suspect that random without accounting for abundance is just as good as all the rest
* Select 20%(?) of data to be the reference by each of the following methods...
    - [x] Random, without accounting for abundance (10 random seeds)
    - [x] Random, skewed by abundance (10 random seeds)
    - [x] Sequences with the most connections
    - [x] Sequences with the fewest connections
* For each method of generating reference
    - [x] Cluster reference sequences using OptiClust, get MCC value
    - [x] With remaining 80%(?) of data
        - [x] OptiClust sequences, get MCC values
        - [x] OptiFit sequences
            - [x] Capture MCC values for U sequences with and without R sequences included
            - [x] % of sequences that map to references
            - [x] Run 10 times to see level of variation

## 3. Use dataset as its own reference
* Questions:
    - How much of the dataset is needed to generate the reference?
        - Need to fix the sample size and vary the reference size to answer this.
    - How does OptiFit perform on fraction of dataset relative to entire dataset?
    - How sensitive is OptiFit to random variation?
* For each dataset (human, marine, mouse, soil):
    - [ ] Use of 10, 20, 30... 90% of data to create reference databases (R)
        - [ ] Subset with different weighted methods:
            - none
            - abundance
            - distances (connectivity)
        - [ ] Cluster reference sequences using OptiClust
        - [ ] Get MCC
        - [ ] Run 10 times to get the seed that gives the best clustering for each set of references
    - [ ] On remaining fraction of data (U)
        - [ ] OptiClust on its own
        - [ ] Use OptiFit/open to fit U sequences into R OTUs
            - [ ] Capture MCC values for U sequences with and without R sequences included
            - [ ] % of sequences that map to references
            - [ ] Run 10 times to see level of variation
        - [ ] Use OptiFit/closed to fit U sequences into R OTUs
            - [ ] Capture MCC values for U sequences with and without R sequences included
            - [ ] % of sequences that map to references
            - [ ] Run 10 times to see level of variation

## 4. Deploy (separate paper)
* Questions:
    - Can we use OptiFit to fit new data to existing data and model to predict presence of lesion?
* Use Baxter dataset
    - [ ] Cluster all samples (N~490) using OptiClust
    - [ ] Train RF model to differentiate normals from lesions
    - [ ] Get AUC, sens/spec
* Generate N references that each lack one sample
    - [ ] Use OptiClust to generate N reference clusterings
    - [ ] Train RF model to differentiate normals from lesions
    - [ ] OptiFit sequences from missing sample to respective reference
    - [ ] Classify sample against model
    - [ ] Get AUC, sens/spec
