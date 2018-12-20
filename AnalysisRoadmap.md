### To-Do
#### In progress
- [ ] What percent of sequences map to the reference
	- Rule written, need to test
- [ ] Finish implementing support for using external databases as the reference for optifit.
    - On branch `silva-ref-db`. Merge into `master` when complete.
    
#### Not in progress
- [ ] Determine how much of the dataset is needed to generate the reference (when using the dataset as its own reference
).
    - Need to fix the sample size and vary the reference size to answer this.
    - Ask Sarah to implement a sample.accnos parameter in OptiFit
- [ ] Replace the `{dataset}.batch` scripts in `code/data_processing/` with snakemake workflow(s).
    - Call the data processing workflows from the main Snakefile with the `include` statement.
- [ ] Download the human sample dataset.


# Analysis Roadmap

## 0. Relevant repositories to get data & see how modeling was done
* De novo (peerj and msystems papers): [(GitHub repo)](https://github.com/SchlossLab/Schloss_Cluster_PeerJ_2015)
	- Includes code for open/closed referencing using VSEARCH
* OptiClust: [(GitHub repo)](https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017)
	- Includes human, murine, marine, and soil datasets
* Baxter: [(GitHub repo)](https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015)
	- Describes RF modeling and CRC

## 1. Independent database
* Questions:
	- How well does fitting into a reference perform vs. de novo
	- How does fitting into a reference using OptiFit perform vs. VSEARCH-based open/closed clustering
	- How important is the composition of the reference database?
 	- How sensitive is OptiFit to random variation?
* First, use SILVA (without greengenes or rdp):
    - [x] Download silva v132
    - [ ] Cluster the silva reference (OptiClust)
    - [ ] For test dataset (soil_1000):
        - [ ] De novo clustering on samples (OptiClust)
        - [ ] Open reference clustering against silva reference
        - [ ] Closed reference clusering against silva reference
* Create reference OTUs and find representative sequence for each OTU using SILVA and greengenes databases
	- [ ] Align greengenes database to SILVA reference alignment (talk to me about what these should be)
	- [ ] Trim to V4 region of 16S rRNA gene
	- [ ] Use OptiClust to cluster reference sequences
* For each dataset (human, murine, marine, and soil), with 10 replicates, perform & compare MCCs:
    - [ ] De novo clustering on  samples
    - [ ] Open reference clustering (OptiFit) against reference OTUs
	- [ ] Closed reference clustering (OptiFit) against reference OTUs
    	- [ ] % of sequences that map to references
	- [ ] Open reference clustering (VSEARCH) against representative sequences from each OTU
	- [ ] Closed reference clustering (VSEARCH) against representative sequences from each OTU
    	- [ ] % of sequences that map to references
	- [ ] Open reference clustering (VSEARCH) against core_gg_97 collection
	- [ ] Closed reference clustering (VSEARCH) against core_gg_97 collection
    	- [ ] % of sequences that map to references

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
	- How does OptiFit on fraction of dataset perform relative to entire dataset?
	- How sensitive is OptiFit to random variation?
* Perform de novo clustering on human, murine, marine, and soil samples (same as above)
	- [x] get MCC values
	- [x] 10 reps
	- datasets:
		- [ ] human
		- [x] murine
		- [x] marine
		- [x] soil
* Repeat the following N times...
	- [x] Use of 10, 20, 30, 40, 50, 60, 70, 80, 90% of data to create reference databases (R)
		- [x] Cluster reference sequences using OptiClust
		- [x] Get MCC
		- [x] Run 10 times to get the seed that gives the best clustering for each set of references
	- [x] On remaining fraction of data (U)
		- [x] OptiClust on its own
			- [x] get MCC values
			- [x] 10 random seeds
		- [x] Use OptiFit/open to fit U sequences into R OTUs
			- [x] Capture MCC values for U sequences with and without R sequences included
			- [x] % of sequences that map to references
			- [x] Run 10 times to see level of variation
		- [x] Use OptiFit/closed to fit U sequences into R OTUs
			- [x] Capture MCC values for U sequences with and without R sequences included
			- [x] % of sequences that map to references
			- [x] Run 10 times to see level of variation

## 4. Deploy
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
