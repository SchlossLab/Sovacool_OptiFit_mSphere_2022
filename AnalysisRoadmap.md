0. Relevant repositories to get data, see how modeling was done
* de novo (peerj and msystems papers): https://github.com/SchlossLab/Schloss_Cluster_PeerJ_2015
	- Includes code for open/closed referencing using VSEARCH
* OptiClust: https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017
	- Includes human, murine, marine, and soil datasets
* Baxter: https://github.com/SchlossLab/Baxter_glne007Modeling_GenomeMed_2015
	- Describes RF modeling and CRC

1. Independent database.
* Questions...
	- How well does fitting into a reference perform vs. de novo
	- How does fitting into a reference using OptiFit perform vs. VSEARCH-based open/closed clustering
	- How important is the composition of the reference database?
 	- How sensitive is OptiFit to random variation?
* Create reference OTUs and find representative sequence for each OTU using SILVA and greengenes databases
	- Align greengenes database to SILVA reference alignment (talk to me about what these should be)
	- Trim to V4 region of 16S rRNA gene
	- Use OptiClust to cluster reference sequences
* Perform de novo clustering on human, murine, marine, and soil samples
	- get MCC values
	- 10 reps
* Perform open reference clustering (OptiFit) against reference OTUs
	- get MCC values
	- 10 reps
* Perform closed reference clustering (OptiFit) against reference OTUs
	- % of sequences that map to references
	- get MCC values
	- 10 reps
* Perform open reference clustering (VSEARCH) against representative sequences from each OTU
	- get MCC values
	- 10 reps
* Perform closed reference clustering (VSEARCH) against representative sequences from each OTU
	- % of sequences that map to references
	- get MCC values
	- 10 reps
* Perform open reference clustering (VSEARCH) against core_gg_97 collection
	- get MCC values
	- 10 reps, randomize reference order each time
* Perform closed reference clustering (VSEARCH) against core_gg_97 collection
	- % of sequences that map to references
	- get MCC values
	- 10 reps, randomize reference order each time


2. Importance of database
* Questions...
	- Is there an optimal way to select sequences to be in the database?
	- Suspect that random without accounting for abundance is just as good as all the rest
* Select 20%(?) of data to be the reference by each of the following methods...
	- Random, without accounting for abundance (10 random seeds)
	- Random, skewed by abundance (10 random seeds)
	- Sequences with the most connections
	- Sequences with the fewest connections
* For each method of generating reference
	- Cluster reference sequences using OptiClust, get MCC value
	- With remaining 80%(?) of data
		- OptiClust sequences, get MCC values
		- OptiFit sequences
 			- Capture MCC values for U sequences with and without R sequences included
			- % of sequences that map to references
			- Run 10 times to see level of variation


3. Use dataset as its own reference
* Questions...
	- How much of the dataset is needed to generate the reference?
	- How does OptiFit on fraction of dataset perform relative to entire dataset?
	- How sensitive is OptiFit to random variation?
* Perform de novo clustering on human, murine, marine, and soil samples (same as above)
	- get MCC values
	- 10 reps
* Repeat the following N times...
	- Use of 10, 20, 30, 40, 50, 60, 70, 80, 90% of data to create reference databases (R)
		- Cluster reference sequences using OptiClust
		- Get MCC
		- Run 10 times to get the seed that gives the best clustering for each set of references
	- On remaining fraction of data (U)
		- OptiClust on its own
			- get MCC values
			- 10 random seeds
		- Use OptiFit/open to fit U sequences into R OTUs
			- Capture MCC values for U sequences with and without R sequences included
			- Run 10 times to see level of variation
		- Use OptiFit/closed to fit U sequences into R OTUs
			- Capture MCC values for U sequences with and without R sequences included
			- % of sequences that map to references
			- Run 10 times to see level of variation


4. Deploy
* Questions...
	- Can we use OptiFit to fit new data to existing data and model to predict presence of lesion?
* Use Baxter dataset
	- Cluster all samples (N~490) using OptiClust
	- Train RF model to differentiate normals from lesions
	- Get AUC, sens/spec
* Generate N references that each lack one sample
	- Use OptiClust to generate N reference clusterings
	- Train RF model to differentiate normals from lesions
	- OptiFit sequences from missing sample to respective reference
	- Classify sample against model
	- Get AUC, sens/spec
