myvatn_charr
========

Data and code for reproducing the main analyses in the manuscript “Opposing trends in survival and recruitment slow recovery of an historically overexploited fishery.”
-------

Joe Phillips, Hólar University and University of Wisconsin-Madison

## Contents

* `model`: Files for fitting the model in Stan and storing output. The files are defined as follows:
	* site_data.csv: yearly survey catch for each site and age class; primary data for analysis. The file contains the following columns:
		- year: sample year
		-stage: age class (“first”, “second”, “third”, “adult”)
		-area: numerical code for sampling site
		-count: number of sampled individuals

* `analyses`: Analyses of the model fits, including code for figures.
