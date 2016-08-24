imPhy
=

Welcome! imPhy is a small [Python][python] and [R][r] pipeline for simulating, imputing, and analyzing phylogenetic trees. The main purpose of imPhy is to assess the effectiveness of methods imputing the positions of missing leaves in bifurcating trees. Effectiveness is measured by the distance between the original simulated trees and the imputed trees, using Robinson-Foulds and Billera-Holmes-Vogtmann (BHV) distances. Unless only a small number of trees are being imputed, it is recommended to run imPhy on a server.

imPhy has many moving parts and will be developed very slowly, if at all, so this repository is intended as a proof of concept and record for other researchers who would like to implement a similar pipeline.

---

Quick Start Guide
-
Each file uses docopt, so running the file with `-h` will bring up a help guide if you need it.
 1. Clone this repository. 
 2. Install [dependencies](#Dependencies).
 3. Switch directories, using `cd imPhy/src/`. It is safest to run the files from this directory.
 4. Edit `imPhy/src/cpp/Makefile` to match your Gurobi installation and C++ compiler, and run `make`.
 5. Run a test experiment, with `python3 experiment.py -tc test`. This will fail if you are missing dependencies. If it works, your data will be in `imPhy/test/`.
 6. Congratulations! You can now run your own experiments by editing the "experimental setup" section in `experiment.py`. A good starting command is `python3 experiment.py -pc my_experiment`.


Features
-
There are 5 main steps in the imPhy pipeline:

1. Tree Generation
   - Creates trees with [DendroPy][dp], a very useful phylogenetic library for Python. ImPhy uses the Yule process to create species trees, and the contained coalescent model to create gene trees. Trees are written to `imPhy/my_experiment/batch_A/nexus/`.
   - Parameters:
	 - Species Depth: Height in generations of the species tree. Calculated by `(Effective Population Size) * (C-ratio)`.
	 - Effective Population Size is equivalent to N, for a population of N haploid individuals, or 2N for a population of N diploid individuals. Defaults to 10000.
	 - C-ratio: Ratio of Species Depth to Effective Population Size. Used as a proxy for species depth.
	 - Number of Species Trees: The number of species trees to create with each depth and population size. This is effectively a "number of trials" adjuster.
	 - Number of Species: Number of leaves in the species tree.
	 - Number of Gene Trees: This many gene trees will be coalesced within each species tree. More gene trees means more information for the imputation software, but higher memory requirements.
	 - Number of Individuals per Species: Controls the number of leaves in each gene tree, which is equal to `(Num Individuals per Species) * (Num Species)`
2. Dropping Leaves
   - To impute leaves (individuals), some must be missing. For this purpose imPhy uses the [APE][ape] package in R. There are two methods for choosing the number of leaves to drop, which are chosen automatically based on the value of the leaf dropping parameter p. The identities of leaves to be dropped are chosen randomly regardless of which method is chosen. In `experiment.py`, p is set using the `probs` list. Trees are written to `imPhy/my_experiment/batch_A/data/` in the data format found in `imPhy/src/cpp/examples/`, which uses vectorized distance matrices.
	 - p < 1 causes p to be interpreted as the probability of success (dropping a leaf) in a binomial distribution with size equal to the number of leaves in the tree. A single draw is made from the distribution, which serves as the number of leaves to drop from the tree.
	 - p >= 1 will result in 1/p of the leaves in a given tree being dropped. If the resulting number is not an integer, it will be rounded.
3. Imputation
   - Imputation is carried out using mutual information from multiple gene trees that have coalesced on the same species tree. [Gurobi][gurobi] and [C++][cpp] are used to impute leaves, so it is necessary to have a valid Gurobi installation and to make the C++ file on the machine it will be run on. The imputation software should be easy to change without damaging the rest of the pipeline, provided its inputs and outputs are maintained.
	 - Inputs go to `src/cpp/data/`, and outputs go to `src/cpp/sol/`.
	 - Examples can be found in `src/cpp/examples/`.
	 - For best results, name your imputation code `missingX.o`, where X is a number. Then, in `experiment.py`, the methods list can be set as:
	 >\# impute using imPhy/src/cpp/missingX.o
	 >methods = [X]
	 >
	 >\# impute using imPhy/src/cpp/missingX.o and 
	 > \# imPhy/src/cpp/missing1.o
	 >methods = [X, 1]
4. Analysis
   - In the analysis section, DendroPy is used to take the Robinson-Foulds distance between original trees and their imputed siblings, while Owen and Provost's GTP code is used to calculate the BHV geodesic. Files containing information about these distances are written to `imPhy/my_experiment/batch_A/stats/`. If the next step is impossible to run, these files can still provide interesting information after a bit of file gymnastics. Headers are included in the CSVs.
5. Statistics
   - This section creates CSV files containing information about the distances between imputed and original leaves and trees. It can also create some diagnostic plots in R, and make heatmaps of some small instances, to solidify imputation quality intuition if the [dependencies](#Dependencies) are installed. Files created in this step are located in `imPhy/my_experiment/`.

Other Features:

1. Multiprocessing:
 imPhy uses Python's multiprocessing library to reduce total computation time. Parallel imputation methods have not been tested with the 'parallel' flag in the code, so caution is advised when using imputation methods that use more than one process.
2. Logging:
When big jobs are run on servers, it can be difficult to identify failure points. For this purpose, `imPhy/my_experiment/output.log` is placed inside the experiment folder, which contains sequentially all the output created by each process during the run.


Dependencies
-
This code has only been tested on Python 3.5 and R 3.3. Other functions *may* break when using previous versions. [Anaconda][conda] is recommended to install Python and its packages. It can also be used to install R, but is not as well supported as Anaconda for Python.

Python Packages:
 - [DendroPy][dp]
 - [docopt][docopt_py]
 - [numpy][np]
 - [matplotlib][plt]
 - [pandas][pd]
 - [seaborn][sns]

R Packages:
 - [APE][ape]
 - [car][car]
 - [corrplot][corrplot]
 - [docopt][docopt_r]
 - [lattice][lattice]
 - [reshape2][reshape2]

Others:
 - [Gurobi][gurobi]
 - [Java][java] (at least SE 1.8.0_31)

File Structure
-
```
imPhy
│   README
│
└───my_experiment
│   │   interleaf_error.csv
│   │   interleaf_error.pdf
│   │   intertree_all.csv
│   │   intertree_error.csv
│   │   intertree_error.pdf
│   │   outlier_counts.pdf
│   │   output.log
│   │
│   └───batch_A
│   │   └───data
│   │   │	│   data_file_1.txt.gz
│   │   │	│   data_file_1_true.txt.gz
│   │   └───nexus
│   │   │	│   gene_trees_1.nex
│   │   │	│   separated.txt
│   │   │	│   species.nex
│   │   └───solutions
│   │   │	│   gene_trees_1.sol.gz
│   │   └───stats
│   │   │	│   gene_trees_1_all.txt.gz
│   │   │	│   gene_trees_1_tree_all.txt.gz
│   │   │	│   gene_trees_1_tree.txt.gz
│   │   │	│   gene_trees_1.txt.gz
│   └───batch_B
│   │
│   ...
│   
└───src
    │
	...
```



  [ape]: https://cran.r-project.org/web/packages/ape/index.html "APE"
  [car]: https://cran.r-project.org/web/packages/car/index.html "car"
  [conda]: https://www.continuum.io/anaconda-overview "Anaconda"
  [corrplot]: https://cran.r-project.org/web/packages/corrplot/index.html "corrplot"
  [cpp]: https://isocpp.org/ "C++"
  [dp]: https://pythonhosted.org/DendroPy/ "DendroPy"
  [docopt_py]: http://docopt.org/ "docopt for Python"
  [docopt_r]: https://github.com/docopt/docopt.R "docopt for R"
  [gurobi]: http://www.gurobi.com/ "Gurobi"
  [java]: https://java.com/en/ "Java"
  [lattice]: https://cran.r-project.org/web/packages/lattice/index.html "lattice"
  [np]: http://www.numpy.org/ "numpy"
  [plt]: http://matplotlib.org/ "matplotlib"
  [pd]: http://pandas.pydata.org/ "pandas"
  [python]: https://www.python.org/ "python"
  [r]: https://www.r-project.org/ "r"
  [reshape2]: https://cran.r-project.org/web/packages/reshape2/index.html "reshape2"
  [sns]: https://stanford.edu/~mwaskom/software/seaborn/ "seaborn"


