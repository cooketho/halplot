Code for reproducing the HAL graph shown in Figure S7 of Cooke et al 2017

C++ dependencies:
HAL: github.com/ComparativeGenomicsToolkit/hal
sonLib: github.com/benedictpaten/sonLib

Python dependencies:
numpy: http://www.numpy.org
scipy: https://www.scipy.org
newick: https://pypi.python.org/pypi/newick

R dependencies:
ggplot2: http://ggplot2.org
plyr: http://plyr.had.co.nz
reshape2: https://github.com/hadley/reshape
RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
IRanges: http://bioconductor.org/packages/release/bioc/html/IRanges.html

###############################################################################################
Installation:

Probably the easiest way to install the C++ dependencies is to install progressive cactus,
which automatically installs them as submodules:

git clone git://github.com/glennhickey/progressiveCactus.git
cd progressiveCactus
git pull
git submodule update --init
make

Next, change the following line of make.sh to the path to your progressive cactus installation:
PGPATH=<path to your progressive cactus installation>

Fetch the code for intervaltree:
git submodule init
git submodule update

Then, to compile haltraverse:
./make.sh

To use haltraverse, be sure to add the following path to the environmental variable LD_LIBRARY_PATH:
progressiveCactus/submodules/hdf5/lib/

###############################################################################################
Use:

Align sequences with progressive cactus (Already done for you--see wga6.hal):
runProgressiveCactus.sh wga6.txt wga6 wga6/wga6.hal --stats &

Traverse HAL graph and print edges:
haltraverse --bed wga6.bed --root Anc0 wga6.hal > wga6.edges

Get unique set of edges:
cat wga6.edges | sort | uniq > wga6.unique.edges
cut -f 1-9 wga6.unique.edges | sort | uniq > wga6.no_features.edges

Get the newick-formatted tree used in the PG alignment:
halStats --tree wga6.hal > wga6.newick

Find the optimal order and orientation of sequences for plotting:
python count_inversions_greedy.py wga6.newick wga6.no_features.edges > wga6.config

Plot the HAL graph by running the R code contained in the following file:
halplot.r
