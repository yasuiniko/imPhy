
'
This file is part of imPhy, a pipeline for evaluating the quality of
phylogenetic imputation software.
Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis

imPhy uses GTP, which is Copyright © 2008, 2009  Megan Owen, Scott Provan

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: outliers.R <exp_folder>

Options:
  -h, --help            Show this help message.
  <exp_folder>          Experiment folder. Must not end in "/".' -> doc

library('car')
library('corrplot')
suppressMessages(library('docopt'))
library('lattice')
opts <- docopt(doc)

data <- read.csv(paste0(opts$exp_folder,"/interleaf_error.csv"), header=TRUE)

get_outliers_of <- function(data) { 
    function(field) {
        main_field <- get("Imputation.NRMSE", data)
        grouping_field <- get(field, data)
        outliers <- Boxplot(main_field, grouping_field)
        data[outliers,][c(2:9, length(data))]
    }
}

extreme_patterns_in <- function (data) {
    function (field) {
        extreme_values <- function (x) {max(x) > 1.5*min(x)}
        Filter(extreme_values, lapply(get_outliers_of(data)(field), table))
    }
}

plot_to <- function(out){
    function(data){
        pdf(out)

        max_ind <- length(data)

        for (i in 1:(max_ind-1)){
            barplot(table(data[,i]),
                    ylab='Counts',
                    main=paste("Outlier Counts by", names(data[i])),
                    xlab=names(data[i]))
        }
        dev.off()
    }
}

suppressWarnings(corrplot(cor(data[,-1][,-4], use="complete.obs"), type='upper'))

d <- Reduce(function (...) merge(..., all=T), 
            lapply(names(data)[-(8:14)][-1], get_outliers_of(data)))

plot_to(paste0(opts$exp_folder, "/outlier_counts.pdf"))(d)

invisible(file.remove("Rplots.pdf"))