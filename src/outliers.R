
'Usage: outliers.R <infolder>

Options:
  -h, --help            Show this help message.
  <infolder>            Folder containing data, nexus, and solutions
                        subfolders. Must not end in "/".' -> doc
library('car')
library('corrplot')
suppressMessages(library('docopt'))
library('lattice')
opts <- docopt(doc)

data <- read.csv(paste0(opts$infolder,"/summary.csv"), header=TRUE)

get_outliers_of <- function(data) { 
    function(field) {
        main_field <- get("Imputation.NRMSE", data)
        grouping_field <- get(field, data)
        outliers <- Boxplot(main_field, grouping_field)
        data[outliers,][-(8:13)][-1]
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
                    main=paste("Number of outliers by", names(data[i])),
                    xlab=names(data[i]))
        }
        dev.off()
    }
}

corrplot(cor(data[,-1][,-4], use="complete.obs"), type='upper')

d <- Reduce(function (...) merge(..., all=T), 
            lapply(names(data)[-(8:14)][-1], get_outliers_of(data)))

plot_to(paste0(opts$infolder, "/summary_plots_outliers.pdf"))(d)

invisible(file.remove("Rplots.pdf"))