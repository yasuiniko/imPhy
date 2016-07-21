'Usage: summary.R <infolder>

Options:
  -h, --help            Show this help message.
  <infolder>            Folder containing data, nexus, and solutions
                        subfolders. Must not end in "/".' -> doc
library('corrplot')
suppressMessages(library('docopt'))
library('lattice')
opts <- docopt(doc)

# folders <- list.files(args["<infolder>"], include.dirs=TRUE)

data <- read.csv(paste0(opts$infolder,"/summary.csv"), header=TRUE)

pdf(paste0(opts$infolder, "/summary_plots.pdf"))

corrplot(cor(data[,-1][,-4], use="complete.obs"), type='upper')

plotdata <- data[,-(8:13)][,-1][,-4]

for (i in 1:length(plotdata)) {
    if (i < length(plotdata)){
        boxplot(data[,14]~plotdata[,i], ylab='NRMSE', main=paste("NRMSE by", names(plotdata[i])), xlab=names(plotdata[i]), type="l")
    } else {
        boxplot(data[,14], ylab='NRMSE', main="NRMSE boxplot")
    }
}

dev.off()

