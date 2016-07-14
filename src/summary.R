'Usage: summary.R <infolder>

Options:
  -h, --help            Show this help message.
  <infolder>            Folder containing data, nexus, and solutions
                        subfolders. Must end in "/".' -> doc
library('corrplot')
library('docopt')
library('lattice')
opts <- docopt(doc)

# folders <- list.files(args["<infolder>"], include.dirs=TRUE)

data <- read.csv(paste0(opts$infolder,"summary.csv"), header=TRUE)

pdf(paste0(opts$infolder, "summary_plots.pdf"))

corrplot(cor(data[,-1][,-4], use="complete.obs"), type='upper')

myboxplots <- function(x){
	boxplot(data$Imputation.RMSE~x)
}

plotdata <- data[,-(8:12)][,-1][,-4]

for (i in 1:length(plotdata)) {
	if (i < length(plotdata)){
		boxplot(data[,13]~plotdata[,i], ylab='RMSE', main=paste("RMSE by", names(plotdata[i])), xlab=names(plotdata[i]), type="l")
    } else {
    	boxplot(data[,13], ylab='RMSE', main="RMSE boxplot for Python Data")
    }
}

dev.off()