'
This file is part of imPhy, a pipeline for evaluating the quality of
phylogenetic imputation software.
Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis

imPhy uses GTP, which is Copyright © 2008, 2009  Megan Owen, Scott Provan

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: summary.R <infolder>

Options:
  -h, --help            Show this help message.
  <infolder>            Folder containing folders containingdata, nexus, 
                        and solutions subfolders. Must not end in "/".
' -> doc

library('corrplot')
suppressMessages(library('docopt'))
library('lattice')
library('reshape2')
opts <- docopt(doc)



data <- read.csv(paste0(opts$infolder,"/interleaf_error.csv"), header=TRUE)

pdf(paste0(opts$infolder, "/interleaf_error.pdf"))

plotdata <- data[,-(10:15)][,-1]

for (i in 1:length(plotdata)) {
    if (i < length(plotdata)){
        boxplot(data[,16]~plotdata[,i], ylab='NRMSE', main=paste("Interleaf NRMSE by", names(plotdata[i])), xlab=names(plotdata[i]), type="l")
    } else {
        boxplot(data[,16], ylab='NRMSE', main="Interleaf NRMSE Overview")
    }
}

suppressWarnings(corrplot(cor(data[,-1], use="complete.obs"), type='upper'))

invisible(dev.off())



data <- read.csv(paste0(opts$infolder,"/intertree_error.csv"), header=TRUE)
data.r <- reshape(data, 
                  direction="long",
                  varying=c("rf_nj", "rf_upgma", "bhv_nj", "bhv_upgma"),
                  v.names="nrmse",
                  sep="")

pdf(paste0(opts$infolder, "/intertree_error.pdf"))

plotdata <- data.r[,-1][-11]

for (i in 1:(length(plotdata) - 6)) {
    if (i < length(plotdata) - 6){
        boxplot(as.formula(paste0("nrmse~",names(plotdata)[i])),
                ylab='NRMSE',
                main=paste("Intertree NRMSE by", names(plotdata[i])),
                xlab=names(plotdata[i]),
                data=plotdata)
    } else {
        boxplot(data.r$nrmse, ylab='NRMSE', main="Intertree NRMSE Overview")
    }
}

invisible(dev.off())

