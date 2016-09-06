formalities <- '
This file is part of imPhy, a pipeline for evaluating the quality of
phylogenetic imputation software.
Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis

imPhy uses GTP, which is Copyright © 2008, 2009  Megan Owen, Scott Provan

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'

two_digits <- function(x) format(round(x, 2), nsmall=2)

valid_prob <- function(x) if (0<=x & x<=1) TRUE else FALSE

invalid_probability_error <- function(x) {
    paste("Error:", x, "is an invalid probability.")
}

assert_valid_prob <- function(x) {
    if (valid_prob(x)) {
        x
    } else {
        exit(invalid_probability_error(x))
    }
}

exit <- function(x) {
    print(x)
    quit(status=1)
}
