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
