#'@title CountChemicalElements.
#'@description \code{CountChemicalElements} will split a character (chemical formula)
#'  into its elements and count their occurence.
#'@details No testing for any chemical alphabet is performed. Elements may occur
#'  several times and will be summed up in this case without a warning.
#'  Information within brackets, i.e. [13]C will be removed prior to counting
#'  together with other symbols not contained in the R set 'alnum'.
#'  The result is filtered and oredered according to parameter 'ele' if provided.
#'@param x Chemical formula.
#'@param ele Character vector of elements to count particularly or counting all contained in string if NULL.
#'@return A named numeric with counts for all contained or specified elements.
#'@examples
#'# count every element
#'CountChemicalElements("C3H7Cl")
#'
#'# remove additional symbols and sum up redundant elements (here 'C')
#'CountChemicalElements("[13]CC2H8Cl+")
#'
#'# count specific elements and return in specified order
#'CountChemicalElements("[13]CC2H8Cl+", ele=c("Cl","O","H"))
#'
#'# apply on a vector of formulas using sapply
#'sapply(c("C3H7Cl", "[13]CC2H8Cl+"), CorMID::CountChemicalElements, ele=c("Cl","O","H"))
#'@export
CountChemicalElements <- function(x = NULL, ele = NULL) {
    # ensure character and length=1
    x <- as.character(x[1])
    # count all elements present within 'x'
    # remove square bracket constructs (e.g. [13]C6 --> C6) upfront
    x <- gsub("[[].+[]]","",x)
    # all elements start with a LETTER...
    p <- gregexpr("[[:upper:]]", x)[[1]]
    # split initial string at the large letter positions
    out <- sapply(1:length(p), function(i) {
        substr(x, p[i], ifelse(i == length(p), nchar(x), p[i + 1] - 1))
    })
    # remove all non letter/digit (e.g. further brackets, charges...)
    out <- gsub("[^[:alnum:]]", "", out)
    count <- as.numeric(gsub("[^[:digit:]]", "", out))
    count[is.na(count)] <- 1
    names(count) <- gsub("[^[:alpha:]]", "", out)
    # sum up in case that elements were found repeatedly
    if (any(duplicated(names(count)))) {
        for (i in rev(which(duplicated(names(count))))) {
            count[which((names(count) == names(count)[i]))[1]] <- count[which((names(count) == names(count)[i]))[1]] + count[i]
            count <- count[-i]
        }
    }
    # reorder or limit output vector according to 'ele' and 'order_ele'
    if (!is.null(ele)) count <- sapply(ele, function(e) { ifelse(e %in% names(count), count[names(count)==e], 0) })
    return(count)
}
