## c:/Dropbox/Rpackages/crossplotr/R/crossplot_outliers.R

##    Chandler Lutz
##    Questions/comments: cl.eco@cbs.dk
##    $Revisions:      1.0.0     $Date:  2017-04-25

##To get the data corresponding to outliers in certain variables


#' To get the index of the outliers for a vector x
#'
#' This function will return a logical index where outliers are set to
#' \code{TRUE}. The user can imput either imput the number of desired
#' outliers, \code{num.outliers}, or the percentile corresponding to
#' the outliers \code{percentile.outliers}, but not both.
#'
#' Note: \code{NA} values will be replaced with the mean value.
#'
#' @param x a numeric vector
#' @param num.outliers the number of outliers. Defaults to \code{NULL}
#' @param percentile.outliers the percentile used to measure
#'     outliers. Defaults to \code{NULL}
#' @return a logical index where outliers are set to \code{TRUE}
get_outlier_logicals <- function(x, num.outliers = NULL, percentile.outliers = NULL) {

   if (!is.numeric(x)) stop("Error: x must be numeric")

   if (!is.null(num.outliers) & !is.null(percentile.outliers))
      stop("Error: get_outlier_logical() only uses either num.outliers or percentile.outliers")

   ##Replace NAs with the mean value
   x[is.na(x)] <- rep(mean(x, na.rm = TRUE), sum(is.na(x)))

   ##Sort x
   x.sorted <- sort(x)
   ##the length of x
   nn <- length(x)

   if (!is.null(num.outliers)) {
      ##Use the number of outliers
      outliers <- c(x.sorted[1:num.outliers], x.sorted[(nn - num.outliers + 1):nn])
      which(x %in% outliers)
   } else if (!is.null(percentile.outliers)) {
      ##use the percentile to get the outliers

      ##Make sure percentile.outliers is between 0 and 1
      if (percentile.outliers < 0 | percentile.outliers > 1 )
         stop("Error: percentile.outliers must be between 0 and 1")

      ##The number of outliers to be used if the percentile is set
      num.percentile.outliers <- round(nn * percentile.outliers)
      outliers <- c(x.sorted[1:num.percentile.outliers],
                    x.sorted[(nn - num.percentile.outliers):nn])

   } else {
      stop("Error: get_outlier_logical() requires either num.outliers or percentile.outliers to be set")
   }

   ##Return the logicals for the outliers
   return(x %in% outliers)

}


#' To get the outlier data form a ggplot
#'
#' This function is useful for getting outliers and labeling outlying
#' points on a plot. Based a \code{data.frame} or \code{ggplot} and
#' variables from the user, \code{crossplot_outliers} returns the
#' corresponding outliers. The user can input either the number of
#' desired outliers, \code{num.outliers}, or the percentile
#' corresponding to the outliers \code{percentile.outliers}, but not
#' both. Note: \code{NA} values will be replaced with the mean value.
#'
#' @param x a \code{ggplot2} plot or a \code{data.frame}. Can be
#'     generated from \code{crossplot}.
#' @param vars the variables for which the outliers will be
#'     calculated. If \code{x} is a \code{ggplot2} object, use the
#'     variable names in the underlying data
#' @param num.outliers the number of outliers for each variable in
#'     \code{vars}. Defaults to \code{NULL}
#' @param percentile.outliers the percentile used to measure outliers
#'     for each variable in \code{vars}. Defaults to \code{NULL}
#' @return a \code{data.frame} with the plot data for the outliers
#' @examples
#' data(mtcars)
#' ##Get the outliers for hp and and mpg from mtcars
#' crossplot_outliers(mtcars,vars = c("hp","mpg"), percentile.outliers = 0.05)
#' crossplot_outliers(mtcars,vars = c("hp","mpg"), num.outliers = 2)
#'
#' ##Outliers from a ggplot
#' p <- ggplot(mtcars, aes(x = hp, y = mpg)) + geom_point(aes(size = wt))
#' crossplot_outliers(p, vars = c("hp", "mpg"), num.outliers = 2)
#' crossplot_outliers(p, vars = c("mpg", "hp", "wt"), percentile.outliers = .05)
#'
#' ## -- Full Example using crossplot() and labels -- ##
#' mtcars$name <- rownames(mtcars)
#' mtcars.outliers <- crossplot_outliers(mtcars, vars = c("mpg", "hp"), num.outliers = 2)
#' p <- crossplot(mtcars, x.var = "mpg", y.var = "hp", size.var = "wt",
#'               shapes.var = "cyl", label.var = "name", points.alpha = 0.7) +
#'      geom_text(data = mtcars.outliers)
#' print(p)
#'
#' ##If there are missing values in vars, these will be replaced by the
#' ##mean so they are less likely to affect the outliers
#' data(mtcars)
#' mtcars$mpg[7:12] <- rep(NA, 6)
#' ##Get the outliers for hp and and mpg from mtcars
#' crossplot_outliers(mtcars,vars = c("hp","mpg"), percentile.outliers = 0.05)
#' crossplot_outliers(mtcars,vars = c("hp","mpg"), num.outliers = 2)
#' @export
crossplot_outliers <- function(x, vars = NULL, num.outliers = NULL,
                               percentile.outliers = NULL) {


   if (!is.null(num.outliers) & !is.null(percentile.outliers))
      stop("Error: crossplot_outliers() only uses either num.outliers or percentile.outliers")

   if (is.ggplot(x)) {
      ##ggplot2 object -- get data from the plot
      all.data <- x$data
   } else if (is.data.frame(x)) {
      all.data <- as.data.frame(x)

   } else {
      stop("Error: x must be either a ggplot or a data.frame")
   }

   ##The data variables
   data.vars <- names(all.data)

   ##Keep only the variables requested by the used
   data.vars <- data.vars[data.vars %in% vars] %>% as.character

   if (!is.null(num.outliers)) {
      f_get_outliers <- function(x) get_outlier_logicals(x, num.outliers)
   } else if (!is.null(percentile.outliers)) {
      f_get_outliers <- function(x) get_outlier_logicals(x, percentile.outliers)
   } else {
      stop("Error: crossplot_outliers() requires either num.outliers or percentile.outliers to be set")
   }

   ##Get just the data that we need for the outlier indices
   outliers.logicals <- all.data[, data.vars, drop = FALSE]
   ##apply the f_get_outliers function
   outliers.logicals[] <- lapply(outliers.logicals, f_get_outliers)
   ##sum across rows and get a logical
   outliers.logicals <- rowSums(outliers.logicals) %>% as.logical

   ##Return the data
   all.data.outliers <- all.data[outliers.logicals, ]

   return(all.data.outliers)

}
