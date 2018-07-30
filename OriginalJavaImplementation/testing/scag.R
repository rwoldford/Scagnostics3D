library(rJava)

## update this list when scagnostics are added to the Java code!
.scagnostics3D.names <- c("Outlying", "Skewed", "Clumpy", "Sparse", "Striated", "Convex", "Skinny", "Stringy", "Monotonic")

scagnostics3D <- function(x, ...) UseMethod("scagnostics3D", x)

## calls scagnostics Java code for a given data which must be a Java reference of the class "[[D"
.scagnostics3D.for.data <- function(data, names=NULL) {
  results <- .jcall("scagnostics3D/Main","[[D","computeScagnostics",data)
  r <- lapply(results, .jevalArray)
  n <- length(r)
  s <- length(r[[1]])
  if (n == 1) {
    r <- r[[1]]
    names(r) <- .scagnostics3D.names
  } else {
    r <- matrix(unlist(r), s)
    rownames(r) <- .scagnostics3D.names
    if (!is.null(names)) {
      dn <- length(names)
      l <- data.frame(x=rep(1:dn,(dn*dn)),y=rep(1:dn,each=dn),z=rep(1:dn,each=dn*dn))
	  l <- unique(l[l$x < l$y, ]) # collapse to those satisfying x < y
	  l <- unique(l[l$y < l$z, ]) # AND those satisfying y < z	
      colnames(r) <- paste(names[l$x], "*", names[l$y], "*", names[l$z])
    }
  }
  class(r) <- "scagnostics3D"
  r
}

## --- scagnostics3D methods for different data types ---

scagnostics3D.default <- function(x, y, z,...) {
  if (length(x) != length(y)) stop("x and y must have the same length")
  if (length(y) != length(z)) stop("y and z must have the same length")

  complete <- !is.na(x) & !is.na(y)  
  x <- as.double(x[complete])
  y <- as.double(y[complete])
  z <- as.double(z[complete])
  
  data <- .jarray(list(.jarray(x),.jarray(y),.jarray(z)),"[D")
  .scagnostics3D.for.data(data)
}

scagnostics3D.data.frame <- function(x, ...) {
  if (dim(x)[2] < 2) stop("need at least three variables")
  data <- .jarray(lapply(x, function(x) .jarray(as.double(x))), "[D")
  .scagnostics3D.for.data(data, names(x))
}

scagnostics3D.list <- function(x, ...) {
  if (length(x) < 3) stop("need at least three variables")
  n <- unlist(lapply(x, length))
  if (!all(n == n[1])) stop("all variables must have the same length")
  data <- .jarray(lapply(x, function(x) .jarray(as.double(x))), "[D")
  nam <- names(x)
  if (is.null(nam)) nam <- as.character(1:length(x))
  .scagnostics3D.for.data(data, nam)
}

scagnostics3D.matrix <- function(x, ...) {
  if (dim(x)[2] < 2) stop("need at least three variables")
  data <- .jarray(lapply(1:(dim(x)[2]), function(i) .jarray(as.double(x[,i]))), "[D")
  nam <- colnames(x)
  if (is.null(nam)) nam <- as.character(1:length(x))
  .scagnostics3D.for.data(data, nam)
}

#################################################################################
# testing
require(rggobi)
require(scatterplot3d)

## quote
newwindow <- function(...) {
	gui <- .Platform$GUI
	if (gui=="AQUA") quartz(...) else
	  if (gui=="X11") X11(...) else
	   if (gui=="Rgui") windows(...) else {
	   	 print("Sorry.  Can't tell what the GUI is to pop the correct window. Will just use plot(1)")
	   	plot(1)}
	}
## end


.jinit()

# for Mac
#.jaddClassPath("/tmp/scagnostics3D/scagnostics3D.jar")
.jaddClassPath("Mac/scagnostics3D.jar")

# for Windows
#.jaddClassPath("C:/scagnostics3D/scagnostics3D.jar")

# test... to make sure the JAR can be loaded
stest <- scagnostics3D(rnorm(10), rnorm(10), rnorm(10))
stest
stest <- scagnostics3D(iris[,1:4])
library(cluster)

stest <- scagnostics3D(chorSub[,1:6])