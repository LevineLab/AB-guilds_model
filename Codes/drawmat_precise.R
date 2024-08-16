##' Another helper function to /precisely/ draw the entries of a matrix.
##' @param mat Matrix of interest.
##' @param contour If \code{TRUE}, draw a contour using
##'   \code{lattice::levelplot()}.
##' @param ... Other arguments to \code{lattice::levelplot()}.
##'
##' @return lattice object.
drawmat_precise <- function(mat, contour = FALSE, ...){
  
  if(is.null(colnames(mat))){
    colnames(mat) <- paste(rep("col\n",ncol(mat)),
                           c(1:ncol(mat)) , sep=" ")
  }
  if(is.null(rownames(mat))){
    rownames(mat) <- paste(rep("row",nrow(mat)),
                           c(1:nrow(mat)) , sep=" ")
  }
  
  ## Color function
  colfun = colorRampPalette(c("blue", "red"))
  
  #If need to plot specific boundaries like 0 to 1 use the levels and modify the lattice plot with the following: at=levels
  levels<-seq(from=0,to=0.8,by=.01)
  
  # plot it flipping the axis
  lattice::levelplot(t(mat[c(nrow(mat):1) , ]),
                     col.regions = colfun(100),
                     contour = contour,
                     ## xaxt = 'n',
                     las = 2,at=levels,
                     ...)
}