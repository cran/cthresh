"wr.wd"<-
function(wd, start.level = 0, verbose = FALSE, bc = wd$bc, return.object = FALSE, 
    filter.number = wd$filter$filter.number, family = wd$filter$family, ...)
{
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    if(verbose == TRUE) cat("Argument checking...") #
#
#       Check class of wd
#
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(wd)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")
    if(start.level < 0)
        stop("start.level must be nonnegative")
    if(start.level >= nlevels(wd))
        stop("start.level must be less than the number of levels")
    if(is.null(wd$filter$filter.number))
        stop("NULL filter.number for wd")
    if(bc != wd$bc)
        warning("Boundary handling is different to original")
    if(wd$type == "station")
        stop("Use convert to generate wst object and then AvBasis or InvBasis"
            )
    if(wd$bc == "interval") {
        warning("All optional arguments ignored for \"wavelets on the interval\" transform"
            )
        return(wr.int(wd))
    }
    type <- wd$type
    filter <- filter.select(filter.number = filter.number, family = family)
    LengthH <- length(filter$H) #
#
#   Build the reconstruction first/last database
#
    if(verbose == TRUE)
        cat("...done\nFirst/last database...")
    r.first.last.c <- wd$fl.dbase$first.last.c[(start.level + 1):(wd$
        nlevels + 1),  ]    #
    r.first.last.d <- matrix(wd$fl.dbase$first.last.d[(start.level + 1):(wd$
        nlevels),  ], ncol = 3)
    ntotal <- r.first.last.c[1, 3] + r.first.last.c[1, 2] - r.first.last.c[
        1, 1] + 1
    names(ntotal) <- NULL
    C <- accessC(wd, level = start.level, boundary = TRUE)
    C <- c(rep(0, length = (ntotal - length(C))), C)
    Nlevels <- nlevels(wd)- start.level
    error <- 0  #
#
#   Load object code
#
    if(verbose == TRUE)
        cat("...built\n")
    if(verbose == TRUE) {
        cat("Reconstruction...")
        error <- 1
    }
    ntype <- switch(type,
        wavelet = 1,
        station = 2)
    if(is.null(ntype))
        stop("Unknown type of decomposition")
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary handling")
    if(!is.complex(wd$D)) {
        wavelet.reconstruction <- .C("waverecons",
            C = as.double(C),
            D = as.double(wd$D),
            H = as.double(filter$H),
            LengthH = as.integer(LengthH),
            nlevels = as.integer(Nlevels),
            firstC = as.integer(r.first.last.c[, 1]),
            lastC = as.integer(r.first.last.c[, 2]),
            offsetC = as.integer(r.first.last.c[, 3]),
            firstD = as.integer(r.first.last.d[, 1]),
            lastD = as.integer(r.first.last.d[, 2]),
            offsetD = as.integer(r.first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    else {
        wavelet.reconstruction <- .C("comwr",
            CR = as.double(Re(C)),
            CI = as.double(Im(C)),
            LengthC = as.integer(length(C)),
            DR = as.double(Re(wd$D)),
            DI = as.double(Im(wd$D)),
            LengthD = as.integer(length(wd$D)),
            HR = as.double(Re(filter$H)),
            HI = as.double(Im(filter$H)),
            GR = as.double(Re(filter$G)),
            GI = as.double(Im(filter$G)),
            LengthH = as.integer(LengthH),
            nlevels = as.integer(Nlevels),
            firstC = as.integer(r.first.last.c[, 1]),
            lastC = as.integer(r.first.last.c[, 2]),
            offsetC = as.integer(r.first.last.c[, 3]),
            firstD = as.integer(r.first.last.d[, 1]),
            lastD = as.integer(r.first.last.d[, 2]),
            offsetD = as.integer(r.first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    if(verbose == TRUE)
        cat("done\n")
    error <- wavelet.reconstruction$error
    if(error != 0) {
        cat("Error code returned from waverecons: ", error, "\n")
        stop("waverecons returned error")
    }
    fl.dbase <- wd$fl.dbase
    if(!is.complex(wd$D)) {
        l <- list(C = wavelet.reconstruction$C, D = 
            wavelet.reconstruction$D, fl.dbase = fl.dbase, nlevels
             = nlevels(wd), filter = filter, type = type, bc = bc, 
            date = date())
    }
    else {
        l <- list(C = complex(real = wavelet.reconstruction$CR, im = 
            wavelet.reconstruction$CI), D = complex(real = 
            wavelet.reconstruction$DR, im = wavelet.reconstruction$
            DI), fl.dbase = fl.dbase, nlevels = nlevels(wd), filter
             = filter, type = type, bc = bc, date = date())
    }
    class(l) <- "wd"
    if(return.object == TRUE)
        return(l)
    else return(accessC(l))
    stop("Shouldn't get here\n")
}
