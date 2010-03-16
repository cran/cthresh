"wst"<-
function(data, filter.number = 10, family = "DaubLeAsymm", verbose = FALSE)
{
    if(verbose == TRUE)
        cat("Argument checking...")
    DataLength <- length(data)  #
#
# Check that we have a power of 2 data elements
#
    nlevels <- log(DataLength)/log(2)
    if(round(nlevels) != nlevels)
        stop("The length of data is not a power of 2")  #
    if(verbose == TRUE) {
        cat("There are ", nlevels, " levels\n")
    }
#
# Select the appropriate filter
#
    if(verbose == TRUE)
        cat("...done\nFilter...")
    filter <- filter.select(filter.number = filter.number, family = family)
#
#
# Compute the decomposition
#
    if(verbose == TRUE)
        cat("Decomposing...\n")
    newdata <- c(rep(0, DataLength * nlevels), data)
    Carray <- newdata
    error <- 0  #
#
#   See whether we are using complex wavelets
#
    if(is.null(filter$G)) {
        wavelet.station <- .C("wavepackst",
            Carray = as.double(Carray),
            newdata = as.double(newdata),
            DataLength = as.integer(DataLength),
            levels = as.integer(nlevels),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            error = as.integer(error), PACKAGE  = "wavethresh")
    }
    else {
        wavelet.station <- .C("comwst",
            CaR = as.double(Re(Carray)),
            CaI = as.double(Im(Carray)),
            newdataR = as.double(Re(newdata)),
            newdataI = as.double(Im(newdata)),
            DataLength = as.integer(DataLength),
            levels = as.integer(nlevels),
            HR = as.double(Re(filter$H)),
            HI = as.double( - Im(filter$H)),
            GR = as.double(Re(filter$G)),
            GI = as.double( - Im(filter$G)),
            LengthH = as.integer(length(filter$H)),
            error = as.integer(error), PACKAGE = "wavethresh")
                }
    if(wavelet.station$error != 0)
        stop(paste("Memory error in wavepackst (or comwst). Code ", 
            wavelet.station))
    if(is.null(filter$G)) {
        wpm <- matrix(wavelet.station$newdata, ncol = DataLength, byrow
             = TRUE)
        Carray <- matrix(wavelet.station$Carray, ncol = DataLength, 
            byrow = TRUE)
    }
    else {
        newdata <- complex(real = wavelet.station$newdataR, im = 
            wavelet.station$newdataI)
        Carray <- complex(real = wavelet.station$CaR, im = 
            wavelet.station$CaI)
        wpm <- matrix(newdata, ncol = DataLength, byrow = TRUE)
        Carray <- matrix(Carray, ncol = DataLength, byrow = TRUE)
    }
    wp <- list(wp = wpm, Carray = Carray, nlevels = nlevels, filter = 
        filter, date = date())
    class(wp) <- "wst"
    wp
}
