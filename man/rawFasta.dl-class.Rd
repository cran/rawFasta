\name{rawFasta.dl-class}
\Rdversion{1.1}
\docType{class}
\alias{rawFasta.dl-class}

\title{Class \code{"rawFasta.dl"}}
\description{
  This \code{\linkS4class{rawFasta}} implementation stores sequence elements on two layers of distinct bit counts, allowing odd bit counts in alphabet length (storing elements on 3 bits is actually managed by storing the two first bits and the last one is separate containers, to simplify extraction). This implementation is the slowest, but provides more precision in alphabet allocation size to achieve better memory savings.
}
\section{Objects from the Class}{
  Objects can be created by \code{\link{rawFasta}}, with \code{nbits} = 3 (2+1), 5 (4+1) or 6 (4+2).
}
\section{Slots}{
  \describe{
    \item{\code{alpha}:}{Character vector holding the allowed characters in the used alphabet.}
    \item{\code{table}:}{List of two integer matrixes, each of them is used to convert octets into elements sequences (the 256 octets are used as column inputs to extract the indexes in \code{alpha} of the elements compacted in it).}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{rawFasta}"}, directly.
}
\author{
  Sylvain Mareschal
}
\seealso{
  \code{link{extract}}, \code{rawFasta}, \code{\link{rawFasta-class}}, \code{\link{rawFasta.8-class}}, \code{\link{rawFasta.sl-class}}
}
\keyword{classes}
