\name{rawFasta.sl-class}
\Rdversion{1.1}
\docType{class}
\alias{rawFasta.sl-class}

\title{Class \code{"rawFasta.sl"}}
\description{
  This \code{\linkS4class{rawFasta}} implementation stores sequence elements on less than 8 bits, to achieve memory savings. This implementation is slower than \code{\linkS4class{rawFasta.8}}, but large sequences can be stored in smaller memory spaces.
}
\section{Objects from the Class}{
  Objects can be created by \code{\link{rawFasta}}, with \code{nbits} = 1, 2 or 4.
}
\section{Slots}{
  \describe{
    \item{\code{alpha}:}{Character vector holding the allowed characters in the used alphabet.}
    \item{\code{table}:}{An integer matrix, used to convert octets into elements sequences (the 256 octets are used as column inputs to extract the indexes in \code{alpha} of the elements compacted in it).}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{rawFasta}"}, directly.
}
\author{
  Sylvain Mareschal
}
\seealso{
  \code{link{extract}}, \code{rawFasta}, \code{\link{rawFasta-class}}, \code{\link{rawFasta.8-class}}, \code{\link{rawFasta.dl-class}}
}
\keyword{classes}
