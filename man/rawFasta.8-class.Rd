\name{rawFasta.8-class}
\Rdversion{1.1}
\docType{class}
\alias{rawFasta.8-class}

\title{Class \code{"rawFasta.8"}}
\description{
  This \code{\linkS4class{rawFasta}} implementation stores sequence elements on 8 bits, using usual ASCII alphabet. It is the faster implementation, but also the less memory efficient as the 256 available characters are hardly used.
}
\section{Objects from the Class}{
  Objects can be created by \code{\link{rawFasta}}, with \code{nbits} = 8.
}
\section{Extends}{
  Class \code{"\linkS4class{rawFasta}"}, directly.
}
\author{
  Sylvain Mareschal
}
\seealso{
  \code{link{extract}}, \code{rawFasta}, \code{\link{rawFasta-class}}, \code{\link{rawFasta.sl-class}}, \code{\link{rawFasta.dl-class}}
}
\keyword{classes}
