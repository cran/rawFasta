\name{rawFasta-class}
\Rdversion{1.1}
\docType{class}
\alias{rawFasta-class}
\alias{show,rawFasta-method}

\title{Class \code{"rawFasta"}}
\description{
  This interface is dedicated to hold in memory the compacted content of a FASTA sequence file, in order to make several extractions. It is implemented by several classes, like \code{\linkS4class{rawFasta.8}}, \code{\linkS4class{rawFasta.sl}} and \code{\linkS4class{rawFasta.dl}}.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{nbits}:}{Integer vector, the amount of bits allocated to each sequence element. Notice \code{\linkS4class{rawFasta.dl}} have two integer values in this vector, as the elements are stored on two distinct layers.}
    \item{\code{header}:}{Single character string, the comment line extracted from the FASTA file (without '>').}
    \item{\code{source}:}{Single character string, the absolute path to the FASTA file parsed to produce the object.}
    \item{\code{content}:}{An environment holding the content variable(s), to avoid memory duplication on copy.}
    \item{\code{size}:}{Single integer value, the amount of elements in the FASTA sequence parsed.}
  }
}
\author{
  Sylvain Mareschal
}
\seealso{
  \code{link{extract}}, \code{rawFasta}, \code{\link{rawFasta.8-class}}, \code{\link{rawFasta.sl-class}}, \code{\link{rawFasta.dl-class}}
}
\examples{
showClass("rawFasta")
}
\keyword{classes}
