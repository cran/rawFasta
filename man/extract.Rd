\name{extract}
\alias{extract}
\alias{extract-methods}
\alias{extract,rawFasta.8-method}
\alias{extract,rawFasta.dl-method}
\alias{extract,rawFasta.sl-method}
\title{
  Extracts a sub-sequence
}

\description{
  This generic function is dedicated to extract a sub-sequence from a large sequence object, such as \code{\linkS4class{rawFasta}} objects.
}

\usage{
  extract(object, start, end, check = TRUE, ...)
  
  \S4method{extract}{rawFasta.8}(object, start, end, check = TRUE, ...)
  \S4method{extract}{rawFasta.sl}(object, start, end, check = TRUE, ...)
  \S4method{extract}{rawFasta.dl}(object, start, end, check = TRUE, ...)
}

\arguments{
  \item{object}{
    An object to extract from, typically a \code{\linkS4class{rawFasta}}-inheriting object.
}
  \item{start}{
    Single integer value, the starting position of the sub-sequence to extract in the whole sequence object. First element of the sequence is '1', and \code{start} is part of the sub-sequence. Negative indexes are taken from the end of the sequence : '0' is the last element, '-1' the element before and so on.
}
  \item{end}{
    Single integer value, the ending position of the sub-sequence to extract in the whole sequence object. \code{end} is part of the sub-sequence. Negative indexes are taken from the end of the sequence : '0' is the last element, '-1' the element before and so on.
}
  \item{check}{
    Single logical value, whether to check \code{start} and \code{end} validaty and perform negative conversion or not. It can be time-saving to turn checking off for large batches of extractions, at your own risks.
}
  \item{\dots}{
    Further arguments to be dispatched to methods.
}
}

\value{
  Returns the sub-sequence as a variable length vector. Objects with explicit alphabets (such as \code{\linkS4class{rawFasta.sl}} and \code{\linkS4class{rawFasta.dl}}) return a \code{factor}, while \code{\linkS4class{rawFasta.8}} returns a vector of \code{character}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{rawFasta}}
}

\examples{
  ### Typical DNA FASTA file
  fastaFile <- system.file("extdata/dna.fa", package="rawFasta")
  f8 <- rawFasta(fastaFile, nbits=8)
  f4 <- rawFasta(fastaFile, alpha="ACGTN-", nbits=4)
  f3 <- rawFasta(fastaFile, alpha="ACGTN-")
  
  # Forward extraction
  print(extract(f8, 18, 71))
  print(extract(f4, 18, 71))
  print(extract(f3, 18, 71))
  
  # Backward extraction (10 last characters)
  print(extract(f8, -10, 0))
  print(extract(f4, -10, 0))
  print(extract(f3, -10, 0))
  
  # Consistency checking
  for(from in 1:10) {
    for(to in from:20) {
      x8 <- extract(f8, from, to)
      x4 <- as.character(extract(f4, from, to))
      x3 <- as.character(extract(f3, from, to))
      if(!identical(x8, x4)) stop("Inconsistency between x8 and x4")
      if(!identical(x8, x3)) stop("Inconsistency between x8 and x3")
    }
  }
}

