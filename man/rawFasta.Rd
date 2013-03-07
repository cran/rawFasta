\name{rawFasta}
\alias{rawFasta}
\title{
  Parses a FASTA file
}

\description{
  This function reads a FASTA sequence file and return its content as a \code{\linkS4class{rawFasta}} complient S4 object.
}

\usage{
  rawFasta(fileName, maxBuffer = 8192L, alpha = "ACGTN-", nbits = NA)
}

\arguments{
  \item{fileName}{
    Single character name, path and name of the FASTA file to read. See the 'Note' section for details on the FASTA format.
}
  \item{maxBuffer}{
    Single integer value, maximal size of the buffer to use for parsing (in R \code{\link{raw}} element sizes). Notice the buffer must be large enough to hold at least 8 lines of the parsed FASTA file, and large buffers are not necessary faster. The default value has been optimised, and should only be changed by advanced users.
}
  \item{alpha}{
    Single character string, the concatenated letters to allow in the sequence alphabet. Unexpected letters will raise errors, so define this argument with care. This argument is ignored when \code{nbits} is 8, as the full ASCII alphabet is used instead.
}
  \item{nbits}{
    Single numeric value, the amount of bits to allow to each letter (can be 1, 2, 3, 4, 5, 6 or 8). Use \code{NA} to use the smallest possible value according to \code{alpha}. See the 'Details' section for further details.
}
}

\details{
  The purpose of this function is to use \code{nbits} smaller than 8, in order to minimize memory consumption for sequences of few letters such as DNA or amino-acid sequences. Taking DNA as example, the four possible letters (A, C, G and T) can be stored on only 2 bits, as two bits leads to 4 possible values. This is achieved by storing 4 letters in each octet, which lead to a 4 fold memory saving according to full ASCII storage (8 bits per letter).
  
  In a more general way, n bits can hold 2^n distinct letters, so :
  \tabular{ccl}{
    \strong{nbits}	\tab	\strong{Values}	\tab	\strong{Example}	\cr
    1	\tab	2	\tab	Binary alphabet : 0 or 1	\cr
    2	\tab	4	\tab	DNA alphabet without ambiguity : A, C, G or T	\cr
    3	\tab	8	\tab	DNA alphabet with ambiguities : A, C, G, T, N, '-'	\cr
    4	\tab	16	\tab	DNA alphabet with all IUPAC ambiguity codes	\cr
    5	\tab	32	\tab	Amino-acid sequence and ambiguity (23 letters)	\cr
    6	\tab	64	\tab	Amino-acids with case (23x2)	\cr
    7	\tab	128	\tab	\strong{Not implemented}	\cr
    8	\tab	256	\tab	Full ASCII alphabet	\cr
  }
  
  And memory saving according to full ASCII can be avaluated by 8 / \code{nbits}.
}

\value{
  A \code{\linkS4class{rawFasta}} compliant object, whose actual implementation differs according to \code{nbits}. This can be an object of class \code{\linkS4class{rawFasta.8}}, \code{\linkS4class{rawFasta.sl}} or \code{\linkS4class{rawFasta.dl}}.
}

\note{
  The file to be parsed is supposed to be a full-text uncompressed ASCII file, begin with a comment line (from a '>' sign to the line break), and contain sequence letters on a random amount of fixed-width lines (width will be guessed from the first sequence line). The line breaking can be either '\\n', '\\r', '\\n\\r' or '\\r\\n', the file ending can hold an arbitrary amount of line breakings.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{extract}}, \code{\link{rawFasta-class}}
}

\examples{
  ### Typical DNA FASTA file
  fastaFile <- system.file("extdata/dna.fa", package="rawFasta")
  print(scan(fastaFile, what="", sep="\n"))
  
  # DNA, full ASCII
  f8 <- rawFasta(fastaFile, nbits=8)
  print(f8)
  
  # DNA, optimised
  f3 <- rawFasta(fastaFile, alpha="ACGTN-")
  print(f3)
  
  
  ### FASTA file, 2 letters
  fastaFile <- system.file("extdata/binary.fa", package="rawFasta")
  print(scan(fastaFile, what="", sep="\n"))
  
  # Binary, full ASCII
  f8 <- rawFasta(fastaFile, nbits=8)
  print(f8)
  
  # Binary, optimised
  f1 <- rawFasta(fastaFile, alpha="01")
  print(f1)
}
