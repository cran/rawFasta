# S4 class definitions
# Author : Sylvain Mareschal <maressyl at gmail.com>
# License : GPL 3 http://www.gnu.org/licenses/gpl.html

# Virtual interface
setClass(
	Class = "rawFasta",
	contains = "VIRTUAL",
	representation = list(
		nbits = "integer",		# Bits per letter, "8", "4" or "21" (2+1)
		header = "character",		# Comment line of the FASTA file, without ">"
		source = "character",		# Path and file of the parsed FASTA file
		content = "environment",	# Environment holding content variable(s)
		size = "integer"		# Length of the sequence
	)
)

# 8 bits special case
setClass(
	Class = "rawFasta.8",
	contains = "rawFasta"
)

# Single layer
setClass(
	Class = "rawFasta.sl",
	contains = "rawFasta",
	representation = list(
		alpha = "character",
		table = "matrix"
	)
)

# Double layer
setClass(
	Class = "rawFasta.dl",
	contains = "rawFasta",
	representation = list(
		alpha = "character",
		table = "list"
	)
)



# Show method
setMethod(
	f = "show",
	signature = signature(object="rawFasta"),
	definition = function(object) {
		# Sequence storage size
		seqSize <- ceiling(sum(object@size / (8/object@nbits)))
		units <- c("o", "kio", "Mio", "Gio", "Tio")
		u <- 1L
		while(seqSize > 1024) {
			seqSize <- seqSize / 1024
			u <- u +1L
		}
		
		# Print
		cat("An object of class \"rawFasta\"\n")
		cat("Source file   : ", object@source, "\n", sep="")
		cat("Sequence size : ", format(object@size, big.mark=" "), " (", round(seqSize,1), " ", units[u], ")\n", sep="")

		if(is(object, "rawFasta.8")) { cat("Alphabet      : ASCII (8 bits)\n", sep="")
		} else                       { cat("Alphabet      : [", object@alpha, "] (", sum(object@nbits), " bit", if(!identical(object@nbits, 1L)) { "s" }, ")\n", sep="")
		}
		
		cat("Methods       : extract, show\n")
	}
)

