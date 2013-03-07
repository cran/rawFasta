# Extraction methods
# Author : Sylvain Mareschal <maressyl at gmail.com>
# License : GPL 3 http://www.gnu.org/licenses/gpl.html

setGeneric(
	name = "extract",
	def = function(object, start, end, check=TRUE, ...) {
		# Checks
		if(isTRUE(check)) {
			# Conversion
			start <- as.integer(start)
			end <- as.integer(end)
			
			# Negative positions
			if(start <= 0L) start <- object@size + start
			if(end <= 0L)   end <- object@size + end
			
			# Checks
			if(length(start) != 1L) stop(call.=FALSE, "'start' must be a single value")
			if(length(end) != 1L)   stop(call.=FALSE, "'end' must be a single value")
			if(is.na(start))        stop(call.=FALSE, "'start' must not be NA")
			if(is.na(end))          stop(call.=FALSE, "'end' must not be NA")
			if(end > object@size)   stop(call.=FALSE, "'end' must be lesser or equal to 'object' size")
			if(start > end)         stop(call.=FALSE, "'end' must be greater or equal to 'start'")
		}
		
		standardGeneric("extract")
	}
)

setMethod(
	f = "extract",
	signature = signature(object="rawFasta.8"),
	definition = function(object, start, end, check=TRUE, ...) {
		return(rawToChar(get("pack8", envir=object@content)[start:end], multiple=TRUE))		
	}
)

setMethod(
	f = "extract",
	signature = signature(object="rawFasta.sl"),
	definition = function(object, start, end, check=TRUE, ...) {
		indexes <- .rawFasta.extract(nbits=object@nbits, contentEnv=object@content, table=object@table, start=start, end=end)
		levels(indexes) <- object@alpha
		class(indexes) <- "factor"
		return(indexes)		
	}
)

setMethod(
	f = "extract",
	signature = signature(object="rawFasta.dl"),
	definition = function(object, start, end, check=TRUE, ...) {
		layer1 <- .rawFasta.extract(nbits=object@nbits[1], contentEnv=object@content, table=object@table[[1]], start=start, end=end)
		layer2 <- .rawFasta.extract(nbits=object@nbits[2], contentEnv=object@content, table=object@table[[2]], start=start, end=end)
		indexes <- (layer2-1L)*as.integer(2^object@nbits[1]) + layer1
		levels(indexes) <- object@alpha
		class(indexes) <- "factor"
		return(indexes)		
	}
)

# Extracts alphabet indexes from an environment (unique layer)
.rawFasta.extract <- function(nbits, contentEnv, table, start, end) {
	# Vector to return
	out <- integer(end-start+1L)
	
	# Byte holding 'start' and 'end' positions
	startByte <- ceiling(start / (8/nbits))
	endByte <- ceiling(end / (8/nbits))

	if(startByte == endByte) {
		
		### In-byte extraction
		
		# Expand concerned byte
		byte <- get(sprintf("pack%i", nbits), envir=contentEnv)[ startByte ]
		expandedByte <- table[, as.integer(byte)+1L ]
		
		# Bits to extract
		rangeStart <- start %% (8/nbits)
		rangeEnd <- end %% (8/nbits)
		if(rangeStart == 0) rangeStart <- (8/nbits)
		if(rangeEnd == 0) rangeEnd <- (8/nbits)
		
		# Extract
		out <- expandedByte[ rangeStart:rangeEnd ]
	} else {
		
		### Multi-byte extraction (before + [core] + after)

		# Before the core
		beforeCore <- (8/nbits) - (start-1L) %% (8/nbits)
		if(beforeCore == (8/nbits)) beforeCore <- 0L
		if(beforeCore != 0) {
			# Splitted byte
			previousByte <- get(sprintf("pack%i", nbits), envir=contentEnv)[ startByte ]
			expandedByte <- table[, as.integer(previousByte)+1L ]
			out[ 1:beforeCore ] <- tail(expandedByte, beforeCore)
		}

		# After the core
		afterCore <- end %% (8/nbits)
		if(afterCore != 0) {
			# Splitted byte
			nextByte <- get(sprintf("pack%i", nbits), envir=contentEnv)[ endByte ]
			expandedByte <- table[, as.integer(nextByte)+1L ]
			out[ length(out) - ((afterCore-1):0) ] <- head(expandedByte, afterCore)
		}

		# Core (if any)
		if(length(out) > beforeCore + afterCore) {
			coreIndexes <- (beforeCore + 1) : (length(out) - afterCore)
			coreBytes <- (ceiling((start-1) / (8/nbits)) + 1L) : (ceiling((end+1) / (8/nbits)) - 1L)
			out[ coreIndexes ] <- as.vector(table[, as.integer(get(sprintf("pack%i", nbits), envir=contentEnv)[ coreBytes ])+1L ])
		}
	}
	
	return(out)	
}

