# Memory efficient handling of FASTA files, pure R implementation
# Author : Sylvain Mareschal <maressyl at gmail.com>
# License : GPL 3 http://www.gnu.org/licenses/gpl.html

# Initializes environment and 'nbits' related variables
.rawFasta.initialize <- function(nbits, contentSize, contentEnv, rawAlpha) {
	# Alpha check
	if(length(rawAlpha) > 2^sum(nbits)) stop("'alpha' can hold ", 2^sum(nbits), " characters at most on ", sum(nbits), " bits")
	
	# Some modes needs two layers
	layers <- 1:length(nbits)
	
	# Storage sizes (buffers have to be divisible by 8 for compactability)
	contentStorageSizes <- (ceiling(contentSize/8)*8) / (8/nbits)
	
	# Allocate storage vectors
	for(l in layers) assign(sprintf("pack%i", nbits[l]), raw(contentStorageSizes[l]/(8/nbits[l])), envir=contentEnv)
	
	# Storage expressions
	storeExpr <- parse(text=sprintf("pack%i[ range ] <- content[[l]]", nbits))
	
	# Bit base representation
	bits <- list()
	for(l in layers) bits[[l]] <- t(matrix(as.integer(intToBits(1:2^nbits[l]-1)), nrow=2^nbits[l], byrow=TRUE)[,1:nbits[l]])
	
	return(list(bits=bits, storeExpr=storeExpr, contentStorageSizes=contentStorageSizes))
}

# Stores a buffer (complete buffer or partial final buffer)
.rawFasta.store <- function(nbits, buffer, bufferLength, contentStorageSizes, i, contentEnv, storeExpr, rawAlpha, bits) {
	if(identical(nbits, 8L)) {
		# Use raw buffer in 8 bits mode
		content <- list(buffer)
	} else {
		# Compact the buffer
		content <- match(buffer, rawAlpha)
		if(any(is.na(content))) stop("Unexpected character(s) encountered : ", paste(unique(buffer[ is.na(content) ]), collapse=", "))
		if(length(nbits) == 1) {
			# Single layer
			content <- list(packBits(as.vector(bits[[1]][, content ])))
		} else {
			# Double layer
			content <- list(
				packBits(as.vector(bits[[1]][, (content - 1L) %% 2^nbits[1]  + 1L ])),
				packBits(as.vector(bits[[2]][, (content - 1L) %/% 2^nbits[1] + 1L ]))
			)
		}
	}
	
	layers <- 1:length(nbits)
	for(l in layers) {
		# Upper limit
		if(length(buffer) < bufferLength) {
			# Last buffer : untill the end
			upper <- contentStorageSizes[l]
		} else {
			# Full buffers
			upper <- (i*bufferLength/(8/nbits[l]))
		}
	
		# Store the buffer content
		range <- ((i-1L)*bufferLength/(8/nbits[l])+1L) : upper
		eval(
			expr = do.call(substitute, list(storeExpr[[l]])),
			envir = contentEnv
		)
	}
}

# Resizes a filled partial final buffer
.rawFasta.resize <- function(buffer, rawAlpha, lastBufferLength) {
	# Resize buffer (divisible by 8 assure the compactability)
	lastBufferResize <- ceiling(lastBufferLength/8)*8
	buffer <- buffer[1:lastBufferResize]
	
	# Fill with last base letter
	if(lastBufferResize > lastBufferLength) {
		buffer[ (lastBufferLength+1) : lastBufferResize ] <- tail(rawAlpha, 1L)
	}
	
	return(buffer)
}

# Builds the correspondence table for future extractions
.rawFasta.table <- function(nbits, bits) {
	layers <- 1:length(nbits)
	intVersions <- list()
	for(l in layers) {
		# 256 bit combinations
		bitVersion <- t(matrix(as.integer(intToBits(1:2^8-1)), nrow=2^8, byrow=TRUE)[,1:8])

		# Base index sequence corresponding to each
		intVersions[[l]] <- matrix(as.integer(NA), nrow=(8/nbits[l]), ncol=2^8)
		for(r in 1:nrow(intVersions[[l]])) {
			for(k in 1:ncol(intVersions[[l]])) {
				charBits <- bitVersion[ ((r-1)*nbits[l]+1):(r*nbits[l]) , k ]
				charIndex <- which(apply(bits[[l]] == charBits, 2, all))
				intVersions[[l]][r,k] <- charIndex
			}
		}
	}
	
	return(intVersions)
}

# Constructor
rawFasta <- compiler::cmpfun(
	function(
			fileName,
			maxBuffer = 8192L,
			alpha = "ACGTN-",
			nbits = NA
		) {
		### Assuming first line is always header
		### Assuming fixed common width for all sequence lines
		### Final octet completion is done with 'alpha' final character
		
		# nbits checks
		if(length(nbits) != 1) stop("'nbits' must be a single integer value")
		nbits <- as.integer(nbits)
		
		# Auto nbits
		if(identical(nbits, NA_integer_)) {
			nbits <- as.integer(ceiling(log(nchar(alpha), 2)))
			if(nbits == 7L) nbits <- 8L
		}
		
		# nbits decomposition
		if(!nbits %in% c(1L, 2L, 3L, 4L, 5L, 6L, 8L)) stop("'nbits' must be 1, 2, 3, 4, 5, 6 or 8")
		if(nbits == 3L)        { nbits <- c(2L, 1L)
		} else if(nbits == 5L) { nbits <- c(4L, 1L)
		} else if(nbits == 6L) { nbits <- c(4L, 2L)
		} 
		
		# Environment to store the content vector
		contentEnv <- new.env()
		
		# Opening file
		con <- file(fileName, open="rb")
		on.exit(close(con))
		
		# Header line
		header <- readLines(con, n=1L)
		if(substr(header,1,1) != ">") stop("'fileName' is not a regular FASTA file (no '>' header)")
		headerSize <- nchar(header)
		headerEnd <- seek(con)
		
		# Line breaking detection
		seek(con, where=-2L, origin="current")
		breaks <- readBin(con, what="raw", n=2L)
		if(all(breaks == charToRaw("\r\n")))    { lineBreaks <- 2L
		} else if(breaks[2] == charToRaw("\n")) { lineBreaks <- 1L
		} else                                  { stop("Unknown line breaking character pattern")
		}
		
		# First line (to guess the FASTA width)
		firstLine <- readLines(con, n=1L)
		lineLength <- nchar(firstLine)
		
		# Ending detection
		seek(con, where=-20L, origin="end")
		ending <- readBin(con, "raw", n=20L)
		i <- 0L
		while(ending[20L-i] %in% charToRaw("\r\n")) i <- i + 1L
		endingSize <- i
		
		# Sequence size estimation from file size
		contentSize <- as.integer(file.info(fileName)[1,'size']) - (headerSize + lineBreaks) - endingSize
		contentBreaks <- as.integer(floor(contentSize / (lineLength + lineBreaks)))
		contentSize <- contentSize - contentBreaks*lineBreaks
		
		# Sequence allocation
		rawAlpha <- charToRaw(alpha)
		ini <- .rawFasta.initialize(nbits=nbits, contentSize=contentSize, contentEnv=contentEnv, rawAlpha=rawAlpha)
		bits <- ini$bits
		storeExpr <- ini$storeExpr
		contentStorageSizes <- ini$contentStorageSizes
		
		# Buffer initialization (multiple of 8 lines, thus divisible by 1,2,4,8)
		bufferChunks <- floor(maxBuffer / (8*lineLength)) * 8
		if(bufferChunks == 0) stop("'maxBuffer' is too small to hold at least 8 lines of the file")
		bufferLength <- bufferChunks * lineLength
		buffer <- raw(bufferLength)
		bufferFilled <- 0L
		
		# Go back to sequence start
		seek(con, where=headerEnd, origin="start")
		
		# Looping on complete buffers (if file is large enough)
		i <- 0L
		if(floor(contentBreaks / bufferChunks) > 0) {
			bufferLoop <- 1:floor(contentBreaks / bufferChunks)
			bufferFillingLoop <- 1:bufferChunks
			for(i in bufferLoop) {
				# Bufferize line by line
				for(j in bufferFillingLoop) {
					buffer[ ((j-1L)*lineLength+1L) : (j*lineLength) ] <- readBin(con, what="raw", n=lineLength)
					readBin(con, what="raw", n=lineBreaks)
				}
			
				# Store buffer
				.rawFasta.store(nbits=nbits, buffer=buffer, bufferLength=bufferLength, contentStorageSizes=contentStorageSizes, i=i, contentEnv=contentEnv, storeExpr=storeExpr, rawAlpha=rawAlpha, bits=bits)
				}
		}
		
		# Last partial buffer
		lastBufferLength <- contentSize - i*bufferLength
		if(lastBufferLength > 0) {
			# Bufferize complete lines
			lastBufferLoops <- floor(lastBufferLength / lineLength)
			if(lastBufferLoops > 0) {
				for(j in 1:lastBufferLoops) {
					buffer[ ((j-1L)*lineLength+1L) : (j*lineLength) ] <- readBin(con, what="raw", n=lineLength)
					readBin(con, what="raw", n=lineBreaks)
				}
			}
			
			# Bufferize last partial line, if any
			lastBufferRemains <- lastBufferLength - lastBufferLoops * lineLength
			if(lastBufferRemains > 0) {
				buffer[ (j*lineLength) + (1:lastBufferRemains) ] <- readBin(con, what="raw", n=lastBufferRemains)
			}
			
			# Resize last buffer
			buffer <- .rawFasta.resize(buffer=buffer, rawAlpha=rawAlpha, lastBufferLength=lastBufferLength)
			
			# Store last buffer
			i <- i +1L
			.rawFasta.store(nbits=nbits, buffer=buffer, bufferLength=bufferLength, contentStorageSizes=contentStorageSizes, i=i, contentEnv=contentEnv, storeExpr=storeExpr, rawAlpha=rawAlpha, bits=bits)
		}
		
		# S4 Object
		if(length(nbits) > 1)  { object <- new("rawFasta.dl")
		} else if(nbits == 8L) { object <- new("rawFasta.8")
		} else                 { object <- new("rawFasta.sl")
		}
		object@nbits <- nbits
		object@header <- substr(header, 2, nchar(header))
		object@source <- normalizePath(fileName)
		object@content <- contentEnv
		object@size <- contentSize
		if(class(object) %in% c("rawFasta.sl", "rawFasta.dl")) object@alpha <- strsplit(alpha, split="")[[1]]
		if(class(object) == "rawFasta.sl") object@table <- .rawFasta.table(nbits=nbits, bits=bits)[[1]]
		if(class(object) == "rawFasta.dl") object@table <- .rawFasta.table(nbits=nbits, bits=bits)
		
		return(object)
	},
	options = list("optimize"=3)
)

