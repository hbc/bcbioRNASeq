# Load up the featureCounts aligned counts matrix.
# Use the genes argument to dynamically resize the matrix. This is necessary
# when slotting this data into assays along with pseudoaligned counts.
.featureCounts <-
    function(projectDir, samples, genes = NULL) {
        assert(
            isADirectory(projectDir),
            isCharacter(samples),
            isCharacter(genes, nullOK = TRUE)
        )
        message("Importing aligned counts from featureCounts.")
        counts <- import(file = file.path(projectDir, "combined.counts"))
        assert(is.matrix(counts))
        colnames(counts) <- makeNames(colnames(counts))
        assert(isSubset(samples, colnames(counts)))
        counts <- counts[, samples, drop = FALSE]
        # We have to coerce to data frame to resize the matrix by row names.
        # You can't do this on a matrix directly, otherwise you'll hit a
        # "subscript out of bounds" error.
        if (is.character(genes)) {
            counts <- as.data.frame(counts)
            counts <- counts[genes, , drop = FALSE]
            rownames(counts) <- genes
            counts <- as.matrix(counts)
        }
        mode(counts) <- "integer"
        counts
    }
