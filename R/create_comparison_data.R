create_comparison_data <- function(records, types, breaks){
    r <- nrow(records)
    record_pairs <- t(utils::combn(r, 2))
    rps1 <- record_pairs[, 1]
    rps2 <- record_pairs[, 1]
    FF <- ncol(records)

    comparisons <- matrix(0, nrow = ncol(record_pairs), ncol = (2 + FF))
    comparisons[, 1:2] <- record_pairs
    colnames(comparisons) <- c("i", "j", colnames(records))

    for(f in 1:FF){
        if(types[f] == "bi"){
            raw_comp <- records[rps1, f] != records[rps2, f]
            comparisons[, 2 + f] <- as.numeric(raw_comp) + 1
        }
        if(types[f] == "nu"){
            raw_comp <- abs(records[rps1, f] - records[rps2, f])
            temp_breaks <- unique(c(-Inf, breaks[[f]], Inf))
            comparisons[, 2 + f] <- cut(raw_comp, breaks = temp_breaks,
                                        labels = seq_len(length(temp_breaks) -
                                                             1))
        }
        if(types[f] == "lv"){
            raw_comp <- as.numeric(utils::adist(records[, f], records[, f]))/
                pmax(nchar(records[, f])[rps1], nchar(records[, f])[rps2])
            temp_breaks <- unique(c(-Inf, breaks[[f]], Inf))
            comparisons[, 2 + f] <- cut(raw_comp, breaks = temp_breaks,
                                        labels = seq_len(length(temp_breaks) -
                                                             1))
        }
    }

    return(as.data.frame(comparisons))
}
