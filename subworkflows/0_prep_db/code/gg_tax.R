# adapted from https://mothur.org/blog/2014/greengenes-v13_8_99-reference-files/

tax <- scan(file = "data/gg/gg_13_8_99.gg.tax", sep = "\t", what = "", quiet = T)
nlines <- length(tax)
ids <- tax[1:nlines %% 2 == 1] # the sequence names will be in the odd slots of the vector
tax.strings <- tax[1:nlines %% 2 == 0] # the taxonomy strings will be in the even slots of the vector

# tax.strings <- paste0(tax.strings, ";")      #make sure every sequence ends in a semi-colon
tax.strings <- gsub("; ", ";", tax.strings) # remove the spaces from between the taxonomic levels
tax.strings <- gsub(".__;", "", tax.strings) # remove the unclassified taxonomic levels (e.g. "c__")
tax.strings <- gsub(" ", "_", tax.strings) # if a taxonomic name has a space in it, replace it with an underscore
# tax.strings <- paste0("Root;", tax.strings)  #make every line start with the Root designation
tax.strings <- gsub(".__", "", tax.strings) # this line added by KLS: remove taxon level designations
updated.tax <- cbind(ids, tax.strings)
updated.tax <- updated.tax[order(as.numeric(updated.tax[, "ids"])), ] # order taxonomy file numerically
write.table(updated.tax, file = "data/gg/gg.tax", row.names = F, col.names = F, quote = F, sep = "\t")
