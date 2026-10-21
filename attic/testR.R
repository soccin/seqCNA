# This will only work if you call this script with SDIR set properly
# like from a minimal bash wrapper script

SDIR=Sys.getenv("SDIR")
source(file.path(SDIR,"include/misc.R"))
library(seqDNAcopy)
