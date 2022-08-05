###########################################################
#
# An example launch script, which will process an
# arglist file, and launch cluster jobs on the DIDE
# cluster, each of which will call the batch file
# to run a single job.

# For first time setup:

first_time <- function() {
  install.packages("drat")
  drat:::add("mrc-ide")
  install.packages(c("didehpc", "context"))
}

library(didehpc)

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.template = "GeneralNodes",
        didehpc.cores = 1,
        didehpc.wholenode = FALSE)

# Set working directory to the network share.
# Below is specific to me (Z: is the arbitrary mapping on windows
# I have set up for the share. On mac/linux, set to
# "//wpia-hpc-hn/Fisher/bwa_mem_picard_align"

setwd("Z:/bwa_mem_picard_align/")

# Set up very basic cluster queue - no packages or sources needed
# as we're just using this as a batch runner

root <- "contexts"
didehpc::web_login()
ctx <- context::context_save(root)
queue <- didehpc::queue_didehpc(ctx)

lookup_file <- function(f) {
  files <- list.files(dirname(f), basename(f))
  files <- files[!grepl("md5", files)]
  file.path(dirname(f), files)
}

launch_arglist <- function(jobpath, batpath, argfile) {
  txt <- read.csv(argfile, sep = " ", header = FALSE,
                  col.names = c("R1", "R2", "pre"))
  for (i in seq_len(nrow(txt)))  {
    job <- txt[i, ]
    
    # Lookup files from the wildcards
    # Build into arg list
    # And launch the job.
    
    args <- c(jobpath, lookup_file(job$R1), lookup_file(job$R2), job$pre)
    queue$enqueue(system2(command = batpath, args = args))
  }
}

jobpath <- "\\\\wpia-hpc-hn\\Fisher\\bwa_mem_picard_align"
batpath <- paste0(jobpath, "\\bwa_mem_picard_align.bat")
argfile <- "arglist_run7.txt"
launch_arglist(jobpath, batpath, argfile)
