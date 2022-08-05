###########################################################
#
# An example launch script for raxml - which is just one
# line, so we don't need the arglist file. 

# For first time setup:

first_time <- function() {
  install.packages("drat")
  drat:::add("mrc-ide")
  install.packages(c("didehpc", "context"))
}

library(didehpc)

threads <- 16

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.template = "20Core",
        didehpc.cores = threads,
        didehpc.wholenode = FALSE)

# Set working directory to the network share.
# Below is specific to me (Z: is the arbitrary mapping on windows
# I have set up for the share. On mac/linux, it depends where the
# mount is - perhaps set to "//wpia-hpc-hn/Fisher/global" or
# on Mac it might be somewhere in "/Volumes"

setwd("/Volumes/Fisher/global/")

# Set up very basic cluster queue - no packages or sources needed
# as we're just using this as a batch runner

root <- "contexts"
didehpc::web_login()
ctx <- context::context_save(root)
queue <- didehpc::queue_didehpc(ctx)

jobpath <- "\\\\wpia-hpc-hn\\Fisher\\global"
batpath <- paste0(jobpath, "\\raxml_global.bat")

x <- queue$enqueue(system2(command = batpath, args = threads))
x$status()
