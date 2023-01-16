if (!require(archive)) install.packages("archive")
library(archive)

if (!require(readr)) install.packages("readr")
library(readr)

importedData = read_csv(archive_read("data/grade10math.tar.gz", file=1))
