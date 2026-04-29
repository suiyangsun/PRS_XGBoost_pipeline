'
Notes:
    1. Read data from stdin and output results to stdout.
    2. Normolize column by scale.
    3. Add on more column with header as scale.

Usage:
    Scale.R -c column_name [-t title]

Options:
    -c string   Column name to scale
    -t string   Output name.
    -h --help   Show this screen.
    --version   Show version.
' -> doc

## Auto-detect and install needed packages.
## Written by someone
options (warn = -1)
.libPaths(c("rpackage_score",.libPaths()))
#if (!require("pacman")) install.packages("pacman")
suppressMessages(library(pacman))
pacman::p_load(docopt,data.table)

opts <- docopt(doc,version="V1_2023.11.02")


title  <- if(is.null(opts$t)) 'scale' else opts$t
df = read.table(file("stdin"),header = T,check.names=F)

df[[title]] <- scale(df[opts$c])

write.table(df ,"",row.names=F,col.names=T,quote=F)
