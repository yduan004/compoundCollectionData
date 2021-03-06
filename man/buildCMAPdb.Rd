% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manageCMAP.R
\name{buildCMAPdb}
\alias{buildCMAPdb}
\title{Build CMAP Database}
\usage{
buildCMAPdb(dest_dir = ".", rerun)
}
\arguments{
\item{dest_dir}{character(1), destination direcotry under which the result 
SQLite database named as 'cmap.db' stored}

\item{rerun}{TRUE or FALSE indicating whether to rerun the function}
}
\value{
cmap.db
}
\description{
This function builds the 'cmap.db' that contains id mappings of cmap names
to PubChem/DrugBank IDs as well as compound structure information.
}
\details{
For about 2/3 of the CMAP drungs, one can obtain their PubChem/DrugBank IDs from 
the DMAP site here: http://bio.informatics.iupui.edu/cmaps. 
The SMILES strings for CMAP entries were obtained from ChemBank. Compounds 
were matched by names using the 'stringdist' library where cmap_name from 
CMAP were mapped to the closest name in ChemBank.
}
\examples{
library(ChemmineR)
## Query database
# conn <- initDb("/inst/scripts/cmap.db")
# results <- getAllCompoundIds(conn)
# sdfset <- getCompounds(conn, results, keepOrder=TRUE)
# sdfset
# as.data.frame(datablock2ma(datablock(sdfset)))[1:4,]
# myfeat <- listFeatures(conn)
# feat <- getCompoundFeatures(conn, results, myfeat)
# feat[1:4,]
}
