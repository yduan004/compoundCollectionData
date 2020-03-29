#' Compound Annotation and Structure Datasets
#'
#' @name compoundCollectionData-package
#' @aliases compoundCollectionData-package compoundCollectionData
#' @docType package
#' @description This package contains annotation and 
#' structure datasets for compounds in DrugAge, DrugBank, CMAP02 and LINCS databases.
#' The actual datasets are stored in \code{AnnotationHub}
#' 
#' @details
#' The description of the 5 datasets in this package is as follows:
#' 
#' \strong{SQLite annotation database:}
#' 
#' It is a SQLite database storing compound annotation tables for DrugAge, DrugBank,
#' CMAP02 and LINCS, respectively. It also contains an ID mapping table of 
#' ChEMBL ID to IDs of individual databases. 
#' 
#' \strong{DrugAge SDF:}
#' 
#' It is an SDF (Structure-Data File) file storing molecular structures of DrugAge
#' comopunds. The source annotation file of DrugAge is downloaded from
#' \url{http://genomics.senescence.info/drugs/dataset.zip}. The extracted csv file
#' only contains drug names, without id mappings to external resources such as 
#' PubChem or ChEMBL. The extracted 'drugage.csv' file was processed by the 
#' \code{processDrugage} function from \pkg{drugTargetInteractions} package available
#' at \url{https://github.com/longevity-consortium/LC_Chemoinformatics/tree/master/Rpackages/drugTargetInteractions}.
#' The missing IDs were then added manually. The result DrugAge annotation table
#' with external ID mappings is available at https://bit.ly/3dCWKWo. This annotation
#' table along with the id-mapping table (ChEMBL ID to DrugAge internal id) were
#' stored in the compoundCollection SQLite database. The drug structures were
#' obtained from PubChem CIDs by \code{getIds} function from \pkg{ChemmineR} package,
#' then wrote to the 'drugage.sdf' file
#' 
#' \strong{DrugBank SDF:}
#' 
#' This SDF file stores structures of compounds in DrugBank \url{https://www.drugbank.ca/}.
#' The full DrugBank xml database was downloaded from \url{https://www.drugbank.ca/releases/latest}.
#' The most recent release version at this time is 5.1.5.  
#' The extracted xml file was procesed by the \pkg{drugbankR} package \url{https://github.com/yduan004/drugbankR}.
#' The result DrugBank annotation table was then stored in the compoundCollection 
#' SQLite database. The DrugBank to ChEMBL id mappings were obtained from 
#' UniChem at \url{ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz}.
#' The DrugBank SDF file was downloaded from \url{https://www.drugbank.ca/releases/latest#structures}
#' and made some validity checks and modifications via utilities in the \pkg{ChemmineR} package.
#' 
#' \strong{CMAP SDF:}
#' 
#' The CMAP compound annotation table and its structures were obtained from the
#' \pkg{longevityDrugs} package \url{https://github.com/tgirke/longevityDrugs}. 
#' Since this table only contains PubChem CID, the ChEMBL ids were added via
#' PubChem CID to ChEMBL id mappings from UniChem 
#' \url{ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz}.
#' The CMAP internal IDs were made for ChEMBL to CMAP id mappings.
#' 
#' \strong{LINCS SDF:}
#' 
#' The LINCS compound annotation table was downloaded from GEO
#' \url{ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz}
#' where only compounds were selected.
#' The LINCS ids were mapped to ChEMBL ids via inchi key. The LINCS compounds 
#' structures were obtained from PubChem CIDs via \code{getIDs} function from
#' \pkg{ChemmineR} package.
#' 
#' The R script of generating the above 5 datasets is available at the 
#' 'inst/scripts/make-data.R' file in this package.
#' 
#' @import AnnotationHub
#' @author
#' \itemize{
#'   \item Yuzhu Duan (yduan004@ucr.edu)
#'   \item Thomas Girke (thomas.girke@ucr.edu)
#' } 
#' @examples 
#' library(AnnotationHub)
#' \dontrun{
#' 
#'     ah <- AnnotationHub()
#'     
#'     ## Load compoundCollection annotation SQLite database
#'     query(ah, c("compoundCollectionData", "annot"))
#'     annot_path <- ah[[""]]
#'     library(RSQLite)
#'     conn <- dbConnect(SQLite(), annot_path)
#'     dbListTables(conn)
#'     dbDisconnect(conn)
#'     
#'     ## Load DrugAge SDF file
#'     query(ah, c("compoundCollectionData", "drugage"))
#'     da_path <- ah[[""]]
#'     da_sdfset <- ChemmineR::read.SDFset(da_path)
#'     
#'     ## Load DrugBank SDF file
#'     query(ah, c("compoundCollectionData", "drugbank"))
#'     db_path <- ah[[""]]
#'     db_sdfset <- ChemmineR::read.SDFset(db_path)
#'     
#'     ## Load CMAP SDF file
#'     query(ah, c("compoundCollectionData", "cmap"))
#'     cmap_path <- ah[[""]]
#'     cmap_sdfset <- ChemmineR::read.SDFset(cmap_path)
#'     
#'     ## Load LINCS SDF file
#'     query(ah, c("compoundCollectionData", "lincs"))
#'     lincs_path <- ah[[""]]
#'     lincs_sdfset <- ChemmineR::read.SDFset(lincs_path)
#' }
#' 
#' 

NULL