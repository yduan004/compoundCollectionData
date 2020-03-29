##########################################
## Create SQLite database and SDF files ##
##########################################
## Author: Yuzhu Duan
## Last update: 25-March-2020

# The SQLite database containing id mappings and drug annotation tables in
# DrugAge, DrugBank, CMAP2 and LINCS databases
getwd() # under compoundCollection directory

###### DrugBank ######

# get drugbank annotation table
## go to https://www.drugbank.ca/releases/latest download the drugbank xml database
## as 'drugbank_5.1.5.xml.zip' file to `inst/script` directory
unzip("data-raw/drugbank_5.1.5.xml.zip", exdir="inst/script")
file.rename("inst/script/full database.xml", "inst/script/drugbank_5.1.5.xml")
library(drugbankR)
drugbank_dataframe <- dbxml2df(xmlfile="data-raw/drugbank_5.1.5.xml", version="5.1.5")
df2SQLite(dbdf=drugbank_dataframe, version="5.1.5")
library(filesstrings)
file.move("drugbank_5.1.5.db", "inst/script")
library(RSQLite)
dbpath <- "inst/script/drugbank_5.1.5.db"
conn <- dbConnect(SQLite(), dbpath)
dbListTables(conn)
dbdf <- dbGetQuery(conn, 'SELECT * FROM dbdf')
pryr::object_size(dbdf)
colnames(dbdf)

# get drugbank id to ChEMBL id mapping table
download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz",
              "inst/script/chembl2drugbank.txt.gz")
library(readr)
chem2db <- read_tsv("inst/script/chembl2drugbank.txt.gz") # 7290 X 2
colnames(chem2db) <- c("chembl_id", "drugbank_id")
length(intersect(dbdf$`drugbank-id`, chem2db$drugbank_id))
## 7290 / 13475 (54.1%) drugs in drugbank have chembl id
library(dplyr)
dbdf <- dbdf %>% rename("drugbank_id"="drugbank-id")

# download drugbank sdf file
## downlaod 'drugbank_all_structures.sdf.zip' file from https://www.drugbank.ca/releases/latest#structures
## under the 'inst/script' directory
unzip("inst/script/drugbank_all_structures.sdf.zip", exdir="inst/script")
file.rename("inst/script/structures.sdf", "inst/script/drugbank_5.1.5.sdf")
db_stru <- read.SDFset("inst/script/drugbank_5.1.5.sdf") # 10,695

# Try to write drugbank id to Molecule_Name
dbids <- sapply(seq_along(db_stru), function(i){
    unlist(datablock(db_stru[[i]]))[1]
})
cid(db_stru) <- dbids
db_valid <- db_stru[validSDF(db_stru)] # 10,569
write.SDF(db_valid, "inst/extdata/drugbank_5.1.5.sdf", cid=TRUE)

######## DrugAge ########

# Generate drugage.db by using buildDrugAgeDB function from drugTargetInteractions package
library(drugTargetInteractions)
chembldb <- "/bigdata/girkelab/shared/lcshared/chemoinformatics/compoundDBs/chembl_24/chembl_24_sqlite/chembl_24.db"
cmp_ids.rds <- "/rhome/tgirke/Projects/longevity/LongevityConsortium/drugTargetInteraction_Analysis/results/cmp_ids.rds"
queryBy <- list(molType="cmp", idType="chembl_id", ids=c("CHEMBL1233058", "CHEMBL1200916", "CHEMBL437765"))
qresult <- drugTargetAnnot(dbpath=chembldb, queryBy, cmpid_file=cmp_ids.rds)
buildDrugAgeDB(dest_path="inst/extdata/drugage.db") # run only once
qresult_da <- drugAgeAnnot(dtAnnotRes=qresult, chemblCol="chembl_id", drugage_dbpath="inst/extdata/drugage.db")

# Add DrugBank, CMAP and LINCS compounds annotations to drugage.db and rename as compoundCollection.db
file.rename("inst/extdata/drugage.db", "inst/extdata/compoundCollection.db")
cc_path <- "inst/extdata/compoundCollection.db"
cc_conn <- dbConnect(SQLite(), cc_path)
dbListTables(cc_conn)
id_mapping <- dbGetQuery(cc_conn, "SELECT * FROM id_mapping")
id_mapping2 <- merge(chem2db, id_mapping, by="chembl_id", all.x=TRUE, all.y=TRUE) # 7442 X 3
dbWriteTable(cc_conn, "id_mapping", id_mapping2, overwrite=TRUE)
dbWriteTable(cc_conn, "DrugBankAnnot", dbdf)
dbListTables(cc_conn)
dbDisconnect(cc_conn)

# Get drugage sdf file
da_annot <- dbGetQuery(cc_conn, "SELECT * FROM drugAgeAnnot")
ida2pubchem <- na.omit(da_annot[,c("drugage_id", "pubchem_cid")])
pid_unique <- gsub(",.*", "", ida2pubchem$pubchem_cid)
library(ChemmineR)
da_cmp <- getIds(as.numeric(pid_unique))
write.SDF(da_cmp, "inst/extdata/drugage.sdf")

####### CMAP2 #########

download.file("https://portals.broadinstitute.org/cmap/cmap_instances_02.xls",
              "inst/script/cmap_instances_02.xls")
## Note, this file required some cleaning in LibreOffice (Excel would work for this too).
## After this it was saved as tab delimited txt file named [cmap_instances_02.txt]
download.file("http://cluster.hpcc.ucr.edu/~tgirke/projects/longevity/cmap/data/cmap_instances_02.txt",
              "inst/script/cmap_instances_02.txt")
cmap_inst <- read.delim("inst/script/cmap_instances_02.txt", check.names=FALSE)
## only contains compound names and catelog number

# The cmap.db has been built by the longevityDrugs package
## cmap.db is a SQLite structure database for drugs from CMAP02. It is built by the
## "buildCMAPdb" funciton in the longevityDrugs https://github.com/tgirke/longevityDrugs package
devtools::install_github("tgirke/longevityDrugs")
library(longevityDrugs)
library(ChemmineR); library(RSQLite)
mypath <- system.file("extdata", "cmap.db", package="longevityDrugs")
conn <- dbConnect(SQLite(), mypath) # or conn <- initDb(mypath)
results <- getAllCompoundIds(conn)
sdfset <- getCompounds(conn, results, keepOrder=TRUE)
sdfset # An instance of "SDFset" with 1309 molecules
plot(sdfset[1:4], print=FALSE)

myfeat <- listFeatures(conn)
feat <- getCompoundFeatures(conn, results, myfeat)
feat[1:4,]

cmp_tb <- dbGetQuery(conn, 'SELECT * FROM compounds')
cmapAnnot <- as.data.frame(datablock2ma(datablock(sdfset))) # 1309 X 27
cmap_id <- paste0("CMAP", sprintf("%04d", as.integer(rownames(cmapAnnot))))
cmapAnnot2 <- data.frame(cmap_id=cmap_id, cmapAnnot, stringsAsFactors=FALSE)
# write CMAP SDF file
cid(sdfset) <- cmap_id
write.SDF(sdfset, file="inst/extdata/cmap.sdf", cid=TRUE)
cmap_sdf <- read.SDFset("inst/extdata/cmap.sdf")

# get pubchem cid to ChEMBL id mapping table
download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz",
              "inst/script/chembl2pubchemcid.txt.gz")
library(readr)
chem2pub <- read_tsv("inst/script/chembl2pubchemcid.txt.gz") # 1837162 X 2
colnames(chem2pub) <- c("chembl_id", "pubchem_cid")
chem2pub$pubchem_cid <- as.character(chem2pub$pubchem_cid)
# get cmap_id (created internal) to pubchem cid mappings
cmap2pub <- cmapAnnot2[,c("cmap_id", "PUBCHEM_ID")]
sum(cmap2pub$PUBCHEM_ID=="NA") # 443
## remove "NA"
cmap2pub <- cmap2pub[cmap2pub$PUBCHEM_ID != "NA",] # 866 X 2
dup <- cmap2pub[cmap2pub$PUBCHEM_ID %in% cmap2pub$PUBCHEM_ID[duplicated(cmap2pub$PUBCHEM_ID)],]
cmapAnnot2[cmapAnnot2$cmap_id %in% dup$cmap_id,]
### two cmap drugs with different names have the same pubchem cid, there are three this situation
# cmap_id to chembl id mapping
cmap2pub$PUBCHEM_ID <- gsub("CID", "", cmap2pub$PUBCHEM_ID)
library(dplyr)
cmap2chem <- cmap2pub %>% left_join(chem2pub, by=c("PUBCHEM_ID"="pubchem_cid"))
cmap2chem <- unique(na.omit(cmap2chem[,c("cmap_id", "chembl_id")])) # 800 X 2

cc_path <- "inst/extdata/compoundCollection.db"
cc_conn <- dbConnect(SQLite(), cc_path)
dbListTables(cc_conn)
id_mapping <- dbGetQuery(cc_conn, "SELECT * FROM id_mapping")
id_mapping2 <- merge(id_mapping, cmap2chem, by="chembl_id", all.x=TRUE, all.y=TRUE) # 7599 X 4
dbWriteTable(cc_conn, "id_mapping", id_mapping2, overwrite=TRUE)
dbWriteTable(cc_conn, "cmapAnnot", cmapAnnot2)
dbListTables(cc_conn)
dbDisconnect(cc_conn)

####### LINCS #######

# get LINCS annotation table
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz"
download.file(url, "./inst/script/GSE92742_Broad_LINCS_pert_info.txt.gz")
pertIDs <- read.delim("inst/script/GSE92742_Broad_LINCS_pert_info.txt.gz")
pertIDs <- pertIDs[pertIDs$pert_type=="trt_cp", ] # 20,413 X 8
lincsAnnot <- pertIDs
colnames(lincsAnnot)[1] <- "lincs_id"

# get lincs internal id (e.g. BRD-A00100033) to chembl id mapping via inchi_key
chembldb <- "~/insync/project/ChEMBL_data/chembl_25/chembl_25_sqlite/chembl_25.db"
library(RSQLite)
conn <- dbConnect(SQLite(), chembldb)
cmp_str <- dbGetQuery(conn, 'SELECT chembl_id, molregno, standard_inchi_key
                             FROM compound_structures, chembl_id_lookup
                             WHERE compound_structures.molregno = chembl_id_lookup.entity_id') # 3,099,738 X 3
lincs_inchi <- pertIDs[,c("pert_id", "inchi_key")]
lincs_inchi <- lincs_inchi[lincs_inchi$inchi_key != "-666",] # 20350 X 2
length(unique(lincs_inchi$inchi_key)) # 20307
library(dplyr)
lincs_inchi2 <- lincs_inchi %>% left_join(cmp_str[,c("chembl_id", "standard_inchi_key")],
                                          by=c("inchi_key"="standard_inchi_key"))
lincs2chem <- unique(na.omit(lincs_inchi2[,c("pert_id", "chembl_id")])) # 10,930 X 2
colnames(lincs2chem) <- c("lincs_id", "chembl_id")
length(unique(lincs2chem$lincs_id)) # 5774 unique lincs internal id
length(unique(lincs2chem$chembl_id)) # 10,878 unique chembl id

# write annotation table and id_mapping table to SQLite db
cc_path <- "inst/extdata/compoundCollection.db"
cc_conn <- dbConnect(SQLite(), cc_path)
dbListTables(cc_conn)
id_mapping <- dbGetQuery(cc_conn, "SELECT * FROM id_mapping")
id_mapping2 <- merge(id_mapping, lincs2chem, by="chembl_id", all.x=TRUE, all.y=TRUE) # 17,047 X 5
dbWriteTable(cc_conn, "id_mapping", id_mapping2, overwrite=TRUE)
dbWriteTable(cc_conn, "lincsAnnot", lincsAnnot)
dbListTables(cc_conn)
dbDisconnect(cc_conn)

# write LINCS sdf file
library(ChemmineR)
library(ChemmineOB) # make sure openbabel module is loaded!!
pertIDs <- pertIDs[pertIDs$canonical_smiles!=-666, ]
smiset <- as.character(pertIDs$canonical_smiles)
pertids <- as.character(pertIDs$pert_id)
names(smiset) <- pertids
lincs_sdfset <- smiles2sdf(smiset)
valid <- validSDF(lincs_sdfset); lincs_sdfset <- lincs_sdfset[valid]
## note, original smiset contained 20350 and final valid SDFset 20333
saveRDS(lincs_sdfset, "inst/script/lincs_sdfset.rds")
lincs_sdfset <- readRDS("inst/script/lincs_sdfset.rds")
write.SDF(lincs_sdfset, file="inst/extdata/lincs.sdf")
