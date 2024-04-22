# this script is used to formated the .info file correctly.
# Sicne the imputed file(.pvar) from Sanger has different INFO value from the Minchigan one;
# and the format cannot be adjust simply by linux command.
# So, this script will get rid of the TYPE column for some of the SNPs and make sure 
# the INFO value is the last column.

# input: chr*_filtered.info
# for example
#              ID                TYPED                RPAF     ACAN INFO
# chr1:752566_A_G                TYPED RefPanelAF=0.819865  AC=2422    1
# chr1:752721_A_G  RefPanelAF=0.813505             AN=2838 0.985754   NA
# chr1:753405_A_C  RefPanelAF=0.852772             AN=2838 0.968644   NA

# output: chr*_filtered_formatted.info
# for example
#              ID                 RPAF     INFO
# chr1:752721_A_G  RefPanelAF=0.813505 0.985754
# chr1:753405_A_C  RefPanelAF=0.852772 0.968644
# chr1:753541_A_G  RefPanelAF=0.151355 0.957479

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

arguments <- commandArgs(T)
infofile <- arguments[1]
outputfile <- arguments[2]

# read in .info file
message("Reading in ", infofile)
df <- fread(infofile, data.table=F, fill=TRUE)
head(df)
# subset the SNPs with TYPED
type_index <- grep("^T", df$TYPED, ignore.case=TRUE)
type_table <- df[type_index,]
type_table_formated <- type_table[,c(1,3,5)]
colnames(type_table_formated) = c("ID", "RPAF", "INFO")
type_table_formated[c('RPAFname', 'RPAF')] <- str_split_fixed(type_table_formated$RPAF, '=', 2)
type_table_formated = type_table_formated[,c("ID", "RPAF", "INFO")]
head(type_table_formated)


# subset the SNPs without TYPED
nontype_index <- grep("^Ref", df$TYPED, ignore.case = T)
nontype_table <- df[nontype_index,]
nontype_table_formated <- nontype_table[,c(1,2,4)]
colnames(nontype_table_formated) = c("ID", "RPAF", "INFO")
nontype_table_formated[c('RPAFname', 'RPAF')] <- str_split_fixed(nontype_table_formated$RPAF, '=', 2)
nontype_table_formated = nontype_table_formated[,c("ID", "RPAF", "INFO")]
head(nontype_table_formated)


# combination
df_formated <- rbind(nontype_table_formated, type_table_formated)


write.table(df_formated,file=outputfile, sep = " ", quote = F)
