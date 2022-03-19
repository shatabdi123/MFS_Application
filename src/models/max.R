# ----------------
# Filename: max.R
# ----------------
# Fetch command line arguments
#install.packages("protr")
myArgs <- commandArgs(trailingOnly = TRUE)
#print(myArgs)
peptide_list =list()

peptide_list=append(peptide_list,as.character(myArgs))
len = nchar(as.character(myArgs)) -1


names(peptide_list) <- 'protein_seq'

.libPaths("C:\\Users\\Shatabdi\\Documents\\R\\win-library\\4.0")

library(protr)

peptide_list <- peptide_list[(sapply(peptide_list, protcheck))]

a1 = t(sapply(peptide_list, extractTC))
colnames(a1) <- paste("TC", colnames(a1), sep = "_")

a2 = t(sapply(peptide_list, extractCTriad))
colnames(a2) <- paste("CTriad", colnames(a2), sep = "_")

a3 = t(sapply(peptide_list, extractCTDD))
colnames(a3) <- paste("CTDD", colnames(a3), sep = "_")

a4 = t(sapply(peptide_list, extractPAAC,lambda=len))
colnames(a4) <- paste("PAAC", colnames(a4), sep = "_")

x1=cbind(a1,a2)
x1=cbind(x1,a3)
x1=cbind(x1,a4)

write.table(x1,"C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/GFF_File_Extractor/trying_protr.txt",quote=FALSE,sep="\t")

output = peptide_list

# Write the result to output stream
cat(as.character(output))
#cat(myArgs)
####################################
