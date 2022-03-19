# ----------------
# Filename: max.R
# ----------------
# Fetch command line arguments
#install.packages("protr")

DNAArgs <- commandArgs(trailingOnly = TRUE)
#print(myArgs)
coding_list =list()

coding_list=append(coding_list,as.character(DNAArgs))
DNAlen = nchar(as.character(DNAArgs)) -1


names(coding_list) <- 'DNA_seq'

.libPaths("C:\\Users\\Shatabdi\\Documents\\R\\win-library\\4.0")

####################################
library("rDNAse")
a41 = t(sapply(coding_list, extrPseDNC))
colnames(a41) <- paste("PseDNC", colnames(a41), sep = "_")

a43 = t(sapply(coding_list, extrPseKNC, 3))
colnames(a43) <- paste("PseKNC_3", colnames(a43), sep = "_")

x2=cbind(a41,a43)

write.table(x2,"C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/GFF_File_Extractor/DNA.txt",quote=FALSE,sep="\t")

output2 = coding_list

# Write the result to output stream
cat(as.character(output2))