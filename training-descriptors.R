library(protr)
library(purrr)

setwd(dir ="D:\\projects\\cbh21-protein-solubility-challenge\\data\\training\\crystal_structs")
pdbs <- list.files(pattern = ".pdb")
pdbs_seq <- map(head(pdbs), readPDB)
peptides <- modify(pdbs_seq, ~ .x $chain_A)

aac <- map(peptides, extractAAC)
paac <- map(peptides, extractPAAC)
apaac <- map(peptides, extractAPAAC)
mo_bro <- map(peptides, extractMoreauBroto)
moran <- map(peptides, extractMoran)
geary <- map(peptides, extractGeary)
ctdc <- map(peptides, extractCTDC)
ctdt <- map(peptides, extractCTDT)
ctdd <- map(peptides, extractCTDT)
ctriad <- map(peptides, extractCTriad)
socn <- map(peptides, extractSOCN)
qso <- map(peptides, extractSOCN)









