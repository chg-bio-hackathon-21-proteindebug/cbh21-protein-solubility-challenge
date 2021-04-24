library(protr)
library(purrr)

setwd(dir ="D:\\projects\\cbh21-protein-solubility-challenge\\data\\training\\crystal_structs")
pdbs <- list.files(pattern = ".pdb")
pdbs_seq <- map(head(pdbs), readPDB)

peptide <- pdbs_seq[[1]]$chain_A


pdbaac <- extractAAC(peptide)
accpaac <- extractPAAC(peptide)
apaac <- extractAPAAC(peptide)
mo_bro <- extractMoreauBroto(peptide)
moran <- extractMoran(peptide)
geary <- extractGeary(peptide)
ctdc <- extractCTDC(peptide)
ctdt <- extractCTDT(peptide)
ctdd <- extractCTDD(peptide)
ctriad <- extractCTriad(peptide)
socn <- extractSOCN(peptide)
qso <- extractQSO(peptide)



