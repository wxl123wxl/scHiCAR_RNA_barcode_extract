a = read.csv("sciHiCARM_barcode.csv2", header=T)
a[1]
matrix_length=matrix(nrow=56*96*96*33, ncol=4)
matrix_length[,1]=seq(0,56*96*96*33-1)

b=a[matrix_length[,1]%%56+1,4]
c=a[floor(matrix_length[,1]%%(96*56)/56)+1,2]
d=a[floor(matrix_length[,1]%%(96*96*56)/(96*56))+1,3]

b2=as.character(b)
c2=as.character(c)
d2=as.character(d)
f2=cbind(b2,c2,d2)
f3=f2[1:(96*96*56),]
write.table(f3,file="barcode2.txt",sep="", quote=F, row.names = FALSE, col.names = FALSE)

