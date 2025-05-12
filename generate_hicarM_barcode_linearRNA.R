a = read.csv("sciHiCARM_barcode.csv", header=T)
a[1]
matrix_length=matrix(nrow=48*96*96*33, ncol=4)
matrix_length[,1]=seq(0,48*96*96*33-1)

b=a[matrix_length[,1]%%48+1,4]
c=a[floor(matrix_length[,1]%%(96*48)/48)+1,2]
d=a[floor(matrix_length[,1]%%(96*96*48)/(96*48))+1,3]

b2=as.character(b)
c2=as.character(c)
d2=as.character(d)
f2=cbind(b2,c2,d2)
f3=f2[1:(96*96*48),]
write.table(f3,file="barcode2.txt",sep="", quote=F, row.names = FALSE, col.names = FALSE)

