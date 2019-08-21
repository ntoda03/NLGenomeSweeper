
options(stringsAsFactors = FALSE)

m_opt = c(commandArgs()[3:length(commandArgs())])

infile=m_opt[1]
outfile=m_opt[2] 
gap=as.numeric(m_opt[3])

mytab=read.table(infile,sep="\t")

check=(nrow(mytab)-1)
j=1
while( j<=check ){
  if( (mytab$V1[j] == mytab$V1[j+1]) && (abs(mytab$V2[j+1] - mytab$V3[j]) <= gap) ){
      mytab$V3[j] = mytab$V3[j+1]
      mytab$V4[j] = mytab$V4[j] + mytab$V4[j+1]
      mytab = mytab[-(j+1),]
  }
  else{
    j=j+1}
  check=(nrow(mytab)-1)
}
write.table(mytab,outfile,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

