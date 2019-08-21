
options(stringsAsFactors = FALSE)

m_opt = c(commandArgs()[3:length(commandArgs())])

infile=m_opt[1]
outfile=m_opt[2] 
gap=as.numeric(m_opt[3])
blasttype=m_opt[4]
method=m_opt[5]

mytab=read.table(infile,sep="\t")

check=(nrow(mytab)-1)
j=1
while( j<=check ){
#while( mytab$V2[j] != "000094F_pilon"){
  if( mytab$V2[j]==mytab$V2[j+1] && (mytab$V1[j]==mytab$V1[j+1] || blasttype=='any') ){
    if( (method=="overlap") && ( (abs(mytab$V9[j+1]-mytab$V9[j])<gap) && (abs(mytab$V10[j+1]-mytab$V10[j])<gap) ) ){
      mytab$V10[j]=mytab$V10[j+1]
      mytab=mytab[-(j+1),]
    } else if( (method=="flanking") && ( (abs(mytab$V9[j+1]-mytab$V10[j])<gap) ) ){
      mytab$V10[j]=mytab$V10[j+1]
      mytab$V3[j]=((mytab$V3[j]*mytab$V4[j])+(mytab$V3[j+1]*mytab$V4[j+1]))/(mytab$V4[j]+mytab$V4[j+1])
      mytab$V4[j]=mytab$V4[j]+mytab$V4[j+1]
      mytab=mytab[-(j+1),]
    } else if( (method=="anyoverlap") && (mytab$V9[j+1] <= mytab$V10[j]) ){
      mytab$V10[j]=mytab$V10[j+1]
      mytab=mytab[-(j+1),]
    } else{
      j=j+1}
  }
  else{
    j=j+1}
  check=(nrow(mytab)-1)
}
write.table(mytab,outfile,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

