library('stringr')

######plot the sturcture of one gene from genome annotation file (gtf or gff3).##########################


#####-----------------------------------functions------------------------------------------------------

##format each line from file into data frame.
private_line_format<-function(line){
  line_fmt<-strsplit(line,split='\t')
  line_fmt<-as.data.frame(t(as.data.frame(line_fmt)))
  line_fmt$V4<-as.integer(as.character(line_fmt$V4))
  line_fmt$V5<-as.integer(as.character(line_fmt$V5))
  line_fmt$V1<-as.character(line_fmt$V1)
  line_fmt$V2<-as.character(line_fmt$V2)
  line_fmt$V3<-as.character(line_fmt$V3)
  line_fmt$V6<-as.character(line_fmt$V6)
  line_fmt$V7<-as.character(line_fmt$V7)
  line_fmt$V8<-as.character(line_fmt$V8)
  line_fmt$V9<-as.character(line_fmt$V9)
  return(line_fmt)
}
##read the first 2000 lines,judge the filetype and analyse feature.
private_file_type_features_analysis<-function(filename){
  type_list=list()
  #read the first 2000 lines
  gene_info<-as.data.frame(list())
  n=1
  gtf<- file(filename, "r")
  line=readLines(gtf,n=1)
  while( n<500| length(line)==0 ) {
    if (!str_detect(line,'^#')) {
      line_fmt<-private_line_format(line)
      gene_info<-rbind(gene_info,line_fmt)
      n<-n+1
      }
    line=readLines(gtf,n=1)
    }
  close(gtf)
  rownames(gene_info)<-c()
  
  #judge_file_type
  if (str_detect(gene_info$V9[1],'ID=|Parent=')) {type_list[1]<-'gff3'}
  else if (str_detect(gene_info$V9[1],'gene_id')) {type_list[1]<-'gtf'}
  else {return(message("The file type is unclear based on the attributes in collumn 9."))}
  
  #feature analysis
  features<-c()
  n=1
  if ("gene"%in% levels(as.factor(gene_info$V3))) {
    features[n]<-'gene'
    n<-n+1
    }
  if ("mRNA"%in% levels(as.factor(gene_info$V3))) {
    features[n]<-'mRNA'
    n<-n+1
    }  
  if ("transcript"%in% levels(as.factor(gene_info$V3))) {
    features[n]<-'transcript'
    n<-n+1
    }
  if ("exon"%in% levels(as.factor(gene_info$V3))) {
    features[n]<-'exon'
    n<-n+1
    }  
  if ("CDS"%in% levels(as.factor(gene_info$V3))) {
    features[n]<-'CDS'
    n<-n+1
    } 
  type_list[2]<-list(features)
  names(type_list)<-c('file_type','features')
  return(type_list)
}
##judge the last exon of each transcript for plotting direciton, marked in col 6.
private_judge_the_last_exon<-function(transcript_info){
  
  strand_flag<-transcript_info$V7[1]
  
  line_mark=0
  if (strand_flag=='+'){line_mark=0}
  else {line_mark=1e10}
  
  for (i in 1:length(transcript_info$V4)){
    if (transcript_info$V3[i]=="exon"){
      if (strand_flag=="-") {
        if (line_mark>transcript_info$V4[i]){line_mark<-transcript_info$V4[i]}
      }
      else {
        if (line_mark<transcript_info$V5[i]){line_mark<-transcript_info$V5[i]}
      }
    }
  }
  
  for (i in 1:length(transcript_info$V4)) {
    if (transcript_info$V4[i]==line_mark | transcript_info$V5[i]==line_mark){
      transcript_info$V6[i]<-'TTS'
    }
  }
      
  return(transcript_info)
}
##obtain the information of specific gene,the last exon was judged.
gene_info_extract<-function(filename,gene_id){
  
  type_list<-private_file_type_features_analysis(filename)
  ##for_test---------------------------------------
  #gene_id<-"gene-stambpl1"
  #filename<-"./gene_structure/gff3_sample_5000.gff3"
  #type_list<-private_file_type_features_analysis("./gene_structure/gff3_sample_5000.gff3")
  ##for_test_end-----------------------------------
  
  #obtain gene information
  gene_info<-as.data.frame(list())
  if (type_list$file_type=='gtf'){
    gtf<- file(filename, "r")
    line=readLines(gtf,n=1)
    while( length(line) != 0 ) {
      if (str_detect(line,paste('\"',gene_id,'\"',sep=''))
          & (str_detect(line,"\texon\t") 
             | str_detect(line,"\tCDS\t"))){
        line_fmt<-private_line_format(line)
        li_t1<-unlist(strsplit(line_fmt$V9,split=';'))
        li_t1[2]<-str_trim(li_t1[2])
        li_t2<-unlist(strsplit(li_t1[2],split=' '))
        line_fmt$V9<-gsub('\"','',li_t2[2])
        gene_info<-rbind(gene_info,line_fmt)
        }
      line=readLines(gtf,n=1)
    }
    close(gtf) 
    rownames(gene_info)<-c() 
  }
  
  
  if (type_list$file_type=='gff3'){
    gene_flag=F
    gff<- file(filename, "r")
    line=readLines(gff,n=1)
    while( length(line) != 0 ) {
      if (str_detect(line,paste('ID=',gene_id,';',sep=''))){
        gene_flag<-T
        line=readLines(gff,n=1)
        next
        }
      if (gene_flag){
        if (str_detect(line,"\tgene\t") |str_detect(line,"\tpseudogene\t")){
          gene_flag<-F
          break
        }
        if (str_detect(line,"\texon\t") 
               | str_detect(line,"\tCDS\t")){
          line_fmt<-private_line_format(line)
          li_t1<-unlist(strsplit(line_fmt$V9,split=';'))
          li_t2<-unlist(strsplit(li_t1[2],split='='))
          line_fmt$V9<-li_t2[2]
          gene_info<-rbind(gene_info,line_fmt)
        }
      }
      line=readLines(gff,n=1)
    }
    close(gff) 
    rownames(gene_info)<-c()  
  }
  
  #re assemble gene information.
  transcript_id<-levels(as.factor(gene_info$V9))
  gene_info_reassembly<-data.frame()
  gene_start<-1e10
  gene_end<-0
  for (i in 1:length(transcript_id)){
    transcript_info<-data.frame()
    for (j in 1:length(gene_info$V1)){
      if (gene_info$V9[j]==transcript_id[i]){
        transcript_info<-rbind(transcript_info,gene_info[j,])
      }
    }
    transcript_info<-private_judge_the_last_exon(transcript_info)#judge the last exon.
    
    transcript_start<-1e10
    transcript_end<-0
    for (i in 1:length(transcript_info$V4)){
      if (transcript_info$V4[i]<transcript_start) {transcript_start<-transcript_info$V4[i]}
      if (transcript_info$V5[i]>transcript_end) {transcript_end<-transcript_info$V5[i]}
    }
    if (transcript_start<gene_start) {gene_start<-transcript_start}
    if (transcript_end>gene_end) {gene_end<-transcript_end}
    transcript_line<-transcript_info[1,]
    transcript_line$V3[1]<-"transcript"
    transcript_line$V4[1]<-transcript_start
    transcript_line$V5[1]<-transcript_end
    gene_info_reassembly<-rbind(gene_info_reassembly,transcript_line,transcript_info)
  }
  gene_line<-gene_info_reassembly[1,]
  gene_line$V3[1]<-'gene'
  gene_line$V4[1]<-gene_start
  gene_line$V5[1]<-gene_end
  gene_line$V9[1]<-gene_id
  gene_info_reassembly<-rbind(gene_line,gene_info_reassembly)
  rownames(gene_info_reassembly)<-c()  
  return(gene_info_reassembly)
}
##obtain the number and IDs of transcripts from the specific gene.
trancsript_analysis_file<-function(filename,gene_id){ #from an annotation file.
  gene_info<-gene_info_extract(filename,gene_id)
  transcript_num=0
  transcript_name<-c()
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=="transcript"){
      transcript_num<-transcript_num+1
      transcript_name[transcript_num]<-gene_info$V9[i]
    }
  }
  transcripts<-unlist(list(transcript_num,transcript_name))
  return(transcripts)
}
##obtain the number and IDs of transcripts from the specific gene.
trancsript_analysis_df<-function(gene_info){ #from a data frame.
  transcript_num=0
  transcript_name<-c()
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=="transcript"){
      transcript_num<-transcript_num+1
      transcript_name[transcript_num]<-gene_info$V9[i]
    }
  }
  transcripts<-unlist(list(transcript_num,transcript_name))
  return(transcripts)
}
##plot the gene structure.
plot_gene_structure<-function (filename,gene_id){

  ##for_test---------------------------------------
  #gene_id<-"stambpl1"
  #filename<-"./gene_structure/gtf_sample_5000.gtf"
  ##for_test---------------------------------------
  
  gene_info<-gene_info_extract(filename,gene_id)

    #get the length of the gene
  gene_start=gene_info[1,]$V4
  gene_end=gene_info[1,]$V5
  transcript_num=0
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=="transcript"){
      transcript_num<-transcript_num+1
    }
  }
  gene_info<-gene_info[-1,]
  rownames(gene_info)<-c()

  #create an empty plot
  x<-0
  plot(x,type='n',axes=F,
       ylim=c(0,transcript_num*5+3),xlim=c(gene_start-200,gene_end),
       xlab='Gene Position (bp)',ylab='Transcripts',
       main=paste('Gene Structure:',gene_id,sep=' '),
       cex.axis=1,cex.main=2,cex.lab=1
  )
  axis(1,seq(gene_start-500,gene_end+500,500),labels=seq(gene_start-500,gene_end+500,500),cex.axis=1)

  #set colors for exon and CDS
  col_CDS='red'
  col_exon=''
  if ("CDS" %in% gene_info$V3){col_exon='pink'}
  else {col_exon='red'}
  
  #draw the structure of the gene--bone lines.
  cor_Y_start<-3
  for (i in 1:length(gene_info$V1)){
    if (gene_info$V3[i]=="transcript"){
      lines(c(gene_info$V4[i],gene_info$V5[i]),c(cor_Y_start,cor_Y_start),col=rgb(0, 0, 1,0.5),lwd=0.2)
      text(gene_info$V4[i],cor_Y_start+2,labels=gene_info$V9[i],pos=4,cex=0.5)
      cor_Y_start<-cor_Y_start+5
    }
  }
  
  #draw the exon body.
  cor_Y_start<-3
  start_flag<-F
  for (i in 1:length(gene_info$V1)){
    if (start_flag==F){
      start_flag<-T
      next
    }
    if (gene_info$V3[i]=="transcript"){
      cor_Y_start<-cor_Y_start+5
      next
    }
    if (gene_info$V3[i]=="exon") {
      if (gene_info$V6[i]=="TTS"){
        if (gene_info$V7[i]=='-'){
          polygon(c(gene_info$V4[i]-50,gene_info$V4[i],gene_info$V5[i],gene_info$V5[i],gene_info$V4[i]),
                  c(cor_Y_start,cor_Y_start+1.5,cor_Y_start+1.5,cor_Y_start-1.5,cor_Y_start-1.5),
                  col = col_exon, border = NaN)
        }
        else {
          polygon(c(gene_info$V4[i],gene_info$V5[i],gene_info$V5[i]+50,gene_info$V5[i],gene_info$V4[i]),
                  c(cor_Y_start+1.5,cor_Y_start+1.5,cor_Y_start,cor_Y_start-1.5,cor_Y_start-1.5),
                  col = col_exon, border = NaN)
        }

      }
      else {
        polygon(c(gene_info$V4[i],gene_info$V4[i],gene_info$V5[i],gene_info$V5[i]),
                c(cor_Y_start-1.5,cor_Y_start+1.5,cor_Y_start+1.5,cor_Y_start-1.5),
                col = col_exon, border = NaN)
      }
    }
  }
  
  #draw the CDS body.
  if ("CDS" %in% gene_info$V3){
    cor_Y_start<-3
    start_flag<-F
    for (i in 1:length(gene_info$V1)){
      if (start_flag==F){
        start_flag<-T
        next
      }
      if (gene_info$V3[i]=="transcript"){
        cor_Y_start<-cor_Y_start+5
        next
      }
      if (gene_info$V3[i]=="CDS") {
        polygon(c(gene_info$V4[i],gene_info$V4[i],gene_info$V5[i],gene_info$V5[i]),
                c(cor_Y_start-1.5,cor_Y_start+1.5,cor_Y_start+1.5,cor_Y_start-1.5),
                col = col_CDS, border = NaN)
      }
    }
  }

  
}
##plot the gene structure from a 9-collumns table file.
plot_gene_structure_table<-function (gene_info){
  #get the length of the gene
  gene_start<-gene_info[1,]$V4
  gene_end<-gene_info[1,]$V5
  gene_id<-gene_info[1,]$V9
  transcript_num=0
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=="transcript"){
      transcript_num<-transcript_num+1
    }
  }
  gene_info<-gene_info[-1,]
  rownames(gene_info)<-c()

  #create an empty plot
  x<-0
  plot(x,type='n',axes=F,
       ylim=c(0,transcript_num*5+3),xlim=c(gene_start-200,gene_end),
       xlab='Gene Position (bp)',ylab='Transcripts',
       main=paste('Gene Structure:',gene_id,sep=' '),
       cex.axis=1,cex.main=2,cex.lab=1
  )
  axis(1,seq(gene_start-500,gene_end+500,500),labels=seq(gene_start-500,gene_end+500,500),cex.axis=1)
  
  #set colors for exon and CDS
  col_CDS='red'
  col_exon=''
  if ("CDS" %in% gene_info$V3){col_exon='pink'}
  else {col_exon='red'}
  
  #draw the structure of the gene-bone lines.
  cor_Y_start<-3
  for (i in 1:length(gene_info$V1)){
    if (gene_info$V3[i]=="transcript"){
      lines(c(gene_info$V4[i],gene_info$V5[i]),c(cor_Y_start,cor_Y_start),col=rgb(0, 0, 1,0.5),lwd=0.2)
      text(gene_info$V4[i],cor_Y_start+2,labels=gene_info$V9[i],pos=4,cex=0.5)
      cor_Y_start<-cor_Y_start+5
    }
  }
  
  #draw the exon body.
  cor_Y_start<-3
  start_flag<-F
  for (i in 1:length(gene_info$V1)){
    if (start_flag==F){
      start_flag<-T
      next
    }
    if (gene_info$V3[i]=="transcript"){
      cor_Y_start<-cor_Y_start+5
      next
    }
    if (gene_info$V3[i]=="exon") {
      if (gene_info$V6[i]=="TTS"){
        if (gene_info$V7[i]=='-'){
          polygon(c(gene_info$V4[i]-50,gene_info$V4[i],gene_info$V5[i],gene_info$V5[i],gene_info$V4[i]),
                  c(cor_Y_start,cor_Y_start+1.5,cor_Y_start+1.5,cor_Y_start-1.5,cor_Y_start-1.5),
                  col = col_exon, border = NaN)
        }
        else {
          polygon(c(gene_info$V4[i],gene_info$V5[i],gene_info$V5[i]+50,gene_info$V5[i],gene_info$V4[i]),
                  c(cor_Y_start+1.5,cor_Y_start+1.5,cor_Y_start,cor_Y_start-1.5,cor_Y_start-1.5),
                  col = col_exon, border = NaN)
        }
        
      }
      else {
        polygon(c(gene_info$V4[i],gene_info$V4[i],gene_info$V5[i],gene_info$V5[i]),
                c(cor_Y_start-1.5,cor_Y_start+1.5,cor_Y_start+1.5,cor_Y_start-1.5),
                col = col_exon, border = NaN)
      }
    }    
  }
  
  #draw the CDS body.
  if ("CDS" %in% gene_info$V3){
    cor_Y_start<-3
    start_flag<-F
    for (i in 1:length(gene_info$V1)){
      if (start_flag==F){
        start_flag<-T
        next
      }
      if (gene_info$V3[i]=="transcript"){
        cor_Y_start<-cor_Y_start+5
        next
      }
      if (gene_info$V3[i]=="CDS") {
        polygon(c(gene_info$V4[i],gene_info$V4[i],gene_info$V5[i],gene_info$V5[i]),
                c(cor_Y_start-1.5,cor_Y_start+1.5,cor_Y_start+1.5,cor_Y_start-1.5),
                col = col_CDS, border = NaN)
      }
    }
  }
  
  
}
##plot the structure of multi genes.
plot_multi_gene_structure<-function(gene_info,color_flag=TRUE ){
  
  #summary genes
  genes<-data.frame()
  len_shortest_gene<-1e10 #the length of the shortest gene
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=='gene'){
      genes<-rbind(genes,gene_info[i,])
      if (gene_info$V5[i]<len_shortest_gene){len_shortest_gene<-gene_info$V5[i]}
    }
  }
  
  #set color
  col_exon<-c()
  col_cds<-c()
  col_legend<-c()
  if (color_flag) {
    for (i in 1:length(genes$V3)){
      if (i%%6==1){
        col_cds[i]<-rgb(sample(161:255, 1),0,0, alpha=255,maxColorValue=255)
        col_exon[i]<-rgb(sample(100:120, 1),0,0, alpha=255,maxColorValue=255)
      }
      else if (i%%6==2) {
        col_cds[i]<-rgb(0,sample(161:255, 1),0, alpha=255,maxColorValue=255)
        col_exon[i]<-rgb(0,sample(100:120, 1),0, alpha=255,maxColorValue=255)
      }
      else if (i%%6==3) {
        col_cds[i]<-rgb(0,0,sample(161:255, 1), alpha=255,maxColorValue=255)
        col_exon[i]<-rgb(0,0,sample(100:120, 1), alpha=255,maxColorValue=255)
      }
      else if (i%%6==4) {
        col_cds[i]<-rgb(255,sample(161:255, 1),0, alpha=255,maxColorValue=255)
        col_exon[i]<-rgb(255,sample(100:120, 1),0, alpha=255,maxColorValue=255)
      }
      else if (i%%6==5) {
        col_cds[i]<-rgb(0,sample(190:255, 1),255, alpha=255,maxColorValue=255)
        col_exon[i]<-rgb(0,sample(130:150, 1),255, alpha=255,maxColorValue=255)
      }
      else if (i%%6==0) {
        col_cds[i]<-rgb(255,sample(0:60, 1),255, alpha=255,maxColorValue=255)
        col_exon[i]<-rgb(255,sample(100:120, 1),255, alpha=255,maxColorValue=255)
      }
    }
    col_legend<-col_exon
  }
  else {
    for (i in 1:length(genes$V3)){
      col_exon[i]<-"pink"
      col_cds[i]<-"red"
    }
    col_legend<-col_cds
  }


  
  #obtain transcript numbers and gene_id, calculate the length of the genes
  transcript_num=0
  gene_start<-1
  gene_end<-0
  gene_ids<-c()
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=="transcript"){transcript_num<-transcript_num+1}
    if (gene_info$V3[i]=="gene"){
      if (gene_info$V5[i]>gene_end){gene_end<-gene_info$V5[i]}
      gene_ids<-c(gene_ids,gene_info$V9[i])
    }
  }
  
  #create an empty plot
  x<-0
  plot(x,type='n',axes=F,
       ylim=c(0,transcript_num*5+3),xlim=c(gene_start-100,gene_end),
       xlab='Gene Position (bp)',ylab='Transcripts',
       main='Gene Structure',
       cex.axis=1,cex.main=2,cex.lab=1
  )
  axis(1,seq(gene_start-100,gene_end+100,500),labels=seq(gene_start-100,gene_end+100,500),cex.axis=1)
  
  #draw the structure of the gene-bone lines.
  cor_Y_start<-3
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=="transcript"){
      lines(c(gene_info$V4[i],gene_info$V5[i]),c(cor_Y_start,cor_Y_start),col=rgb(0, 0, 1,0.5),lwd=0.2)
      text(gene_info$V4[i],cor_Y_start+2,labels=gene_info$V9[i],pos=4,cex=0.5)
      cor_Y_start<-cor_Y_start+5
    }
  }
  
  #draw exon and cds for each gene
  cor_Y_start_exon<-3
  cor_Y_start_cds<-3
  for (i in 1:length(genes$V9)){
    
    #obtain single gene information
    single_gene_info<-private_single_gene_info(gene_info,genes$V9[i])
    single_gene_info<-single_gene_info[-1,]
    
    #set colors for exon and CDS
    color_cds=col_cds[i]
    color_exon=col_exon[i]
 
    #draw the exon body.
    start_flag<-F
    for (i in 1:length(single_gene_info$V3)){
      if (start_flag==F){
        start_flag<-T
        next
      }
      if (single_gene_info$V3[i]=="transcript"){
        cor_Y_start_exon<-cor_Y_start_exon+5
        next
      }
      if (single_gene_info$V3[i]=="exon") {
        if (single_gene_info$V6[i]=="TTS"){
          if (single_gene_info$V7[i]=='-'){
            polygon(c(single_gene_info$V4[i]-50,single_gene_info$V4[i],single_gene_info$V5[i],single_gene_info$V5[i],single_gene_info$V4[i]),
                    c(cor_Y_start_exon,cor_Y_start_exon+1.5,cor_Y_start_exon+1.5,cor_Y_start_exon-1.5,cor_Y_start_exon-1.5),
                    col = color_exon, border = NaN)
          }
          else {
            polygon(c(single_gene_info$V4[i],single_gene_info$V5[i],single_gene_info$V5[i]+50,single_gene_info$V5[i],single_gene_info$V4[i]),
                    c(cor_Y_start_exon+1.5,cor_Y_start_exon+1.5,cor_Y_start_exon,cor_Y_start_exon-1.5,cor_Y_start_exon-1.5),
                    col = color_exon, border = NaN)
          }
          
        }
        else {
          polygon(c(single_gene_info$V4[i],single_gene_info$V4[i],single_gene_info$V5[i],single_gene_info$V5[i]),
                  c(cor_Y_start_exon-1.5,cor_Y_start_exon+1.5,cor_Y_start_exon+1.5,cor_Y_start_exon-1.5),
                  col = color_exon, border = NaN)
        }
      }    
    }
    cor_Y_start_exon<-cor_Y_start_exon+5
    
    #draw the CDS body.
    if ("CDS" %in% single_gene_info$V3){
      start_flag<-F
      for (i in 1:length(single_gene_info$V3)){
        if (start_flag==F){
          start_flag<-T
          next
        }
        if (single_gene_info$V3[i]=="transcript"){
          cor_Y_start_cds<-cor_Y_start_cds+5
          next
        }
        if (single_gene_info$V3[i]=="CDS") {
          polygon(c(single_gene_info$V4[i],single_gene_info$V4[i],single_gene_info$V5[i],single_gene_info$V5[i]),
                  c(cor_Y_start_cds-1.5,cor_Y_start_cds+1.5,cor_Y_start_cds+1.5,cor_Y_start_cds-1.5),
                  col = color_cds, border = NaN)
        }
      }
      cor_Y_start_cds<-cor_Y_start_cds+5
    }
  }
  #legend
  if (color_flag) {
    legend_y<-3
    for (i in 1:length(gene_info$V3)){
      if (gene_info$V3[i]=='transcript') {legend_y<-legend_y+5}
      if (gene_info$V3[i]=='transcript' & gene_info$V5[i]==len_shortest_gene) {break}
    }
    
    legend(gene_end-3000,legend_y,legend=rev(genes$V9),
           fill=rev(col_legend),col=rev(col_legend),
           text.col=rev(col_legend),
           bty='n',
           border=NaN)
  }

}
##obtain information of a single gene from gene_info data frame.
private_single_gene_info<-function(gene_info,gene_id){
  single_gene_info<-data.frame()
  gene_flag<-F
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V9[i]==gene_id & gene_info$V3[i]=='gene') {
      single_gene_info<-rbind(single_gene_info,gene_info[i,])
      gene_flag=T
      next
    }
    if (gene_info$V9[i]!=gene_id & gene_info$V3[i]=='gene') {
      gene_flag=F
      next
    }
    if (gene_flag) {
      single_gene_info<-rbind(single_gene_info,gene_info[i,])
    }
  }
  
  return(single_gene_info)
}
##format the table from 9-collumn file.
gene_structure_table_read<-function(filename){
  ##for_test---------------------------------------
  #filename<-"./gene_structure/gene_structure_table_stambpl1.txt"
  ##for_test---------------------------------------
  
  #obtain gene information
  gene_info<-as.data.frame(list())
  gene_file<- file(filename, "r")
  line=readLines(gene_file,n=1)
  while( length(line) != 0 ) {
    line_fmt<-private_line_format(line)
    gene_info<-rbind(gene_info,line_fmt)
    line=readLines(gene_file,n=1)
  }
  close(gene_file) 
  rownames(gene_info)<-c() 
  
  #judge the last exon
  line_mark<-0
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=='transcript'){
      if (gene_info$V7[i]=='-'){line_mark<-gene_info$V4[i]}
      else {line_mark<-gene_info$V5[i]}
    }
    if (gene_info$V3[i]=='exon'){
      if (gene_info$V7[i]=='-'){
        if (gene_info$V4[i]==line_mark) {gene_info$V6[i]<-'TTS'}
      }
      else {
        if (gene_info$V5[i]==line_mark) {gene_info$V6[i]<-'TTS'}
      }
    }
  }
  return(gene_info)
}
##data format of genes information: calculate relative position.
gene_data_format<-function(gene_info){
  
  gene_start_flag<-0
  gene_start<-0
  gene_end<-0
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=='gene'){
      if (gene_info$V7[i]=='-'){
        gene_start_flag<-gene_info$V5[i]
        gene_start<-1
        gene_end<-gene_start_flag-gene_info$V4[i]+1
        gene_info$V4[i]<-gene_start
        gene_info$V5[i]<-gene_end
        gene_info$V7[i]<-'+'
      }
      else {
        gene_start_flag<-gene_info$V4[i]
        gene_start<-1
        gene_end<-gene_info$V5[i]-gene_start_flag+1
        gene_info$V4[i]<-gene_start
        gene_info$V5[i]<-gene_end
      }
    }
    else {
      if (gene_info$V7[i]=='-'){
        unit_start<-gene_start_flag-gene_info$V5[i]+1
        unit_end<-gene_start_flag-gene_info$V4[i]+1
        gene_info$V5[i]<-unit_end
        gene_info$V4[i]<-unit_start
        gene_info$V7[i]<-'+'
      }
      else {
        gene_info$V4[i]<-gene_info$V4[i]-gene_start_flag+1
        gene_info$V5[i]<-gene_info$V5[i]-gene_start_flag+1
      }
    }
  }
  return(gene_info)
}
##merge and format data of multi genes. 
multi_genes_data_merge<-function(filename,gene_ids_file){
  
  #for_test-------------------------------------------------
  #filename<-"./gene_structure/GCF_901000725.2_fTakRub1.2_genomic.gff"
  #gene_ids_file<-"./gene_structure/gene_ids.txt"
  #for_test_end------------------------------------------------  
  
  #merge data 
  multi_genes_info<-data.frame()
  gene_ids<- read.table(gene_ids_file, sep='\t', header=F, as.is=T)  #read gene id
  for (i in 1:length(gene_ids$V1)){
    assign(paste('gene_info',gene_ids$V1[i],sep='_'),gene_data_format(gene_info_extract(filename,gene_ids$V1[i])))
    tmp<-get(paste('gene_info',gene_ids$V1[i],sep='_'))
    multi_genes_info<-rbind(multi_genes_info,tmp)
  }
  return(multi_genes_info)
}
##data format of multi genes information: calculate relative position.
multi_genes_data_format<-function(gene_info,sort_flag=TRUE){
  
  #summary genes
  genes<-data.frame()
  for (i in 1:length(gene_info$V3)){
    if (gene_info$V3[i]=='gene'){
      genes<-rbind(genes,gene_info[i,])
    }
  }
  #sort genes based on the length.
  if (sort_flag){genes<-genes[order(genes$V5),]}
  
  
  #extract each gene and format
  multi_genes_info_fmt<-data.frame()
  for (i in 1:length(genes$V3)) {
    multi_genes_info_fmt<-rbind(multi_genes_info_fmt,gene_data_format(private_single_gene_info(gene_info,genes$V9[i])))
  }
 return(multi_genes_info_fmt)
}

#####-----------------------------------functions---end------------------------------------------------
