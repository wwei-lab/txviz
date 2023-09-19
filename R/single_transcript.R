#' @title single_transcript
#' @description Drawing the exon and intron structure of interest transcript isoform.
#' @param annotation A gtf format data,such as Homo_sapiens.GRCh38.109.chr.gtf.
#' @param read A bed format data after alignment with the genome,having the GRanges object and the metadata of name column.
#' @param readmap A tsv format data,mapping the transcript_id to read name;The column name needs to contain the name and transcript id.
#' @param interest_transcript_id A string data as the interest transcript_id.
#' @param read_color Displays the color of the read structure.
#' @param isoform_color Displays the color of the transcript isoform structure
#' @param limit To limit the number of read displays. If the read count exceeds the limit value,the limit number of reads will be randomly selected.
#' @param title_color The font color  of tittle name.
#' @param title_size The font size of tittle name.
#' @param title_face The  font shape of tittle name.
#' @param num_size The font size of read number.
#' @param x_text_size The text size of x axis.
#' @param show_transcript_name Whether the transcript name is displayed in the plot.
#'
#' @return return a ggplot2 object
#' @export
#'
#' @examples
#' single_transcript(annotation,MUT1_read,readmap,"ENST00000405375",show_transcript_name=TRUE)
single_transcript=function(annotation,read,readmap,interest_transcript_id,
                           read_color="#ADD8E6",isoform_color="#00DD00",
                           limit=30,title_color="#000000",title_size=10,title_face="bold",num_size=10,x_text_size=5,
                           show_transcript_name=FALSE){


  if (is.null(annotation) || is.null(read) || is.null(readmap) || is.null(interest_transcript_id)) {
    stop("One or more  parameter are NULL !")
  }
  if (class(annotation) != "GRanges") {
    stop("annotation object is not of GRanges class.")
  }
  if (class(read) != "GRanges") {
    stop("read object is not of GRanges class.")
  }
  if (!(is.data.frame(readmap))) {
    stop("readmap object is not of DataFrame class.")
  }
  #check the rad column name
  read_colname = colnames(mcols(read))
  if(length(read_colname)==0 || !"name" %in% read_colname ){
    stop("There is no name column in the data ！")
  }
  if("blocks" %in% read_colname){
    read=unlist(rtracklayer::blocks(read))
    read$name = names(read)
    names(read) = NULL
  }
  #check the annotation column name
  anno_colname = colnames(mcols(annotation))
  required_columns <- c("type", "gene_id", "gene_name", "transcript_id")
  missing_columns <- !required_columns %in% anno_colname
  if (any(missing_columns)) {
    stop("The following columns are missing:", paste(required_columns[missing_columns], collapse = ", "))
  }
  if(!("transcript_name" %in% anno_colname) & show_transcript_name){
    stop("There is not the column of transcript_name in the annotation data !")
  }
  transcript_ids=c(interest_transcript_id)

  read$transcript_id = readmap$transcript_id[match(read$name,readmap$name)]
  read = read[!is.na(mcols(read)$transcript_id)]
  read = read[read$transcript_id %in% transcript_ids]
  if(length(read)==0){
    stop(paste("There is no data in ",interest_transcript_id))
  }
  anno=annotation
  annotation = annotation[annotation$type=="exon" & annotation$transcript_id %in% transcript_ids]
  if(length(annotation)==0){
    stop(paste("There is no data in ",interest_transcript_id))
  }
  gene_name=unique(annotation$gene_name)
  gene_id=unique(annotation$gene_id)
  if(show_transcript_name){
    transcript_name=unique(annotation$transcript_name)
  }else{
    transcript_name=interest_transcript_id
  }

  anno = anno[anno$type=="gene" & anno$gene_id==gene_id,]
  read = subsetByOverlaps(read,anno,ignore.strand = FALSE)
  readname=unique(read$name)
  num=length(readname)
  if(num>limit){
    keep_indices <- sample(readname, limit)
    read <- read[read$name %in% keep_indices]
  }
  exons=as.data.frame(read)%>%select(seqnames,start,end,width,strand,name,transcript_id)
  exons = exons %>% arrange(start)
  exons$name <- factor(exons$name, levels = unique(exons$name))
  exons$isoform="read"
  exons$position =  match(exons$name, levels(exons$name))+3
  anno_data = as.data.frame(annotation) %>% select(seqnames,start,end,width,strand,transcript_id)
  anno_data$position=2
  anno_data$name="annotation"
  anno_data$isoform="transcript isoform"
  exons=rbind(exons,anno_data)
  text="tittle"

  introns=to_intron(exons,"name")
  exons$size=1
  introns$size=0.4
  xlim=c(min(exons$start),max(exons$end))
  ylim=c(max(exons$position),min(exons$position)-1)
  p=ggplot()+
    geom_range(data = exons,aes(xstart = start,xend = end,y = position,col=isoform,fill=isoform,height=size))+
    geom_intron(data = introns,aes(xstart = start,xend = end,y = position,col=isoform,size=size),arrow=NULL)+
    theme_classic() +
    scale_y_reverse()+
    scale_color_manual(values =c("transcript isoform"=isoform_color,"read"=read_color))+
    scale_fill_manual(values =c("transcript isoform"=isoform_color,"read"=read_color))+
    coord_cartesian(xlim = xlim,ylim=ylim)+scale_size_identity()+scale_x_continuous(position = "top")+
    theme(
     # axis.text.y=element_blank(),
      #axis.ticks.y=element_blank(),
      axis.title.y = element_text(size=num_size,angle = 0, hjust = 0.5,vjust=0.5),
      axis.text.x = element_text(size = x_text_size,face = "bold"),
     # axis.ticks.x = element_blank(),
      plot.title = element_text(size = title_size,face = title_face,color=title_color,hjust = 0.5),
    )+labs(y = num,x=NULL,title=transcript_name,fill=NULL,color=NULL)+guides(y="none")
  return(p)
}
