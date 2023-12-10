#' @title single_gene
#' @description  Drawing the transcript exon and intron structure of a gene in the range.
#' @param annotation A gtf format data,such as Homo_sapiens.GRCh38.109.chr.gtf.
#' @param read A bed format data after alignment with the genome,having the GRanges object and the metadata of name column.
#' @param readmap A tsv format data,mapping the transcript_id to read name;The column name needs to contain the name and transcript id.
#' @param chr The chromosome of interest,such as chr="6".
#' @param range The range of interest,such as range=c(36000000,36700000).
#' @param strand The strand of interest,such as strand="+" or "-" or "both";the default value is "both".
#' @param read_color Displays the color of the read structure.
#' @param isoform_color Displays the color of the transcript isoform structure
#' @param limit To limit the number of read displays. If the read count exceeds the limit value, the limit+log2(read count-limit)  number of reads will be randomly selected.
#' @param text_color The font color of transcript name.
#' @param text_size The font size of transcript name.
#' @param text_alpha The font transparency of transcript name.
#' @param text_fontface The  font shape of transcript name.
#' @param title_color The font color  of tittle name.
#' @param title_size The font size of tittle name.
#' @param title_face The  font shape of tittle name.
#' @param x_text_size The text size of x axis.
#' @param show_transcript_name Whether the transcript name is displayed in the plot.
#' @param num_color The font color of read number.
#' @param num_size The font size of read number.
#' @param num_alpha The font transparency of read number.
#' @param num_fontface The  font shape of read number.
#' @param num_margin The read numbe moves to the left.
#'
#' @return return a list object with ggplot2.
#' @export
#'
#' @examples
#' p_list=single_gene(annotation,MUT1_read,readmap,chr='6',strand="+",range=c(36050000,36700000),show_transcript_name=TRUE)
#' p_list[[1]]
single_gene=function(annotation,read,readmap,chr,range,strand="+",
                     read_color="#ADD8E6",isoform_color="#00DD00",limit=30,
                     text_color="#000000",text_size=3,text_alpha=1,text_fontface="bold",
                     title_color="#000000",title_size=10,title_face="bold",x_text_size=5,show_transcript_name=FALSE,
                     num_color="#000000",num_size=3,num_alpha=1,num_fontface="bold",num_margin=100){

  if (is.null(annotation) || is.null(read) || is.null(readmap) || is.null(chr) || is.null(range)) {
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

  if(!(strand %in% c("+","-","both"))){
    stop("Check the strand parameter !")
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
  grange = GRanges(seqnames = chr,ranges = IRanges(start = range[1],end = range[2]),strand = strand)
  read = single_gene_process_read(read,grange,readmap)
  annotation = subsetByOverlaps(annotation,grange,ignore.strand = FALSE)
  if(length(annotation)==0){
    err=paste("There are't annotation data in ","chr",chr," : ",range[1],"--",range[2],"!")
    stop(err)
  }
  transcript_ids = single_gene_get_transcript_id(annotation,read)
  if(length(transcript_ids)==0){
    err=paste("There are't data in ",strand," : chr",chr," : ",range[1],"--",range[2],"!")
    stop(err)
  }

  #get the common transcript_id's annotation
  anno=annotation[annotation$type=="gene"]
  anno_map=as.data.frame(anno)%>%select(gene_id,gene_name)%>%distinct()
  annotation = single_gene_filter_transcriptID(annotation,transcript_ids)
  #get the common transcript_id's read data
  read= single_gene_filter_transcriptID(read,transcript_ids)
  genemap= single_gene_get_genemap(annotation,transcript_ids,show_transcript_name)
  #genemap$gene_name=anno_map$gene_name[match(anno_map$gene_id,genemap$gene_id)]
  genemap=left_join(genemap,anno_map, by = "gene_id")
  genemap$gene_name[is.na(genemap$gene_name)] <- genemap$gene_id[is.na(genemap$gene_name)]
  gene_ids=unique(genemap$gene_id)
  #add the transcript_name for read data
  read= single_gene_add_name(read,genemap)

  gene_list=split(read,read$gene_id)
  gene_list=lapply(gene_ids,single_gene_process_grange,gene_list,limit)
  gene_data=do.call(rbind,gene_list)
  rownames(gene_data)=NULL
  anno_data=single_gene_process_annotation(annotation,show_transcript_name)
  data=rbind(gene_data,anno_data)
  split_data=split(data,data$gene_id)
  exon_data=setNames(lapply(gene_ids, single_gene_process_data,split_data,genemap),gene_ids)
  p_list=lapply(exon_data, single_gene_draw,read_color,isoform_color,text_color,text_size,text_alpha,text_fontface,title_color,title_size,title_face,x_text_size,num_color,num_size,num_alpha,num_fontface,num_margin)
  return(p_list)
}

#' @keywords internal
#' @noRd
single_gene_process_read=function(read,plus_grange,readmap){
  read = subsetByOverlaps(read,plus_grange,ignore.strand = FALSE)

  colname=colnames(mcols(read))
  if(length(colname)==0 || !"name" %in% colname ){
    stop("There is no name column in the data ！")
  }
  if("blocks" %in% colname){
    read=unlist(rtracklayer::blocks(read))
    read$name = names(read)
    names(read) = NULL
  }
  read$transcript_id = readmap$transcript_id[match(read$name,readmap$name)]
  read <- read[!is.na(mcols(read)$transcript_id)]
  return(read)
}
#' @keywords internal
#' @noRd
single_gene_get_transcript_id=function(grange1,grange2){
  common_ids <- GenomicRanges::intersect(grange1$transcript_id, grange2$transcript_id)
  return(common_ids)
}
#' @keywords internal
#' @noRd
single_gene_filter_transcriptID=function(grange_obj, transcript_ids){
  # Filter the Grange object based on matching transcript_ids
  filtered_grange <- grange_obj[grange_obj$transcript_id %in% transcript_ids]
  return(filtered_grange)
}
#' @keywords internal
#' @noRd
single_gene_get_genemap=function(grange,transcript_ids,show_transcript_name){
  indices <- match(transcript_ids, grange$transcript_id)
  grange <- grange[indices]
  if(show_transcript_name){
    gene_transcript_map=grange %>% as_tibble() %>% dplyr::select(gene_id, transcript_id,transcript_name)
  }
  else{
    gene_transcript_map=grange %>% as_tibble() %>% dplyr::select(gene_id, transcript_id)
    gene_transcript_map$transcript_name=gene_transcript_map$transcript_id
  }
  return(gene_transcript_map)
}
#' @keywords internal
#' @noRd
single_gene_add_name=function(grange_obj,genemap){
  grange_obj$transcript_name=genemap$transcript_name[match(grange_obj$transcript_id,genemap$transcript_id)]
  grange_obj$gene_id=genemap$gene_id[match(grange_obj$transcript_id,genemap$transcript_id)]
  grange_obj$gene_name=genemap$gene_name[match(grange_obj$transcript_id,genemap$transcript_id)]
  return(grange_obj)
}
#' @keywords internal
#' @noRd
single_gene_delete_read=function(grange,limit){
  readname=unique(grange$name)
  num=length(readname)
  if (num > limit) {
    n=ceiling(limit+log2(num-limit))
    keep_indices <- sample(readname, n)
    grange <- grange[grange$name %in% keep_indices, ]
  }
  data=as.data.frame(grange)%>%select(seqnames,start,end,width,strand,name,gene_id,gene_name,transcript_id,transcript_name)
  data$num=num
  data$sort=2
  data$type="exon"
  return(data)
}
#' @keywords internal
#' @noRd
single_gene_process_grange=function(gene_id,gene_list,limit){
  grange=gene_list[[gene_id]]
  split_grange=split(grange,grange$transcript_id)
  split_data=lapply(split_grange,single_gene_delete_read,limit)
  data=do.call(rbind,split_data)
  return(data)
}
#' @keywords internal
#' @noRd
single_gene_process_annotation=function(annotation,show_transcript_name){
  annotation=annotation[annotation$type=="exon"]
  if(show_transcript_name){
    anno_data=as.data.frame(annotation)%>%select(seqnames,start,end,width,strand,gene_id,gene_name,transcript_id,transcript_name)
  }
  else{
    anno_data=as.data.frame(annotation)%>%select(seqnames,start,end,width,strand,gene_id,gene_name,transcript_id)
    anno_data$transcript_name=anno_data$transcript_id
  }
  anno_data$name=paste(anno_data$transcript_id,"annotation")
  anno_data$sort=0
  anno_data$num=0
  anno_data$type="exon"
  new_data=anno_data
  new_data$name=paste(anno_data$transcript_id,"blank")
  new_data$sort=1
  merge_data=rbind(anno_data,new_data)
  return(merge_data)
}
#' @keywords internal
#' @noRd
single_gene_process_data=function(gene_id,data_list,genemap){
  data=data_list[[gene_id]]
  data$gene_name=unique(genemap$gene_name[genemap$gene_id==gene_id])
  transcript_data=split(data,data$transcript_id)
  n=length(transcript_data)
  transcript_data=lapply(transcript_data,function(x){
    x = x %>% arrange(sort,start)
    x$name <- factor(x$name, levels = unique(x$name))
    x$position =  match(x$name, levels(x$name))+3
    max_num=max(x$num)
    x$num=max_num
    x$min_start=min(x$start)
    x$max_end=max(x$end)
    return(x)
  })
  transcript_data[[1]]$text_position = 2
  transcript_data[[1]]$num_position = ceiling((max(transcript_data[[1]]$position)+min(transcript_data[[1]]$position))/2)
  if(n>1){
    for(i in (2:n)){
      position_value=max(transcript_data[[i-1]]$position)
      transcript_data[[i]]$position=transcript_data[[i]]$position + position_value
      transcript_data[[i]]$text_position=min(transcript_data[[i]]$position)-1
      transcript_data[[i]]$num_position = ceiling((max(transcript_data[[i]]$position)+min(transcript_data[[i]]$position))/2)
    }
  }

  merge_data=do.call(rbind,transcript_data)
  rownames(merge_data)=NULL
  return(merge_data)
}
#' @keywords internal
#' @noRd
single_gene_draw=function(exons,read_color,isoform_color,
                          text_color,text_size,text_alpha,text_fontface,
                          title_color,title_size,title_face,x_text_size,
                          num_color,num_size,num_alpha,num_fontface,num_margin,show_transcript_name){

  text_data = exons %>% select(transcript_id,transcript_name,min_start,max_end,text_position,num_position,num,transcript_name) %>% distinct()
  read_exons=exons[exons$sort==2,]
  read_exons$isoform="read"
  anno_exons=exons[exons$sort==0,]
  anno_exons$isoform="isoform"
  exons=rbind(read_exons,anno_exons)
  introns=to_intron(exons,"name")
  exons$size=1
  introns$size=0.4
  title=unique(exons$gene_name)
  xlim=c(min(exons$start),max(exons$end))
  ylim=c(max(exons$position),min(exons$position)-1)
  p = ggplot()+
    geom_text(data=text_data,aes(x = (min_start + max_end) / 2,y=text_position,label = transcript_name),
              color = text_color, size = text_size,alpha=text_alpha,fontface=text_fontface, vjust = 0.5,hjust=0.5,check_overlap = TRUE)+
    geom_text(data=text_data,aes(x = min_start-num_margin ,y=num_position,label = num),
              color = num_color, size = num_size,alpha=num_alpha,fontface=num_fontface, vjust = 0.5,hjust = 1,check_overlap = TRUE)+
    geom_range(data = exons,aes(xstart = start,xend = end,y = position,col =isoform, fill = isoform,height = size))+
    geom_intron(data = introns,aes(xstart = start,xend = end,y = position,col = isoform,size = size),arrow=NULL)+
    theme_classic() +
    scale_y_reverse()+
    scale_color_manual(values =c("isoform"=isoform_color,"read"=read_color))+
    scale_fill_manual(values =c("isoform"=isoform_color,"read"=read_color))+
    coord_cartesian(xlim = xlim,ylim=ylim)+
    scale_size_identity()+
    scale_x_continuous(position = "top")+
    theme(
      axis.text.x = element_text(size = x_text_size,face = "bold"),
      plot.title = element_text(size = title_size,face = title_face,color=title_color,hjust = 0.5),
    )+labs(y = NULL,x=NULL,title=title,color=NULL,fill=NULL,size=NULL)+guides(y="none")
  return(p)
}
