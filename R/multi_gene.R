#' @title multi_gene
#' @description  Drawing the transcript exon and intron structure of multiple genes at a large scale.
#' @param annotation A gtf format data,such as Homo_sapiens.GRCh38.109.chr.gtf.
#' @param read A bed format data after alignment with the genome,having the GRanges object and the metadata of name column.
#' @param readmap A tsv format data,mapping the transcript_id to read name;The column name needs to contain the name and transcript id.
#' @param novel_transcript_ids A dataframe with the column name needs to contain novel_transcript_id to mark the color of novel isoform.
#' @param chr The chromosome of interest,such as chr="6".
#' @param range The range of interest,such as range=c(36000000,36700000).
#' @param strand The strand of interest,such as strand="+" or "-" or "both";the default value is "both".
#' @param limit  To limit the number of read displays. If the read count exceeds the limit value, the limit+log2(read count-limit)  number of reads will be randomly selected.
#' @param read_low_color The color is the minimum number of reads.
#' @param read_high_color The color is the maximum number of reads.
#' @param annotation_color The annotation exon color in the gtf format data.
#' @param novel_isoform_color Mark the color of the novel transcript isoform.
#' @param known_isoform_color Mark the color of the known transcript isoform.
#' @param text_color The font color of transcript name.
#' @param text_size The font size of transcript name.
#' @param text_alpha The font transparency of transcript name.
#' @param text_fontface The font shape of transcript name.
#' @param title_color The font color  of tittle name.
#' @param title_size The font size of tittle name.
#' @param title_face The font shape of tittle name.
#' @param num_color The font color of read number.
#' @param num_size The font size of read number.
#' @param num_alpha The font transparency of read number.
#' @param num_fontface The font shape of read number.
#' @param num_margin The read numbe moves to the left.
#' @param show_transcript_name Whether the transcript name is displayed in the plot.
#' @param x_text_size The text size of x-axis.
#' @param x_text_just Adjust the x-axis text position.
#' @return return a ggplot2 object
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom ggtranscript to_intron geom_range geom_intron
#' @examples
#' multi_gene(annotation,MUT1_read,readmap,chr='6',strand="both",range=c(36500000,36750000),limit=30,show_transcript_name=TRUE)
#' @export
#'
multi_gene <- function(annotation,read,readmap,novel_transcript_ids=NULL,chr,range,strand="both",limit=30,
                       read_low_color="#ADD8E6",read_high_color="#4682B4",annotation_color="#990099",novel_isoform_color="#CC0000",known_isoform_color="#00DD00",
                       text_color="#000000",text_size=3,text_alpha=1,text_fontface="bold",
                       title_color="#000000",title_size=10,title_face="bold",x_text_size=2,x_text_just=0.25,
                       num_color="#000000",num_size=3,num_alpha=1,num_fontface="bold",num_margin=100,
                       show_transcript_name=FALSE){

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
  if(strand=="+"){
    data_list=multi_gene_process_data(annotation,read,readmap,chr,range,strand,limit,show_transcript_name)
    p=multi_gene_draw(data_list,novel_transcript_ids,chr,strand,range,read_low_color,read_high_color,annotation_color,novel_isoform_color,known_isoform_color,text_color,text_size,text_alpha,text_fontface,title_color,title_size,title_face,x_text_size,x_text_just,num_color,num_size,num_alpha,num_fontface,num_margin)
    return(p)
  }else if(strand=="-"){
    data_list=multi_gene_process_data(annotation,read,readmap,chr,range,strand,limit,show_transcript_name)
    p=multi_gene_draw(data_list,novel_transcript_ids,chr,strand,range,read_low_color,read_high_color,annotation_color,novel_isoform_color,known_isoform_color,text_color,text_size,text_alpha,text_fontface,title_color,title_size,title_face,x_text_size,x_text_just,num_color,num_size,num_alpha,num_fontface,num_margin)
    return(p)
  }else{
    plus_data_list=multi_gene_process_data(annotation,read,readmap,chr,range,strand="+",limit,show_transcript_name)
    minus_data_list=multi_gene_process_data(annotation,read,readmap,chr,range,strand="-",limit,show_transcript_name)
    combined_list <- list()

    for (name in names(plus_data_list)) {
      combined_list[[name]] <- rbind(plus_data_list[[name]], minus_data_list[[name]])
    }
    p=multi_gene_draw(combined_list,novel_transcript_ids,chr,strand,range,read_low_color,read_high_color,annotation_color,novel_isoform_color,known_isoform_color,text_color,text_size,text_alpha,text_fontface,title_color,title_size,title_face,x_text_size,x_text_just,num_color,num_size,num_alpha,num_fontface,num_margin)
    return(p)
  }
}


#' @keywords internal
#' @noRd
multi_gene_get_transcript_id <- function(grange1,grange2){
  common_ids <- intersect(grange1$transcript_id, grange2$transcript_id)
  return(common_ids)
}
#' @keywords internal
#' @noRd
multi_gene_filter_transcriptID <-function(grange_obj, transcript_ids){
  # Filter the Grange object based on matching transcript_ids
  filtered_grange <- grange_obj[grange_obj$transcript_id %in% transcript_ids]
  return(filtered_grange)
}
#' @keywords internal
#' @noRd
multi_gene_get_genemap <-function(grange,transcript_ids,show_transcript_name){
  indices <- match(transcript_ids, grange$transcript_id)
  grange <- grange[indices]
  if(show_transcript_name){
    gene_transcript_map=grange %>% as_tibble() %>% dplyr::select(gene_name,gene_id,transcript_name,transcript_id)
  }else{
    gene_transcript_map=grange %>% as_tibble() %>% dplyr::select(gene_name,gene_id,transcript_id)
    gene_transcript_map$transcript_name=gene_transcript_map$transcript_id
  }
  return(gene_transcript_map)
}
#' @keywords internal
#' @noRd
multi_gene_add_transcript_name <-function(grange_obj,genemap){
  grange_obj$transcript_name=genemap$transcript_name[match(grange_obj$transcript_id,genemap$transcript_id)]
  return(grange_obj)
}
#' @keywords internal
#' @noRd
multi_gene_filter_geneID <- function(grange_obj, gene_ids){
  filtered_grange <- grange_obj[grange_obj$gene_id %in% gene_ids]
  return(filtered_grange)
}
#' @keywords internal
#' @noRd
multi_gene_add_position <- function(datalist){
  n=length(datalist)
  if(n>=2){
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        if(multi_gene_judge_overlap(datalist[[i]],datalist[[j]])){
          strand=unique(datalist[[i]]$strand)
          if(strand=="+"){
            datalist[[j]]$position=datalist[[j]]$position+(max(datalist[[i]]$position)-min(datalist[[i]]$position)+1)+2


          }else{
            datalist[[j]]$position=datalist[[j]]$position+(min(datalist[[i]]$position)-max(datalist[[i]]$position)-1)-2

          }
        }
      }
    }
  }
  return(datalist)
}
#' @keywords internal
#' @noRd
multi_gene_judge_overlap <- function(data1, data2) {
  value1 <- c(min(data1$start), max(data1$end))
  value2 <- c(min(data2$start), max(data2$end))
  position1 <- c(min(abs(data1$position)), max(abs(data1$position)))
  position2 <- c(min(abs(data2$position)), max(abs(data2$position)))

  if (value1[1] <= value2[2] && value1[2] >= value2[1] &&
      position1[1] <= position2[2] && position1[2] >= position2[1]) {
    return(TRUE)
  }

  return(FALSE)
}
#' @keywords internal
#' @noRd
multi_gene_delete_read <- function(grange,num,limit){
  readname=unique(grange$name)
  if (num > limit) {
    n=ceiling(limit+log2(num-limit))
    keep_indices <- sample(readname, n)
    grange <- grange[grange$name %in% keep_indices, ]
  }
  return(grange)
}
#' @keywords internal
#' @noRd
multi_gene_process_data <- function(annotation,read,readmap,chr,range,strand,limit,show_transcript_name){
  grange =GenomicRanges::GRanges(seqnames = chr,ranges = IRanges::IRanges(start = range[1],end = range[2]),strand = strand)
  #filter the read data based on minus_grange
  read = multi_gene_process_read(read,grange,readmap)
  #filter the annotation based on minus_grange
  annotation = subsetByOverlaps(annotation,grange,ignore.strand = FALSE)
  if(length(annotation)==0){
    err=paste("There are't annotation data in ","chr",chr," : ",range[1],"--",range[2],"!")
    stop(err)
  }
  gene_annotation=annotation
  #get the common transcript_id
  transcript_ids = multi_gene_get_transcript_id(annotation,read)
  if(length(transcript_ids)==0){
    err=paste("There are't data in ",strand," : chr",chr," : ",range[1],"--",range[2],"!")
    stop(err)
  }
  #get the common transcript_id's annotation
  annotation = multi_gene_filter_transcriptID(annotation,transcript_ids)
  #get the common transcript_id's read data
  read=multi_gene_filter_transcriptID(read,transcript_ids)
  #get the map of gene_id、gene_name、transcript_id、transcript_name
  genemap=multi_gene_get_genemap(annotation,transcript_ids,show_transcript_name)
  #add the transcript_name for read data
  read=multi_gene_add_transcript_name(read,genemap)
  #get the list of grange object based on transcript_id
  grange_list=read %>% split(.,f = .$transcript_id)
  gene_ids=unique(genemap$gene_id)
  gene_annotation=multi_gene_filter_geneID(gene_annotation,gene_ids)
  gene_data = multi_gene_process_gene(gene_annotation,strand)
  value1=max(gene_data$position)
  if(strand=="-"){
    value1=min(gene_data$position)
  }
  gene_exon = multi_gene_process_gene_exon(gene_annotation,strand,value1,show_transcript_name)
  data_list=lapply(grange_list,multi_gene_process_grange,annotation,strand,value1,limit)
  get_sort_value <- function(df) {
    return(max(df$num))
  }
  data_list <- data_list[order(sapply(data_list, get_sort_value))]
  data_list=multi_gene_add_position(data_list)
  read_exon = do.call(rbind,data_list)
  text_data = read_exon %>% select(transcript_id,transcript_name,min_start,max_end,num,position) %>% distinct() %>% split(.,f=.$transcript_id)
  text_data=lapply(text_data,function(text){
    if(strand=="+"){
      text$name_position=min(text$position)-1
      text$num_position=ceiling((max(text$position)+min(text$position))/2)+1
      text$position=0
    }else{
      text$name_position=max(text$position)+1
      text$num_position=floor((max(text$position)+min(text$position))/2)-1
      text$position=0
    }
    return(text)
  })
  text_data=do.call(rbind,text_data)
  text_data = text_data %>% distinct()
  rownames(gene_data)=NULL
  rownames(gene_exon)=NULL
  rownames(read_exon)=NULL
  rownames(text_data)=NULL
  return(list(gene=gene_data,gene_exon=gene_exon,read_exon=read_exon,text=text_data))
}
#' @keywords internal
#' @noRd
multi_gene_process_read <- function(read,grange,readmap){
  read = subsetByOverlaps(read,grange,ignore.strand = FALSE)
  colname=colnames(mcols(read))
  if(length(colname)==0 || !"name" %in% colname ){
    stop("There is no name column in the data ！")
  }
  if("blocks" %in% colname){
    read=base::unlist(rtracklayer::blocks(read))
    read$name = names(read)
    names(read) = NULL
  }
  read$transcript_id = readmap$transcript_id[match(read$name,readmap$name)]
  read <- read[!is.na(mcols(read)$transcript_id)]
  return(read)
}
#' @keywords internal
#' @noRd
multi_gene_process_gene <- function(gene_annotation,strand){
  gene = gene_annotation[gene_annotation$type=="gene"]
  data=as.data.frame(gene) %>% select(seqnames,start,end,width,strand,type,gene_id,gene_name)

  if(strand=="-"){
    data$position = -3

  }else{
    data$position=3
  }
  datalist=split(data,data$gene_id)
  datalist=multi_gene_add_position(datalist)

  process_data=do.call(rbind,datalist)
  if(strand=="-"){
    process_data$text_position = process_data$position+1

  }else{
    process_data$text_position = process_data$position-1
  }
  return(process_data)
}
#' @keywords internal
#' @noRd
multi_gene_process_gene_exon <- function(gene_annotation,strand,position,show_transcript_name){

  exon = gene_annotation[gene_annotation$type=="exon"]
  if(show_transcript_name){
    exon_data=as.data.frame(exon) %>% select(seqnames,start,end,width,strand,type,gene_id,gene_name,transcript_id,transcript_name)
  }
  else{
    exon_data=as.data.frame(exon) %>% select(seqnames,start,end,width,strand,type,gene_id,gene_name,transcript_id)
    exon_data$transcript_name=exon_data$transcript_id
  }
  if(strand=="+"){
    exon_data$position = position+2
  }
  else{
    exon_data$position = position-2
  }
  return(exon_data)
}
#' @keywords internal
#' @noRd
multi_gene_process_grange <- function(grange,annotation,strand,position,limit){
  num_read = length(unique(grange$name))
  txid=unique(grange$transcript_id)
  txname=unique(grange$transcript_name)
  grange$type = "exon"
  grange$num = num_read
  grange$sort = 2
  grange=multi_gene_delete_read(grange,num_read,limit)
  exon=as.data.frame(grange)
  exon=exon %>% select(seqnames,start,end,width,strand,transcript_id,transcript_name,name,type,num,sort)
  anno_data=multi_gene_process_annotation(annotation,txid,txname)
  exon=rbind(anno_data,exon)
  exon$min_start=min(exon$start)
  exon$max_end=max(exon$end)
  exon = exon %>% arrange(sort,start)
  exon$name <- factor(exon$name, levels = unique(exon$name))
  if(strand=="+"){
    exon$position =  match(exon$name, levels(exon$name))+position+4
  }else{
    exon$position =  -1*match(exon$name, levels(exon$name))+position-4
  }
  exon$num=num_read
  return(exon)
}
#' @keywords internal
#' @noRd
multi_gene_process_annotation <- function(annotation,transcript_id,transcript_name){
  annotation.exon=annotation[annotation$type == "exon" & annotation$transcript_id == transcript_id]
  annotation.exon=annotation.exon %>% GenomicRanges::split(.,f=.$transcript_id)
  grange=rtracklayer::asBED(annotation.exon)
  grange=unlist(blocks(grange))
  names(grange)=NULL
  exon=as.data.frame(grange)
  exon$transcript_id=transcript_id
  exon$transcript_name=transcript_name
  exon$name=paste("annotation_",transcript_id)
  exon$type="exon"
  exon$num=0
  exon$sort=0

  random_index <- sample(nrow(exon), 1)
  new_row <- exon[random_index, ]
  new_row$name <- paste("blank_",transcript_id)
  new_row$sort=1

  exon <- rbind(exon, new_row)
  return(exon)
}
#' @keywords internal
#' @noRd
multi_gene_draw <- function(data_list,novel,chr,strand,xlim,read_low_color,read_high_color,annotation_color,novel_isoform_color,known_isoform_color,
                            text_color,text_size,text_alpha,text_fontface,
                            title_color,title_size,title_face,x_text_size,x_text_just,
                            num_color,num_size,num_alpha,num_fontface,num_margin){

  gene=data_list$gene
  gene$gene_name <- ifelse(is.na(gene$gene_name), "none", gene$gene_name)
  gene_exon=data_list$gene_exon
  gene_exon$isoform="annotation"
  gene_exon$size=1
  exons=data_list$read_exon
  text_data=data_list$text
  text_data$transcript_name <- ifelse(is.na(text_data$transcript_name), "none", text_data$transcript_name)
  anno_exons=exons[exons$sort==0,]
  anno_exons$isoform="known isoform"
  anno_exons$isoform[anno_exons$transcript_id %in% novel$novel_transcript_id]="novel isoform"
  anno_introns=to_intron(anno_exons,"name")
  anno_exons$size=1
  anno_introns$size=0.4
  read_exons=exons[exons$sort==2,]
  read_exons$count=log(read_exons$num)
  read_exons$count <- as.numeric(read_exons$count)
  read_introns=to_intron(read_exons,"name")
  read_exons$size=1
  read_introns$size=0.4
  max_abs_value <- max(abs(read_exons$position))
  #ylim=c(-max_abs_value,max_abs_value)
  range=c(format(xlim[1], scientific = FALSE),format(xlim[2], scientific = FALSE))
  x_lim <- data.frame(
    x=seq(range[1], range[2], length.out = 5),
    label = sapply(seq(range[1], range[2], length.out = 5),function(x) format(x,scientific = FALSE))
  )
  p = ggplot()+
    geom_hline(yintercept = 0,linetype = "solid", color = "black")+
    geom_text(data =x_lim , aes(x = x, y = x_text_just, label = label,size=x_text_size), vjust = 0)+
    #draw gene line
    geom_segment(data=gene,aes(x=start,xend=end,y=position,yend=position),color="#000000",size=1)+
    #draw gene name
    geom_text(data=gene,aes(x = (start + end) / 2,y=text_position,label = gene_name),
              color = text_color, size = text_size,alpha=text_alpha, fontface=text_fontface,vjust = 0.5,hjust=0.5)+
    #draw gene range
    geom_range(data=gene_exon,aes(xstart = start,xend = end,y = position,col=isoform,fill=isoform,height=size))+
    #draw annotation range

    geom_range(data = anno_exons,aes(xstart = start,xend = end,y = position,color=isoform,fill=isoform,height=size))+
    geom_intron(data = anno_introns,aes(xstart = start,xend = end,y = position,color=isoform,size=0.3),arrow=NULL)+

    scale_color_manual(values =c("novel isoform"=novel_isoform_color,"known isoform"=known_isoform_color,"annotation"=annotation_color))+
    scale_fill_manual(values =c("novel isoform"=novel_isoform_color,"known isoform"=known_isoform_color,"annotation"=annotation_color))+
    new_scale_color()+
    new_scale_fill()+

    geom_range(data = read_exons,aes(xstart = start,xend = end,y = position,col=count,fill=count,height=size))+
    geom_intron(data = read_introns,aes(xstart = start,xend = end,y = position,col=count,size=size),arrow=NULL)+
    scale_color_gradient(low = read_low_color, high = read_high_color)+
    scale_fill_gradient(low = read_low_color, high = read_high_color)+
    #draw transcript_name
    geom_text(data=text_data,aes(x = (min_start + max_end) / 2,y=name_position,label = transcript_name),
              color = text_color, size = text_size,alpha=text_alpha,fontface=text_fontface, vjust = 0.5,hjust=0.5,check_overlap = TRUE)+
    #draw read num
    geom_text(data=text_data,aes(x = min_start-num_margin ,y=num_position,label = num),
              color = num_color, size = num_size,alpha=num_alpha,fontface=num_fontface, vjust = 0.5,hjust = 1,check_overlap = TRUE)+
    theme_classic() +
    coord_cartesian(xlim = xlim) +
    scale_size_identity()+
    theme(
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x = element_text(size = x_text_size,face = "bold"),
      plot.title = element_text(size = title_size, face = title_face,color=title_color,hjust = 0.5),
      legend.title.align = 0.5,
    )+labs(x=NULL,y=NULL,title=paste("chr",chr," : ",range[1],"--",range[2]),
           color = "log2(read count)",fill="log2(read count)")+
    guides(x="none",color = guide_legend(title.position = "top"),fill = guide_legend(title.position = "top"),y="none")
  return(p)
}
