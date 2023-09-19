#' @title samples_comp
#' @description Comparison of transcript isoforms structural differences between multiple samples.
#'
#' @param annotation A gtf format data,such as Homo_sapiens.GRCh38.109.chr.gtf.
#' @param readlist A list object with the read data of multiple samples.
#' @param readmap A tsv format data,mapping the transcript_id to read name;The column name needs to contain the name and transcript id.
#' @param interest_transcript_ids A dataframe with the column name needs to contain interest_transcript_id to choose what you want to show transcript isoform.
#' @param novel_transcript_ids A dataframe with the column name needs to contain novel_transcript_id to mark the color of novel isoform.
#' @param limit To limit the number of read displays. If the read count exceeds the limit value, the limit+log2(read count-limit)  number of reads will be randomly selected.
#' @param read_color Displays the color of the read structure.
#' @param novel_isoform_color Mark the color of the novel transcript isoform.
#' @param known_isoform_color Mark the color of the known transcript isoform.
#' @param text_color The font color of transcript name.
#' @param text_size The font size of transcript name.
#' @param text_alpha The font transparency of transcript name.
#' @param text_fontface The  font shape of transcript name.
#' @param title_color The font color  of tittle name.
#' @param title_size The font size of tittle name.
#' @param title_face The  font shape of tittle name.
#' @param x_text_size The text size of x axis.
#' @param num_color The font color of read number.
#' @param num_size The font size of read number.
#' @param num_alpha The font transparency of read number.
#' @param num_fontface The  font shape of read number.
#' @param num_margin The read numbe moves to the left.
#' @param show_transcript_name Whether the transcript name is displayed in the plot.
#'
#' @return return a ggplot2 object
#' @importFrom ggtranscript to_intron geom_range geom_intron
#' @importFrom cowplot plot_grid
#' @export
#'
#' @examples
#' readlist=list("MUT1"=MUT1_read,"MUT2"=MUT2_read,"WT"=WT_read)
#' interest_transcript_ids=data.frame(interest_transcript_id=c("ENST00000405375","ENST00000615513"))
#' samples_comp(annotation=annotation,readlist=readlist,readmap=readmap,limit=100,interest_transcript_ids=interest_transcript_ids)
samples_comp=function(annotation,readlist,readmap,interest_transcript_ids,novel_transcript_ids=NULL,limit=30,
                      read_color="#ADD8E6",novel_isoform_color="#CC0000",known_isoform_color="#00DD00",
                      text_color="#000000",text_size=3,text_alpha=1,text_fontface="bold",
                      title_color="#000000",title_size=10,title_face="bold",x_text_size=5,
                      num_color="#000000",num_size=3,num_alpha=1,num_fontface="bold",num_margin=100,show_transcript_name=FALSE){

  if (is.null(annotation) || is.null(readlist) || is.null(readmap) || is.null(interest_transcript_ids)) {
    stop("One or more arguments are NULL")
  }
  if (class(annotation) != "GRanges") {
    stop("annotation object is not of GRanges class.")
  }
  if (!(is.data.frame(readmap))) {
    stop("readmap object is not of DataFrame class.")
  }
  readlist <- Filter(Negate(is.null), readlist)
  for(read in readlist){
    if (class(read) != "GRanges") {
      stop("read object is not of GRanges class.")
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
  }
  if(length(readlist)==0){
    stop("input read data is NULL")
  }
  if(is.null(names(readlist))) {
    num_elements <- length(readlist)
    new_names <- paste0("sample", 1:num_elements)
    names(readlist) <- new_names
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
  sample_name=c(names(readlist))
  last_name=sample_name[length(sample_name)]


  transcript_ids=unique(interest_transcript_ids$interest_transcript_id)
  anno_data=compare_samples_process_annotation(annotation,transcript_ids)
  transcript_ids=unique(anno_data$transcript_id)
  readlist=setNames(lapply(readlist,compare_samples_process_read,annotation,anno_data,transcript_ids,readmap=readmap,limit),sample_name)
  map = annotation[annotation$transcript_id %in% transcript_ids & annotation$type == "transcript",]

  map = as.data.frame(map)
  transcript_ids=unique(map$transcript_id)
  if(show_transcript_name){
    map = map %>% select(transcript_id,transcript_name) %>% distinct()
  }
  else{
    map = map %>% select(transcript_id) %>% distinct()
    map$transcript_name=map$transcript_id
  }

  data=compare_samples_process_position(readlist,sample_name,interest_transcript_ids)
  data = data %>% dplyr::filter(name != "blank")
  data$transcript_name=map$transcript_name[match(data$transcript_id,map$transcript_id)]
  xlim=c(min(data$start),max(data$end))
  ylim=c(max(data$position),min(data$position)-1)
  sample_data = split(data,f=data$sample)
  p_list=lapply(sample_name,compare_samples_draw,novel_transcript_ids,sample_data,xlim,ylim,last_name,
                read_color,novel_isoform_color,known_isoform_color,
                text_color,text_size,text_alpha,text_fontface,
                title_color,title_size,title_face,
                x_text_size,
                num_color,num_size,num_alpha,num_fontface,num_margin)

  p=plot_grid(plotlist = p_list,nrow = 1,rel_heights= rep(1,length(p_list)),rel_widths = rep(1,length(p_list)),align = 'h',axis="t")

  return(p)

}
#' @keywords internal
#' @noRd
compare_samples_process_annotation=function(annotation,transcript_ids){
  annotation.exon=annotation[annotation$type == "exon" & annotation$transcript_id %in% transcript_ids]
  tx_ids = unique(annotation.exon$transcript_id)
  diff_elements <- setdiff(transcript_ids,tx_ids)
  if(length(diff_elements)!=0){
    str=paste(diff_elements, collapse = " ")
    err=paste0(str, " ", "does not exist in annotation data")
    stop(err)
  }
  annotation.exon=annotation.exon %>% GenomicRanges::split(.,f=.$transcript_id)
  grange=rtracklayer::asBED(annotation.exon)
  grange=unlist(blocks(grange))
  grange$transcript_id=names(grange)
  names(grange)=NULL

  data=as.data.frame(grange)
  data$name=paste(data$transcript_id,"_annotation")
  data$type="exon"
  data$num=0
  data$sort=0

  blank_data=data
  blank_data$name <- "blank"
  blank_data$sort=1
  data <- rbind(data, blank_data)
  return(data)
}
#' @keywords internal
#' @noRd
compare_samples_process_read=function(grange_obj,annotation,anno_data,transcript_ids,readmap,limit){
  colname=colnames(mcols(grange_obj))
  if(length(colname)==0 || !"name" %in% colname ){
    stop("There is no name column in the data ！")
  }
  if("blocks" %in% colname){
    grange_obj=unlist(rtracklayer::blocks(grange_obj))
    grange_obj$name = names(grange_obj)
    names(grange_obj) = NULL
  }
  grange_obj = compare_samples_add_transcript_id(grange_obj,readmap)
  grange = compare_samples_filter_transcriptID(grange_obj,transcript_ids)
  grange = compare_samples_filter_grange(grange,annotation,transcript_ids)
  grange_list = grange %>% split(.,f=.$transcript_id)
  grange_list = compare_samples_delete_read(grange_list,limit)
  grange=Reduce("c",grange_list)
  data=as.data.frame(grange)
  data = data %>% select(seqnames,start,end,width,strand,name,transcript_id,num)
  rownames(data) == NULL
  data$sort=2
  data$type="exon"
  merge_data=rbind(data,anno_data)
  split_data=merge_data %>% split(.,f=.$transcript_id)
  split_data=lapply(split_data,function(data){
    data$num=max(data$num)
    data = data %>% arrange(sort,start)
    data$name <- factor(data$name, levels = unique(data$name))
    data$position =  match(data$name, levels(data$name))+3
    data$min_start=min(data$start)
    data$max_end=max(data$end)
    return(data)
  })

  merge_data=do.call(rbind,split_data)
  rownames(merge_data)=NULL
  return(merge_data)
}
#' @keywords internal
#' @noRd
compare_samples_process_position=function(readlist,sample_name,interest_transcript_ids){
  for(name in sample_name){
    readlist[[name]]$sample=name
  }
  merge_data=do.call(rbind,readlist)
  transcript_ids=unique(merge_data$transcript_id)

  rownames(merge_data)=NULL
  merge_data <- merge_data %>% mutate(transcript_id = factor(transcript_id, levels = interest_transcript_ids$interest_transcript_id))
  transcript_data = merge_data %>% split(.,f=.$transcript_id)
  transcript_data[[1]]$text_position = 2
  transcript_data[[1]]$num_position = ceiling((max(transcript_data[[1]]$position)+min(transcript_data[[1]]$position))/2)
  n=length(transcript_ids)
  if(n > 1){
    for(i in (2:n)){
      position_value=max(transcript_data[[i-1]]$position)
      transcript_data[[i]]$position=transcript_data[[i]]$position + position_value+1
      transcript_data[[i]]$text_position=min(transcript_data[[i]]$position)-2
      transcript_data[[i]]$num_position = ceiling((max(transcript_data[[i]]$position)+min(transcript_data[[i]]$position))/2)
    }
  }
  merge_data=do.call(rbind,transcript_data)
  rownames(merge_data)=NULL
  return(merge_data)
}

compare_samples_add_transcript_id=function(grange_obj,readmap){
  grange_obj$transcript_id=readmap$transcript_id[match(grange_obj$name,readmap$name)]
  return(grange_obj)
}
#' @keywords internal
#' @noRd
compare_samples_filter_transcriptID=function(grange_obj, transcript_ids){
  # Filter the Grange object based on matching transcript_ids
  filtered_grange <- grange_obj[grange_obj$transcript_id %in% transcript_ids]
  return(filtered_grange)
}
#' @keywords internal
#' @noRd
compare_samples_filter_grange=function(grange_obj,annotation,transcript_ids){
  filtered_range <- annotation[annotation$transcript_id %in% transcript_ids & annotation$type=="transcript"]
  filter_read = subsetByOverlaps(grange_obj,filtered_range,ignore.strand = FALSE)
  return(filter_read)
}
#' @keywords internal
#' @noRd

compare_samples_delete_read=function(grange_list,limit){
  processed_granges <- lapply(grange_list, function(grange) {

    readname=unique(grange$name)
    num=length(readname)
    grange$num=num
    if (num > limit){
      n=ceiling(limit+log2(num-limit))
      keep_indices <- sample(readname, n)
      grange = grange[grange$name %in% keep_indices]
      return(grange)
    }else{
      return(grange)
    }
  })
  return(processed_granges)
}
#' @keywords internal
#' @noRd
compare_samples_draw=function(sample_name,novel_transcript_ids,data,xlim,ylim,last_name,
                              read_color,novel_isoform_color,known_isoform_color,
                              text_color,text_size,text_alpha,text_fontface,
                              title_color,title_size,title_face,x_text_size,
                              num_color,num_size,num_alpha,num_fontface,num_margin){
  exons = data[[sample_name]]
  anno_exons=exons[exons$sort==0,]
  anno_exons$isoform="known isoform"
  anno_exons$isoform[anno_exons$transcript_id %in% novel_transcript_ids$novel_transcript_id]="novel isoform"
  read_exons=exons[exons$sort==2,]
  read_exons$isoform="read"
  exons=rbind(anno_exons,read_exons)
  introns=to_intron(exons,"name")
  text_data = exons %>% select(transcript_name,min_start,max_end,text_position,num_position,num) %>% distinct()
  exons$size=1
  introns$size=0.4
  p = ggplot()+
    geom_text(data=text_data,aes(x = (min_start + max_end) / 2,y=text_position,label = transcript_name),
              color = text_color, size = text_size,alpha=text_alpha,fontface=text_fontface, vjust = 0.5,hjust=0.5,check_overlap = TRUE)+
    geom_text(data=text_data,aes(x = min_start-num_margin ,y=num_position,label = num),
              color = num_color, size = num_size,alpha=num_alpha,fontface=num_fontface, vjust = 0.5,hjust = 1,check_overlap = TRUE)+
    geom_range(data = exons,aes(xstart = start,xend = end,y = position,col =isoform,fill=isoform,height=size))+
    geom_intron(data = introns,aes(xstart = start,xend = end,y = position,col=isoform,size=size),arrow=NULL)+
    theme_classic() +
    scale_y_reverse()+
    scale_color_manual(values =c("novel isoform"=novel_isoform_color,"known isoform"=known_isoform_color,"read"=read_color))+
    scale_fill_manual(values =c("novel isoform"=novel_isoform_color,"known isoform"=known_isoform_color,"read"=read_color))+
    coord_cartesian(xlim = xlim,ylim=ylim)+
    scale_size_identity()+
    scale_x_continuous(position = "top")+labs(y = NULL,x=NULL,title=sample_name,color=NULL,fill=NULL,size=NULL)+
    guides(y="none")+
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = x_text_size,face = "bold"),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = title_size,face = title_face,color=title_color,hjust = 0.5)
    )

  return(p)
}
