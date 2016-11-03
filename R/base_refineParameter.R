.refineParameters = function(para){

  if(is.character(para)){
    if(para=="" | tolower(para)=="null" | tolower(para) == "na" ){
      para <- NULL
    }else if(!is.na(as.numeric(para))){
      para <- as.numeric(para)
    }
  }else{
    if(is.na(para))
       para <- NULL
  }
    return(para)

}
