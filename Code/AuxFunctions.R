


readGenericNetwork <- function(nameOfNetwork, method){
  if (method == "Tigramite_"){
    bootscoreName = "d1"
    weightOrDirectionName = "d0"
  }else if (method == "PyCausal_"){
    bootscoreName = "d1"
    weightOrDirectionName = "d0"
  }else{
    bootscoreName = "key_bootscore"
    weightOrDirectionName = "key_weight"
  }
  
  fileR = tryCatch({
    file(nameOfNetwork, "rb")
  }, warning = function(e) {
    NULL
  }, error = function(e) {
    NULL
  })
  network = tryCatch({
    readChar(fileR, file.info(nameOfNetwork)$size) # need to read the terminator explicitly
  }, warning = function(e) {
    ""
  }, error = function(e) {
    ""
  })
  
  if (network == "")
    return (c())

  close(fileR)
  
  
  #READ NETWORK
  network =  gsub("\\(", "", gsub("\\)", "",gsub("\\[", "", gsub("\\]","",gsub("t0", "ti", gsub("tn", "ti+1", tolower(unlist(strsplit(x =network, "\r?\n")))))))))
  network = gsub(paste("\n<data key=\"",bootscoreName,"\">[0-9]*\\.[0-9]*</data>\n",sep=""),"",network)
  newNetwork = network[grep(pattern = paste(".*target=\".*_ti(\\+1)?\">","",sep=""), x = network, ignore.case = T)]
  newNetwork = gsub("_ti\\+1","",gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", gsub("\">", "", gsub("<edge id=\"e?\\d*\" source=\"", "", newNetwork)))))))
  newNetwork = gsub("_ti","",newNetwork)
  
  # newNetwork = newNetwork[setdiff(1:length(newNetwork),grep(".*hg.__.*",newNetwork))]
  
  #See if we need to reverse the names for some pycausal networks
  if (method == "PyCausal_" && length(grep("&",network[grep(paste("<data key=\\\"",bootscoreName,".*",sep=""),network)][1])) == 1){
    bootscoreName = "d0"
    weightOrDirectionName = "d1"
  }
  
  edgeBootscore = network[grep(paste("<data key=\\\"",bootscoreName,".*",sep=""),network)]
  edgeBootscore = as.numeric(gsub("\"", "", gsub(paste("<data key=\\\"",bootscoreName,"\\\">",sep=""), "", gsub("</data>", "", edgeBootscore)))) 
  
  weightOrDirection = network[grep(paste("<data key=\\\"",weightOrDirectionName,".*",sep=""),network)]
  weightOrDirection = gsub("\"", "", gsub(paste(" *<data key=\\\"",weightOrDirectionName,"\\\">",sep=""), "", gsub("</data>", "", weightOrDirection)))
  
  #Extract the information from the network
  pairs = c()
  if (length(newNetwork) == 0){
    return (pairs)
  }
  for (i in 1:length(newNetwork)){
    edge =  c(trimws(strsplit(newNetwork[i],"\" target=\"")[[1]], which ="both"),edgeBootscore[i])
    #We need to take care of the special format of pycausal, witht he direction of the edge as an attribute.
    if (method == "PyCausal_"){
      #If the edge is x--> or x--o, then we add an edge 
      if (substr(weightOrDirection[i],5,5) == "t" || substr(weightOrDirection[i],5,5) == "o"){
        pairs = rbind(pairs,edge)
      } #If the edge is <--x or o--x, then we add an edge but we need to reverse it
      if (substr(weightOrDirection[i],1,1) == "&" || substr(weightOrDirection[i],1,1) == "o"){
        pairs = rbind(pairs,c(edge[2],edge[1],edge[3]))
      }
    }else{
      pairs = rbind(pairs,edge)
    }
  }
  
  return(pairs)
}
