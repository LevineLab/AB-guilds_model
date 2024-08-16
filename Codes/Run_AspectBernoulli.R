#Draw in all necessary user inputs using argparser
library(argparser, quietly=TRUE)

p <- arg_parser("Execute a run of the Aspect Bernoulli")
p <- add_argument(p, "--data", help = "Provide the path to a csv-formatted binary annotation data file with inidividual samples in the row dimension and individual annotations in the column dimension.", short = "-d")
p <- add_argument(p, "--num_guilds", help = "Provide the number of guilds, or aspects, you wish to parse your annotation data into.", default = 10, type = "numeric", short = "-c")
p <- add_argument(p, "--seed", help = "Provide a random seed, otherwise the default seed is 123", default = 123, type = "numeric", short = "-s")
p <- add_argument(p, "--num_init", help = "Number of independent initializations of the stochastic optimization to generate (at least 10 is recommended).", default = 10, type = "numeric", short ="-n")
p <- add_argument(p, "--num_iter", help = "Number of iterations of the optimization problem to run per optimization run (NOTE: your runtime will scale on a power law with the size of your data, number of initializations, and number of iterations).", default = 2000, type = "numeric", short = "-i")
p <- add_argument(p, "--cores", help = "If your machine has multiple cores, you can specify the number to use here (this will improve your runtime).", default = 1, type = "numeric")
p <- add_argument(p, "--output", help = "Specify a folder for all of the output files to be stored in. The default will be a subfolder within the current directory.", default = "Output", short = "-o")
p <- add_argument(p, "--guild_size", help = "Specify the number of functions you want to use to describe each guild.", default = 5)

argv <- parse_args(p)

#Attach all necessary packages. This is listed second so that the --help call is cleaner
library(gtools, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(magrittr, quietly=TRUE)
library(gplots, quietly=TRUE)
library(reshape2, quietly=TRUE)
library(tidyverse, quietly=TRUE)

#Need to source the function that calls the AB method itself
source("Codes/aspect_bernoulli.R")
source("Codes/drawmat_precise.R")

K = argv$num_guilds
data<-read_csv(argv$data)
seed = argv$seed

#Now need to run Aspect Bernoulli on this new dataset
#The following line is a check to see if any of the data columns are fully ones, if they are then they can cause the AB to fail and so those functions will be removed (with a text warning flag)
fullones<-which(colSums(data)==dim(data)[1])
if (is_empty(fullones)==FALSE){
  data<-data[,-fullones]
  print(paste("Function(s)",colnames(data)[fullones],"was/were removed because every genome had it (this can cause methodological error)"))
}

#Set a random seed to use, need to change the RNG kind so that it passes through the aspect bernoulli correctly
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
#Run the Aspect Bernoulli
fit<-aspect_bernoulli(data,k=K,num_init=argv$num_init,num_steps=argv$num_iter,em_alg=em_fast,mc.cores=argv$cores)

#Extract from the models which had the highest value in our objective function
ibest <- fit %>% map_dbl(~ .x$obj[length(.x$obj)]) %>% which.max()
bestfit = fit[[ibest]]

## Format them as before
beta = bestfit$beta 
Gamma =  bestfit$gamma
colnames(beta) <- colnames(data)
FF = ncol(beta)
G = nrow(Gamma)
K = nrow(beta)

## Add row and column names
rownames(Gamma) <- 1:G
colnames(Gamma) <- 1:K
rownames(beta) <- 1:K

# beta %>% t() %>% head(15) %>% drawmat_precise(xlab="Guild", ylab="Function")
# Gamma %>% head(15) %>% drawmat_precise(ylab = "Genome", xlab = "Guild")

## ----fig.width=5, fig.height=10----------------------------------------------------
Vmult = Gamma %*% beta
guild_prob_all = array(NA, dim = c(G, K, FF))
for(f in 1:FF){
  for(g in 1:G){
    numers_by_guild = Gamma[g,] * beta[,f]
    denom = Vmult[g,f]
    if(denom==0) denom = denom + 1E-10
    ratio = numers_by_guild/denom
    guild_prob_all[g,,f] = ratio
  }
}
dimnames(guild_prob_all)[[1]] <- 1:G
dimnames(guild_prob_all)[[2]] <- 1:K
dimnames(guild_prob_all)[[3]] <- colnames(beta)

guild_prob = guild_prob_all %>% apply(c(2,3), mean)
# guild_prob %>%
#   t() %>% head(50) %>% drawmat_precise(xlab = "Guild",
#                                        ylab = "Function")
dimnames(guild_prob) = dimnames(beta)


numer = guild_prob
baseline_probs = Gamma %>% colMeans()
denom = matrix(baseline_probs, nrow=nrow(beta), ncol=ncol(beta), byrow=FALSE)
stopifnot(dim(numer) == dim(denom))
score = numer/denom


## ---- eval=TRUE, results = "asis"--------------------------------------------------
function_names = colnames(data)
num = 30
inds <- apply(Gamma, 1, function(myrow){ any(myrow > 0.95)})
mapbacks_guild <- colSums(round(Gamma[inds, ]))
baseline_probs <- colMeans(Gamma)
ScoreMat<-matrix(nrow=num,ncol=K,data=0)
for(guild in 1:K){
  ord = score[guild,] %>% order(decreasing = TRUE) %>% head(num)
  ScoreMat[,guild]<-score[guild,ord]
}

#--------------------------------------------------------------------------------------------------------------#

#Testing a potential final approach, other previous approaches are listed towards the end of the code for different ways of identifying reps based on Gamma and different score modifications. Here we test a 2x(1/K) rep threshold and potential mean-adjusted abundance score modifier. The 2x(1/K) is a post-check AFTER doing the max column to get reps to make sure there aren't any potential close ties
newadjustedscore<-matrix(data=0,nrow=K,ncol=length(function_names))
Gamma_maxes<-max.col(Gamma,ties.method = "random")
avg_fulldata_abundance<-mean(colSums(data))
store_mapbacks<-matrix(data=0,nrow=0,ncol=2)
filename<-paste0(argv$output,"/guild-function-scores_K",K,".csv")

if (file.exists(filename)){
  file.remove(filename)
}
file.create(filename)
totalhits<-0

#Table to store all guild functions and scores in score order
all.guild.info <- data.frame(index = c(1:dim(data)[2]))

#temporary variables to create publication tables
forpub_guildtables<-matrix(data=NA,nrow=K,ncol=6)
for (q in 1:dim(Gamma)[2]){
  reps<-which(Gamma_maxes==q)
  names(reps)<-c()
  #Remove reps below 2x(1/K) threshold
  low_reps<-which(Gamma[reps,q]<(2*(1/K)))
  if(is_empty(low_reps)==FALSE){
    reps<-reps[-low_reps]
  }
  # print(paste("Number of reps discarded is",length(low_reps)))
  # cat("\n")
  
  #Want to construct a mean-adjusted abundance, I think this might reduce the overall impact of abundance which is perhaps too strong right now.
  avg_ingroup_abundance<-mean(colSums(data[reps,]))
  normalize_repratio<-colSums(data[reps,])/avg_ingroup_abundance
  #fix_NaNs<-which(is.nan(normalize_repratio)==TRUE); normalize_repratio[fix_NaNs]<-0
  #newadjustedscore<-score[q,]*rep_ratio
  # newadjustedscore[q,]<-score[q,]*normalize_repratio
  newadjustedscore[q,]<-score[q,]*normalize_repratio
  
  #Temporary code to generate a paper figure where we compare pre and post abundance adjusted scores to show that top functions have many more reps
  adjust_ord<-newadjustedscore[q,]%>%order(decreasing=TRUE) %>% head(75)
  orig_ord<-score[q,]%>%order(decreasing=TRUE) %>% head(75)
  
  adjust_repfreq<-colSums(data[reps,adjust_ord])/length(reps)
  orig_repfreq<-colSums(data[reps,orig_ord])/length(reps)
  freq_df<-as_tibble(cbind(adjust_repfreq,orig_repfreq))
  freq_df<-mutate(freq_df,x=c(1:length(orig_ord)))
  colnames(freq_df)<-c("Adjusted Score","Original Score","Score Rank")
  freq_melt<-freq_df%>%melt(.,id.vars="Score Rank")
  colnames(freq_melt)[c(2,3)]<-c("Score Type","Within Rep Abundance")
  
  #Plotting the scatters
  #print(ggplot(data=freq_melt,aes(x=`Score Rank`,y=`Within Rep Abundance`,color=`Score Type`)) + geom_point())
  
  #Save the names and scores in decreasing order for export
  new_ord<-newadjustedscore[q,] %>% order(decreasing=TRUE) %>% head(20) %>% function_names[.]
  new_ord_numeric<-newadjustedscore[q,] %>% order(decreasing=TRUE) %>% head(20)
  top_scores<-newadjustedscore[q,new_ord_numeric]%>%round(.,4)
  
  full_ord <- newadjustedscore[q,] %>% order(decreasing=TRUE)%>%function_names[.]
  full_ord_numeric <- newadjustedscore[q,] %>% order(decreasing=TRUE)
  full_scores <- newadjustedscore[q,full_ord_numeric]%>%round(.,4)
  
  all.guild.info <- all.guild.info %>% mutate("Guild_{q}" := full_ord, "Guild{q}_Scores" := full_scores)
  
  # plot(newadjustedscore[new_ord_numeric],type="l",xlab="Function",ylab="Score")
  
  #temporary lines for generating publication tables
  add_guild<-new_ord[c(1:5)]
  add_guild<-c(paste("Guild",q),add_guild)
  forpub_guildtables[q,]<-add_guild
  
  #Want to generate heatmaps of top functions and their abundance in genomes
  repmatrix<-as.matrix(data[reps,new_ord_numeric])
  # heatmap(repmatrix,Rowv=NA,Colv=NA,scale="none",col=c("blue","red"),main="Within reps frequency")
  
  #Generate the information on full mapbacks, I guess start with fixed number of functions (5)
  mapback_num<-argv$guild_size
  mapback_func<-new_ord_numeric %>% head(mapback_num)
  #Comparing number of mapbacks between Gamma reps and full data
  full_mapbacks<-which(rowSums(data[reps,mapback_func])==mapback_num)
  alldata_mapbacks<-which(rowSums(data[,mapback_func])==mapback_num)
  
  store_mapbacks<-rbind(store_mapbacks,c(length(full_mapbacks),length(alldata_mapbacks)))
  
  # print("Number of mapbacks within 'reps'")
  # length(full_mapbacks) %>% print()
  # print("Number of mapbacks within full data")
  # length(alldata_mapbacks) %>% print()
  # cat("\n")
  #Also design an approach where we expand in size until no mapbacks?
  
  #Line to show max number of functions in each genome
  #newadjustedscore %>% order(decreasing=TRUE) %>% head(15) %>% data[,.] %>% rowSums() %>% max()
  
  #Writing new score and functions for each guild to a text file in an append approach. Need to add a break or line to denote different iterations
  printdata<-matrix(data=NA,nrow=length(new_ord),ncol=2)
  printdata[,1]<-new_ord; printdata[,2]<-top_scores
  colnames(printdata)<-c("Function","New_Score")
  # write.table(printdata,filename,append=TRUE,sep=",")
}
#Export the full set of guild functions and scores
all.guild.info <- all.guild.info %>% select(-index)
allguildname <- paste0(argv$output,"/all-guild-info_K",K,".csv")
write.table(all.guild.info,allguildname,sep=",")

#Export the number of mapback genomes per guild
store_mapbacks <- as.data.frame(store_mapbacks)
store_mapbacks <- store_mapbacks%>%mutate(Guild = c(1:K))%>%select(Guild,everything())
colnames(store_mapbacks) <- c("Guild","Max Gamma Mapback Genomes","All Data Mapback Genomes")
mapbacktablename <-paste0(argv$output,"/mapback-genomes-per-guild_K",K,".csv")
write.table(store_mapbacks,mapbacktablename,sep=",")

#temporary line for saving publication table files
guildtables_MAG_df<-as.data.frame(forpub_guildtables)%>%t()
guildtablename <- paste0(argv$output,"/top",argv$guild_size,"-guild-functions_K",K,".csv")
write.table(guildtables_MAG_df,guildtablename,sep=",",col.names=FALSE,row.names=FALSE)

