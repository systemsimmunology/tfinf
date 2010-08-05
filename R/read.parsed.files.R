#import parsed scanning data that is located in:
#/local/thorsson/data/TransScan/allMatrix01ParsedFiles
rm(list=ls()); #start fresh

data.dir <- file.path(Sys.getenv("DATA_DIR"),"TransScan/allMouseUpstreamAugust2010/allMouseUpstreamParsedGFF/")
#data.dir <- '/Users/thorsson/data/TransScan/allMouseUpstreamAugust2010/allMouseUpstreamParsedGFF/'
 
file.set <- dir(data.dir);

TRE.OUT <- NULL;
for(ff in 1:length(file.set)){
##for ( ff in 1:3 ){
 sfile <- file.set[ff]; #pick out a particular file
 if ( FALSE ){
   seq.info <- strsplit(sfile,split='-',fixed=TRUE)[[1]][1:3]; #extract the info from the filename
   names(seq.info) <- c('entrez.id','seq.id','tss');
 }
 toks <- strsplit(sfile,split='-',fixed=TRUE)[[1]]
 ## extract the info from the filename
 ## Sep 2009: Some gene names contain a single hyphen and must be treated separately
 if ( length(toks)==3  ){ ## standard case
   seq.info <- toks[1:3];
 } else if (length(toks)==4){ ## gene name has hyphen
   gname <- paste(toks[1:2],collapse="-")
   seq.info <- c(gname,toks[3:4])
 }
 names(seq.info) <- c('name','entrez.id','ensembl.id');
 seq.info["ensembl.id"] <- strsplit(seq.info["ensembl.id"],split=".",fixed=TRUE)[[1]][1]
   
 ##seq.info <- strsplit(sfile,split='-',fixed=TRUE)[[1]][1:3]; #extract the info from the filename
 ##seq.info["ensembl.id"] <- strsplit(seq.info["ensembl.id"],split=".",fixed=TRUE)[[1]][1]
 
 seq.id <- seq.info["ensembl.id"];
 file.p <- paste(data.dir,sfile,sep=''); 
 tre.hits <- NULL;
 if(file.info(file.p)$size>0){ #test whether there is content to read in or not
   if ( FALSE ){ ## Original code for earlier format
     tre.hits <- read.csv(file=file.p,header=FALSE,sep='\t',comment.char="",
                          row.names=NULL,strip.white=TRUE,blank.lines.skip=TRUE,colClasses=c('real','character',rep('real',5)));
     if(unique(tre.hits[,1])!=seq.id){print('TROUBLE MATCHING FILE WITH CONTENT')}; #make sure working as expected
     tre.hits <- tre.hits[,c(2:4,7)];  #only keep what we are interested in
     colnames(tre.hits) <- c('tre','score','thresh','dist');
   }
   ## AP2_Q6_01       0.873391        -1      +
   tre.hits <- read.csv(file=file.p,header=FALSE,sep='\t',comment.char="",
                        row.names=NULL,strip.white=TRUE,blank.lines.skip=TRUE,colClasses=c('character','real','integer','character'));
   tre.hits <- tre.hits[,c(1:4)];  
   colnames(tre.hits) <- c('tre','score','dist','strand');
 }; 
 tre.out <- list(seq.info,tre.hits);
 names(tre.out) <- c('seq.info','tre.hits');
 TRE.OUT <- c(TRE.OUT,list(tre.out));
 names(TRE.OUT)[ff] <- seq.id;
 rm(sfile, seq.info,file.p,tre.hits,seq.id,tre.out);
 write(ff,file='status.txt');
}; #rm(ff);

seq.dir <- file.path(Sys.getenv("TFINF"),"sequence_data")
save(TRE.OUT,file=paste(seq.dir,'Parsed.Scan.Results.allMouseAug2010.dat',sep="/"))

annot.dir <- file.path(Sys.getenv("TFINF"),"annotations")
load(paste(annot.dir,"expressed.scanned.ensembl.RData",sep="/"))

## September 2009. These steps are not needed
## ese.today <- intersect(names(TRE.OUT),expressed.scanned.ensembl) ## Some ensembl IDs have been lost over the two
## years. Sampled some and they are Riken clones etc.
## TRE.OUT.eset <- TRE.OUT[ese.today]
## write.table(ese.today,file="eset",quote=FALSE,row.name=FALSE)

TRE.OUT.eset <- TRE.OUT[expressed.scanned.ensembl]

save(TRE.OUT.eset,file=paste(seq.dir,"Parsed.Scan.eset.allMouseAug2010.RData",sep="/"))
