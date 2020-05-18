####Libraries####
library("h2o", lib.loc="~/R/win-library/3.5")
library("elmNN", lib.loc="~/R/win-library/3.5")
library("pROC", lib.loc="~/R/win-library/3.5")
library("caret", lib.loc="~/R/win-library/3.5")
library("mclust", lib.loc="~/R/win-library/3.5")
library("kohonen", lib.loc="~/R/win-library/3.5")
library("VennDiagram", lib.loc="~/R/win-library/3.5")
library("rlist", lib.loc="~/R/win-library/3.5")
####WD####
Output_Dir=choose.dir(default="",caption = "Set output folder")
setwd(Output_Dir)
####Choose Files####

Seq_File=choose.files(default = "", caption = "Select peptide sequence table")
MM_File=choose.files(default = "", caption = "Select MM files")
Seq_Filename=basename(Seq_File)
Seq_ID=substr(Seq_Filename,1,(nchar(Seq_Filename)-4))

h2o.init(nthreads=-1, max_mem_size="2G")
h2o.removeAll()
  ####MM import#####

altschul<-read.csv(MM_File[1], header=TRUE,row.names = 1)
blosum62<-read.csv(MM_File[2], header=TRUE,row.names = 1)
dayhoof<-read.csv(MM_File[3], header=TRUE,row.names = 1)
gonnet<-read.csv(MM_File[4], header=TRUE,row.names = 1)
grantham<-read.csv(MM_File[5], header=TRUE,row.names = 1)
henikoff<-read.csv(MM_File[6], header=TRUE,row.names = 1)
johnson<-read.csv(MM_File[7], header=TRUE,row.names = 1)
jones<-read.csv(MM_File[8], header=TRUE,row.names = 1)
levin<-read.csv(MM_File[9], header=TRUE,row.names = 1)
mclachlan<-read.csv(MM_File[10], header=TRUE,row.names = 1)
MMID<- c('altschul','blosum62','dayhoof','gonnet','grantham','henikoff','johnson','jones','levin','mclachlan')
MMIDShort<-c('A','B','D','G','R','H','J','N','L','M','All')
####Sequence import####

Raw_Data<-read.table(Seq_File, header=TRUE)
Raw_Rows<-nrow(Raw_Data)
Pep_Length<-nchar(as.character(Raw_Data[1,1]))
MasterKey=sample(1:5,size = Raw_Rows,replace = TRUE)

####Data Partitioning####



####Signal Version data sample####  
#
#TrueRaw_Data<-read.table(Seq_File, header=TRUE)
#TrueRaw_Rows<-nrow(TrueRaw_Data)
#TruePep_Length<-nchar(as.character(TrueRaw_Data[1,1]))
#
#TrueRawCases=which(TrueRaw_Data[,2]==1)
#TrueRawConts=which(TrueRaw_Data[,2]==0)
#TrueRawMasterkey=sample(1:10,size = length(TrueRawCases),replace = TRUE)
#TrueCases=TrueRaw_Data[TrueRawCases,]
#TrueConts=TrueRaw_Data[TrueRawConts,]
#UsedCases=TrueCases[which(TrueRawMasterkey==1),]
#
#
#Raw_Data<-rbind(UsedCases,TrueConts)
#rownames(Raw_Data)<-NULL
#Raw_Rows<-nrow(Raw_Data)
#Pep_Length<-nchar(as.character(Raw_Data[1,1]))
#MasterKey=sample(1:5,size = Raw_Rows,replace = TRUE)
#
####K fold####
h2ommauc<-list()
h2odeepauc<-list()
h2oorthauc<-list()
TotalSimtraining<-data.frame()
TotalSimtesting<-data.frame()
TestingSequences<-vector(mode = "character",length = 0)
kfold=vector()
TestSequencesKey=data.frame()

for(k in 1:5){
  k=2
##arbitrary for testing (will need to fix for k-fold cross-validaition)
CurrentTesting<- Raw_Data[MasterKey==k,]
CurrentTraining <- Raw_Data[MasterKey!=k,]

####Similarity Calculation#### 
#
#reoragnise for deep learningso all biobases i same row (all fetures in same itteration)

#blank data frame

simnames=vector(length=Pep_Length)
  Orthnames=vector(length=Pep_Length)

 
  Similarity_Training=data.frame("Data_Set","Input_Sequence","BioBase","MM","Class",stringsAsFactors=FALSE)
  Similarity_Testing=data.frame("Data_Set","Input_Sequence","BioBase","MM","Class",stringsAsFactors=FALSE)
  Orth_Training=data.frame("Data_Set","Input_Sequence","Class",stringsAsFactors=FALSE)
  Orth_Testing=data.frame("Data_Set","Input_Sequence","Class",stringsAsFactors=FALSE)
  

  
  
  for (a in 1:Pep_Length) { 
    simnames[a]<-paste0("Similarity ",a) 
    Orthnames[a]<-a
    Similarity_Training<-cbind(Similarity_Training,paste0("Similarity ",a))
    Similarity_Testing<-cbind(Similarity_Testing,paste0("Similarity ",a))
    Orth_Training<-cbind(Orth_Training,a)
    Orth_Testing<-cbind(Orth_Testing,a)
  }
  
 Similarity_Training<-Similarity_Training[FALSE,]
  Similarity_Testing<-Similarity_Testing[FALSE,]
  Orth_Training<-Orth_Training[FALSE,]
  Orth_Testing<-Orth_Testing[FALSE,]
  
  AllMMSimTraining<-data.frame()
  AllMMSimTesting<-data.frame()
  names(Similarity_Training)<-c("Data_Set","Input_Sequence","BioBase","MM","Class",simnames)
  names(Similarity_Testing)<-c("Data_Set","Input_Sequence","BioBase","MM","Class",simnames)
  names(Orth_Training)<-c("Data_Set","Input_Sequence","Class",Orthnames)
  names(Orth_Testing)<-c("Data_Set","Input_Sequence","Class",Orthnames)
  
#####similarity calcuations####
  #for each MM
  for (m in 1:10) {
  CurrentMM<-get(as.character(MMID[m]))
   

  Biobase_Key<-sample(which(CurrentTraining[,2]==1),10,replace = FALSE)
  BioBases<- data.frame("BioBases"=CurrentTraining[Biobase_Key,1])
  Training_Set<- data.frame("BioBases"=CurrentTraining[-Biobase_Key,1],"label"=CurrentTraining[-Biobase_Key,2],stringsAsFactors = FALSE)
 
   AllMMminSimTraining<-data.frame()
   AllMMminSimTesting<-data.frame()
  
   AllMMpartTraining<-data.frame("Data"=rep(Seq_ID,nrow(Training_Set)),"MM"=MMID[m],"Class"=as.integer(Training_Set[,2]))
   AllMMpartTesting<-data.frame("Data"=rep(Seq_ID,nrow(CurrentTesting)),"MM"=MMID[m],"Class"=as.integer(CurrentTesting[,2]))
  
    #for each biobase
    for (b in 1:length(BioBases[,1])) {
    CurrentBB=as.character(BioBases[b,1])
    AllMMminSimTraining<-data.frame()
    AllMMminSimTesting<-data.frame()
      for (s  in 1:nrow(Training_Set)){
      Current_Seq= as.character(Training_Set[s,1])
      TempSimVector=c(Seq_ID,Current_Seq,CurrentBB,MMID[m],as.integer(Training_Set[s,2]))
      CurrentSimSeq<-vector()
        for (i in 1:Pep_Length) {
        CurrentSimNum=CurrentMM[substr(Current_Seq,i,i),substr(CurrentBB,i,i)]
        CurrentSimSeq<-c(CurrentSimSeq,CurrentSimNum)
        TempSimVector<-c(TempSimVector,as.integer(CurrentSimNum))
        }
      Similarity_Training<-rbind(Similarity_Training,setNames(as.list(TempSimVector), names(Similarity_Training)),stringsAsFactors=FALSE)
    AllMMminSimTraining<-rbind(AllMMminSimTraining,CurrentSimSeq)
     names(AllMMminSimTraining)<-as.character(seq(1,Pep_Length,1))
      }
    AllMMpartTraining<-cbind(data.frame(AllMMpartTraining),as.data.frame(AllMMminSimTraining))
   
    
   
  for (t  in 1:nrow(CurrentTesting)){
    Current_Test_Seq= as.character(CurrentTesting[t,1])
    TestTempSimVector=c(Seq_ID,Current_Test_Seq,CurrentBB,MMID[m],as.integer(CurrentTesting[t,2]))
    CurrentTestSimSeq<-vector()
    for (j in 1:Pep_Length) {
      CurrentTestSimNum=CurrentMM[substr(Current_Test_Seq,j,j),substr(CurrentBB,j,j)]
      CurrentTestSimSeq<-c(CurrentTestSimSeq,CurrentTestSimNum)
      TestTempSimVector<-c(TestTempSimVector,CurrentTestSimNum)
    }
    Similarity_Testing<-rbind(Similarity_Testing,setNames(as.list(TestTempSimVector), names(Similarity_Testing)),stringsAsFactors=FALSE)
    AllMMminSimTesting<-rbind(AllMMminSimTesting,CurrentTestSimSeq)
    names(AllMMminSimTesting)<-as.character(seq(1,Pep_Length,1))
    
  } 
    AllMMpartTesting<-cbind(data.frame(AllMMpartTesting),as.data.frame(AllMMminSimTesting))
    }
   AllMMSimTraining<-rbind(AllMMSimTraining,AllMMpartTraining)
   AllMMSimTesting<-rbind(AllMMSimTesting,AllMMpartTesting)
   }


####Orth Data####



#use 20 long binary for input
for (s  in 1:nrow(CurrentTraining)){ 
  Current_Seq= as.character(CurrentTraining[s,1])
  TempOrthVector=c(Seq_ID,Current_Seq,as.integer(CurrentTraining[s,2]))
  for (i in 1:Pep_Length) {
    CurrentOrthAA=gsub("A",00000000000000000001,gsub("C",00000000000000000010,gsub("D",00000000000000000100,gsub("E",00000000000000001000,
                         gsub("F",00000000000000010000,gsub("G",00000000000000100000,gsub("H",00000000000001000000,gsub("I",00000000000010000000,
                                  gsub("K",00000000000100000000,gsub("L",00000000001000000000,gsub("M",00000000010000000000,gsub("N",00000000100000000000,
                                      gsub("P",00000001000000000000,gsub("Q",00000010000000000000,gsub("R",00000100000000000000,gsub("S",00001000000000000000,
                                                                                                                                                                                                                                                                                                                                                                                                              gsub("T",00010000000000000000,gsub("V",00100000000000000000,gsub("W",01000000000000000000,gsub("Y",10000000000000000000,substr(Current_Seq,i,i)))))))))))))))))))))
    TempOrthVector<-c(TempOrthVector,CurrentOrthAA)
  }
  Orth_Training<-rbind(Orth_Training,setNames(as.list(TempOrthVector), names(Orth_Training)),stringsAsFactors=FALSE)
}

for (s  in 1:nrow(CurrentTesting)){ 
  Current_Seq= as.character(CurrentTesting[s,1])
  TempOrthVector=c(Seq_ID,Current_Seq,as.integer(CurrentTesting[s,2]))
  for (i in 1:Pep_Length) {
    CurrentOrthAA=gsub("A",00000000000000000001,gsub("C",00000000000000000010,gsub("D",00000000000000000100,gsub("E",00000000000000001000,
                                                                                                                 gsub("F",00000000000000010000,gsub("G",00000000000000100000,gsub("H",00000000000001000000,gsub("I",00000000000010000000,
                                                                                                                                                                                                                gsub("K",00000000000100000000,gsub("L",00000000001000000000,gsub("M",00000000010000000000,gsub("N",00000000100000000000,
                                                                                                                                                                                                                                                                                                               gsub("P",00000001000000000000,gsub("Q",00000010000000000000,gsub("R",00000100000000000000,gsub("S",00001000000000000000,
                                                                                                                                                                                                                                                                                                                                                                                                              gsub("T",00010000000000000000,gsub("V",00100000000000000000,gsub("W",01000000000000000000,gsub("Y",10000000000000000000,substr(Current_Seq,i,i)))))))))))))))))))))
    
    TempOrthVector<-c(TempOrthVector,CurrentOrthAA)
  }
  Orth_Testing<-rbind(Orth_Testing,setNames(as.list(TempOrthVector), names(Orth_Testing)),stringsAsFactors=FALSE)
}
  #####allsims##
  #TotalSimtraining=rbind(TotalSimtraining,AllMMSimTraining)
  #TotalSimtesting=rbind(TotalSimtesting,AllMMSimTesting)
  
####

  #write.csv(Similarity_Training,paste0(Output_Dir,"\\Training Sim ",k,".csv"))
  #write.csv(Similarity_Testing,paste0(Output_Dir,"\\Testing Sim ",k,".csv"))

####Begin simulations####
  ####elm####

  for (w in 1:10) {
    w=1
    ModelInput<-Similarity_Training[which(Similarity_Training[,4]==MMID[w]),c(5:(5+Pep_Length))]
    ModelInput[]<-lapply(ModelInput, as.integer)
    ModelInput[,1]<-as.logical(ModelInput[,1])
 
    
    ModelTest<-Similarity_Testing[which(Similarity_Testing[,4]==MMID[w]),c(6:(5+Pep_Length))]
    ModelTest[]<-lapply(ModelTest, as.integer)
    
    
   model <- elmtrain(x=ModelInput[,2:(1+Pep_Length)], y=ModelInput[1], nhid=10, actfun="sig")
   modelpredictionElm <- predict(model,newdata=ModelTest[,1:Pep_Length])
    ELMPredictionsDF<- cbind(Similarity_Testing[which(Similarity_Testing[,4]==MMID[w]),],as.data.frame(modelpredictionElm))
    #ELMROC<-auc(roc(as.vector(ELMOrthPredictionsDF$Class),ELMOrthPredictionsDF$Class0))
    write.csv(ELMPredictionsDF,paste0(Output_Dir,"\\b-elm Predictions ",k," - ",MMID[w],".csv"))
  }
 
    
  ####elm all####
  
    
    
    ModelInput<-AllMMSimTraining[,c(3:(3+(10*Pep_Length)))]
    ModelInput[]<-lapply(ModelInput, as.integer)
    ModelInput[,1]<-as.logical(ModelInput[,1])
    
    
    ModelTest<-AllMMSimTesting[,c(3:(3+(10*Pep_Length)))]
    ModelTest[]<-lapply(ModelTest, as.integer)
    
    
    model <- elmtrain(x=ModelInput[,2:(1+(10*Pep_Length))], y=ModelInput[1], nhid=10, actfun="sig")
    modelpredictionElmDBB <- predict(model,newdata=ModelTest[,2:(1+(10*Pep_Length))])
    ELMDBBPredictionsDF<- cbind(AllMMSimTesting,as.data.frame(modelpredictionElmDBB))
    write.csv(ELMDBBPredictionsDF,paste0(Output_Dir,"\\b-elm Predictions ",k," ELM all MM.csv"))
  
    ####elm Orth tester####
    
    
    
    #ModelInput<-Orth_Training[,c(3:13)]
    #ModelInput[,c(2:11)]<-lapply(X=ModelInput[,c(2:11)], FUN= factor,levels=c("1","10","100","1000","10000","1e+05","1e+06","1e+07","1e+08","1e+09","1e+10","1e+11","1e+12","1e+13","1e+14","1e+15","1e+16","1e+17","1e+18","1e+19"))
    #ModelInput[,c(2:11)]<-lapply(X=ModelInput[,c(2:11)],as.numeric)
  # 
  #   
  #  dmy <- dummyVars(" ~ .", data = ModelInput)
  #  trsf <- data.frame(predict(dmy, newdata = ModelInput))
  #  
  #  
  #  ModelTest<-Orth_Testing[,c(3:13)]
  #  
  # 
  #  ModelTest[,c(2:11)]<-lapply(X=ModelTest[,c(2:11)], FUN= factor,levels=c("1","10","100","1000","10000","1e+05","1e+06","1e+07","1e+08","1e+09","1e+10","1e+11","1e+12","1e+13","1e+14","1e+15","1e+16","1e+17","1e+18","1e+19"))
  #   ModelTest[,c(2:11)]<-lapply(X=ModelTest[,c(2:11)],as.numeric)
  #  dmytest <- dummyVars(" ~ .", data = ModelTest)
  #  testrsf <- data.frame(predict(dmytest, newdata = ModelTest))
  #  
  #  
  #  model <- elmtrain(x=trsf[2:ncol(trsf)], y=trsf[1], nhid=10, actfun="sig")
  #  modelpredictionElmOrth <- predict(model,newdata=testrsf[2:ncol(testrsf)])
  #  ELMOrthPredictionsDF<- cbind(Orth_Testing,as.data.frame(modelpredictionElmOrth))
  #  
  #  #auc(roc(as.vector(ELMOrthPredictionsDF$Class),ELMOrthPredictionsDF$Class0))
  #  #plot.roc(roc(as.vector(ELMOrthPredictionsDF$Class),ELMOrthPredictionsDF$Class0))
  #  
  #  
  #  write.csv(ELMOrthPredictionsDF,paste0(Output_Dir,"\\C-Orth-elm Predictions",k,".csv"))
  #  
    
    ####elm Orth####
    
 
    
    ModelInput<-Orth_Training[,c(3:(3+Pep_Length))]
    ModelInput[,c(2:(1+Pep_Length))]<-lapply(X=ModelInput[,c(2:(1+Pep_Length))], FUN= factor,levels=c("1","10","100","1000","10000","1e+05","1e+06","1e+07","1e+08","1e+09","1e+10","1e+11","1e+12","1e+13","1e+14","1e+15","1e+16","1e+17","1e+18","1e+19"))
    ModelInput[,1]<-as.numeric(ModelInput[,1])
    
    
    dmy <- dummyVars(" ~ .", data = ModelInput)
    trsf <- data.frame(predict(dmy, newdata = ModelInput))
    
    ModelTest<-Orth_Testing[,c(3:(3+Pep_Length))]
    ModelTest[,c(2:(1+Pep_Length))]<-lapply(X=ModelTest[,c(2:(1+Pep_Length))], FUN= factor,levels=c("1","10","100","1000","10000","1e+05","1e+06","1e+07","1e+08","1e+09","1e+10","1e+11","1e+12","1e+13","1e+14","1e+15","1e+16","1e+17","1e+18","1e+19"))
    ModelTest[,1]<-as.numeric(ModelTest[,1])
    
    
    dmytest <- dummyVars(" ~ .", data = ModelTest)
    testrsf <- data.frame(predict(dmytest, newdata = ModelTest))
   
   model <- elmtrain(x=trsf[2:ncol(trsf)], y=trsf[1], nhid=10, actfun="sig")
    modelpredictionElmOrth <- predict(model,newdata=testrsf[2:ncol(testrsf)])
    ELMOrthPredictionsDF<- cbind(Orth_Testing,as.data.frame(modelpredictionElmOrth))
    
    
    write.csv(ELMOrthPredictionsDF,paste0(Output_Dir,"\\C-Orth-elm Predictions ",k,".csv"))
    
    
####h20####
  h2o.removeAll()
  for (n in 1:10) {
  ModelInput<-Similarity_Training[which(Similarity_Training[,4]==MMID[n]),c(5:(5+Pep_Length))]
  ModelInput[]<-lapply(ModelInput, as.integer)
  for (q in 1:nrow(ModelInput)) { 
    if(ModelInput[q,1]==1){ModelInput[q,1]<- "Cleaved"} else {ModelInput[q,1]<-"Uncleaved"}}
  ModelInput[,1]<-as.factor(ModelInput[,1])
  
  
  
  ModelTest<-Similarity_Testing[which(Similarity_Testing[,4]==MMID[n]),c(5:(5+Pep_Length))]
   ModelTest[]<-lapply(ModelTest, as.integer)
   for (q in 1:nrow(ModelTest)) { 
     if(ModelTest[q,1]==1){ModelTest[q,1]<- "Cleaved"} else {ModelTest[q,1]<-"Uncleaved"}}
   ModelTest[,1]<-as.factor(ModelTest[,1])
   
   
  
   ModelInputh2o <- as.h2o(ModelInput)
   ModelTesth2o <- as.h2o(ModelTest)
   model_dl <- h2o.deeplearning(x = 2:(1+Pep_Length), 
                                y = 1, 
                                training_frame = ModelInputh2o,
                                seed=123456,
                                hidden = 10,
                                validation_frame = ModelTesth2o)
   #h2ommAUC<-cbind(h2ommAUC,h2o.auc(model_dl, valid = TRUE))
   
    modelpredictionh2o <- h2o.predict(model_dl, ModelTesth2o)
    h2oPredictionsDF<- cbind(Similarity_Testing[which(Similarity_Testing[,4]==MMID[n]),],as.data.frame(modelpredictionh2o))
    #h2oROC<-roc(as.vector(h2oPredictionsDF$Class),h2oPredictionsDF$Cleaved)
   write.csv(h2oPredictionsDF,paste0(Output_Dir,"\\B-h2o Predictions ",k," - ",MMID[n],".csv"))
   #png(paste0("ROC for",k," - ",MMID[n],".csv"))
   plot(h2o.performance(model_dl,newdata = ModelTesth2o))
   #legend("bottomright",legend = paste0("AUC value = ",h2o.auc(model_dl)))
   #dev.off()
  }


##h20 Deep BB####
h2o.removeAll()
  
  ModelInput<-AllMMSimTraining[,c(3:(3+(10*Pep_Length)))]
  ModelInput[]<-lapply(ModelInput, as.integer)
  for (q in 1:nrow(ModelInput)) { 
    if(ModelInput[q,1]==1){ModelInput[q,1]<- "Cleaved"} else {ModelInput[q,1]<-"Uncleaved"}}
  ModelInput[,1]<-as.factor(ModelInput[,1])
  
  ModelTest<-AllMMSimTesting[,c(3:(3+(10*Pep_Length)))]
  ModelTest[]<-lapply(ModelTest, as.integer)
  for (q in 1:nrow(ModelTest)) { 
    if(ModelTest[q,1]==1){ModelTest[q,1]<- "Cleaved"} else {ModelTest[q,1]<-"Uncleaved"}}
  ModelTest[,1]<-as.factor(ModelTest[,1])
  
  ModelInputh2oDBB <- as.h2o(ModelInput)
  ModelTesth2oDBB <- as.h2o(ModelTest)
  model_dl <- h2o.deeplearning(x = 2:(1+(10*Pep_Length)), 
                               y = 1, 
                               training_frame = ModelInputh2oDBB,
                               seed=123456,
                               hidden = 10,
                               validation_frame = ModelTesth2oDBB)

  
  #h2odeepauc<-cbind(h2odeepauc,h2o.auc(model_dl, valid = TRUE))
  modelpredictionh2oDBB <- h2o.predict(model_dl, ModelTesth2oDBB)
  h2oDBBPredictionsDF<- cbind(AllMMSimTesting,as.data.frame(modelpredictionh2oDBB))
  #h2oDBBROC<-auc(roc(as.vector(h2oOrthPredictionsDF$Class),h2oOrthPredictionsDF$predict))
  write.csv(h2oDBBPredictionsDF,paste0(Output_Dir,"\\B-h2o Predictions ",k," deepBB.csv"))

  plot(h2o.performance(model_dl,newdata = ModelTesth2oDBB))
####h2o orth####
  
  ModelInput<-Orth_Training[,c(3:(3+Pep_Length))]
  
  ModelInput[,c(2:(1+Pep_Length))]<-lapply(X=ModelInput[,c(2:(1+Pep_Length))], FUN= factor,levels=c("1","10","100","1000","10000","1e+05","1e+06","1e+07","1e+08","1e+09","1e+10","1e+11","1e+12","1e+13","1e+14","1e+15","1e+16","1e+17","1e+18","1e+19"))
  ModelInput[,1]<-as.numeric(ModelInput[,1])
  
  dmy <- dummyVars(" ~ .", data = ModelInput)
  trsf <- data.frame(predict(dmy, newdata = ModelInput))
  
  
  ModelTest<-Orth_Testing[,c(3:(3+Pep_Length))]
  ModelTest[,c(2:(1+Pep_Length))]<-lapply(X=ModelTest[,c(2:(1+Pep_Length))], FUN= factor,levels=c("1","10","100","1000","10000","1e+05","1e+06","1e+07","1e+08","1e+09","1e+10","1e+11","1e+12","1e+13","1e+14","1e+15","1e+16","1e+17","1e+18","1e+19"))
  ModelTest[,1]<-as.numeric(ModelTest[,1])
  
  dmytest <- dummyVars(" ~ .", data = ModelTest)
  testrsf <- data.frame(predict(dmytest, newdata = ModelTest))
  
  
  ModelInputOrth <- as.h2o(trsf)
  ModelTestOrth <- as.h2o(testrsf)
  model_dl <- h2o.deeplearning(x = 2:ncol(trsf), 
                               y = 1, 
                               training_frame = ModelInputOrth,
                               seed=123456,
                               hidden = 10,
                               validation_frame = ModelTestOrth)
  
  

 # h2oorthauc<-cbind(h2oorthauc,h2o.auc(model_dl, valid = TRUE))
  
  modelpredictionh2oOrth <- h2o.predict(model_dl, ModelTestOrth[2:ncol(testrsf)])
  h2oOrthPredictionsDF<- cbind(Orth_Testing,as.data.frame(modelpredictionh2oOrth))
  
  #auc(roc(as.vector(h2oOrthPredictionsDF$Class),h2oOrthPredictionsDF$predict))
  
  write.csv(h2oOrthPredictionsDF,paste0(Output_Dir,"\\c- Orth h2o Predictions ",k,".csv"))
  
  TestingSequences=append(TestingSequences,as.character(CurrentTesting[,1]))
  kfold=append(kfold,rep(k,length(CurrentTesting[,1])))
} 
TestSequencesKey=data.frame('TestSeq'=TestingSequences,'Fold'=kfold)



####Comparison data ####


library("h2o", lib.loc="~/R/win-library/3.5")
library("elmNN", lib.loc="~/R/win-library/3.5")
library("pROC", lib.loc="~/R/win-library/3.5")
library("caret", lib.loc="~/R/win-library/3.5")
library("mclust", lib.loc="~/R/win-library/3.5")
library("kohonen", lib.loc="~/R/win-library/3.5")
library("VennDiagram", lib.loc="~/R/win-library/3.5")
library("rlist", lib.loc="~/R/win-library/3.5")

Output_Dir=choose.dir(default="",caption = "Set output folder")
setwd(Output_Dir)


Seq_File=choose.files(default = "", caption = "Select peptide sequence table")
MM_File=choose.files(default = "", caption = "Select MM files")
Seq_Filename=basename(Seq_File)
Seq_ID=substr(Seq_Filename,1,(nchar(Seq_Filename)-4))
##
MMID<- c('altschul','blosum62','dayhoof','gonnet','grantham','henikoff','johnson','jones','levin','mclachlan')
MMIDShort<-c('A','B','D','G','R','H','J','N','L','M','All')
MMNum=c(0,1,2,3,4,5,6,7,8,9)

Raw_Data<-read.table(Seq_File, header=TRUE)
Raw_Rows<-nrow(Raw_Data)
Pep_Length<-nchar(as.character(Raw_Data[1,1]))

  FinalData<-list()
    outputFileNames=dir(Output_Dir)
    TestSeqVector=vector()
    TestSequencesKey=data.frame()
     TestSequencesKey <-data.frame(stringsAsFactors = FALSE)
    for(o in 1:5){ path<-paste(Output_Dir,"\\",outputFileNames[(o*10)],sep="")
     OutputDataframe=as.list(read.csv(path,header=TRUE,row.names = 1))
     TestSeqVector=append(TestSeqVector,as.character(OutputDataframe[[2]][1:(0.1*length(OutputDataframe[[2]]))]))
    }
    TestSequencesKey=data.frame("TestSeq"=TestSeqVector)
     
     
     for(i in 1:55){
      path<-paste(Output_Dir,"\\",outputFileNames[i],sep="")
      OutputDataframe=as.list(read.csv(path,header=TRUE,row.names = 1))
      if(i==11|i==22|i==33|i==44|i==55){
      FinalData<-as.list(cbind(FinalData,OutputDataframe[3],OutputDataframe[4+(Pep_Length*10)]))  }else{
          FinalData<-as.list(cbind(FinalData,OutputDataframe[5],OutputDataframe[6+Pep_Length]))
      } }
    
    
 
    for(j in 56:110){
      path<-paste(Output_Dir,"\\",outputFileNames[j],sep="")
      OutputDataframe=as.list(read.csv(path,header=TRUE,row.names = 1))
      if(j==66|j==77|j==88|j==99|j==110){
        FinalData<-as.list(cbind(FinalData,OutputDataframe[3],OutputDataframe[5+(Pep_Length*10)]))}else{
        FinalData<-as.list(cbind(FinalData,OutputDataframe[5],OutputDataframe[7+Pep_Length]))
        } }
      
   
    
    for(l in 111:115){
      path<-paste(Output_Dir,"\\",outputFileNames[l],sep="")
      OutputDataframe=as.list(read.csv(path,header=TRUE,row.names = 1))
      FinalData<-as.list(cbind(FinalData,OutputDataframe[3],OutputDataframe[4+Pep_Length]))
       }
    for(l in 116:120){
      path<-paste(Output_Dir,"\\",outputFileNames[l],sep="")
      OutputDataframe=as.list(read.csv(path,header=TRUE,row.names = 1))
      FinalData<-as.list(cbind(FinalData,OutputDataframe[3],OutputDataframe[4+Pep_Length]))
       }
    
    
    
    
#namevector=c(rep(c('elmaltschul1','elmaltschul2','elmblosum621','elmblosum622','elmdayhoof1','elmdayhoof2','elmgonnet1','elmgonnet2','elmgrantham1','elmgrantham2','elmhenikoff1','elmhenikoff2','elmjohnson1','elmjohnson2','elmjones1','elmjones2','elmlevin1','elmlevin2','elmmclachlan1','elmmclachlan2','elmall1','elmall2'),5),rep(c('h2oaltschul1','h2oaltschul2','h2oblosum621','h2oblosum622','h2odayhoof1','h2odayhoof2','h2ogonnet1','h2ogonnet2','h2ograntham1','h2ograntham2','h2ohenikoff1','h2ohenikoff2','h2ojohnson1','h2ojohnson2','h2ojones1','h2ojones2','h2olevin1','h2olevin2','h2omclachlan1','h2omclachlan2','h2oall1','h2oall2'),5))
#names(FinalData)<-namevector
  
    
#####  Prediction Classification  ####

    ##
    BelmAllMatrixInfo=list()
    B_elmdata=data.frame()
    for(x in 0:9){
      positionelmclass=c(1,23,45,67,89)+(2*x)
      positionelmpred=c(2,24,46,68,90)+(2*x)
      elmAVG<-list()
      elmClass<-list()
      
      for(z in 1:5){
        elmtempAVG=vector()
        elmTClass=vector()
        elmPredsraw=as.vector(FinalData[[positionelmpred[z]]])
        elmClassraw=as.vector(FinalData[[positionelmclass[z]]])
        for (W in 1:(length(elmPredsraw)*0.1)) {
          
          tempavg=mean(elmPredsraw[W+(MMNum*0.1*length(elmPredsraw))]) 
          elmtempAVG=c(elmtempAVG,tempavg)
          elmTClass=c(elmTClass,elmClassraw[W])
        }
        elmAVG[[z]]<-elmtempAVG
        elmClass[[z]]<-elmTClass
      }
      
      
      elmValues<-data.frame()
      MatrixInfo<-NULL
      
      for(z in 1:5){
        
        elmROC<-roc(response=as.factor(elmClass[[z]]),predictor =elmAVG[[z]],direction = "<")
        elmValues<-rbind(elmValues,data.frame(MMIDShort[x+1],auc(elmROC)))
        
        Type=NULL
        CurrentInfo<-list(Posclass=which(elmClass[[z]]==1),NegClass =which(elmClass[[z]]==0),
                          PosPred =which(elmAVG[[z]]>=0.5),NegPred=which(elmAVG[[z]]<0.5))
        
        for(y in 1:length(elmClass[[z]])) {
          #4=true pos, , 3=false pos, 2=true neg, 1= false neg
          Type=c(Type,(if((y %in% CurrentInfo$PosPred)==TRUE){if((y %in% CurrentInfo$Posclass)==TRUE){4}else{3}}else{if((y %in% CurrentInfo$NegClass)==TRUE){2}else{1}}))
          
        }
        MatrixInfo[z]=as.data.frame(Type)
      }
      for(p in 1:length(MatrixInfo)){
        BelmAllMatrixInfo[(5*x)+p]=(as.data.frame(MatrixInfo[[p]]))
      }
      B_elmdata<-rbind(B_elmdata,data.frame(paste('elm',MMIDShort[x+1]),mean(elmValues[,2]),sd(elmValues[,2])))
    } 
    
    names(elmValues)<-c('ID','AUC')
    names(B_elmdata)<-c('ID','Mean AUC', 'SD AUC')
    

    
    
    
    
    
##h2o    
    Bh2oAllMatrixInfo=list()
    B_h2odata=data.frame()
    for(x in 0:9){
      positionh2oclass=c(111,133,155,177,199)+(2*x)
      positionh2opred=c(112,134,156,178,200)+(2*x)
      h2oAVG<-list()
      h2oClass<-list()
      
      for(z in 1:5){
        h2otempAVG=vector()
        h2oTClass=vector()
        h2oPredsraw=as.vector(FinalData[[positionh2opred[z]]])
        h2oClassraw=as.vector(FinalData[[positionh2oclass[z]]])
        for (W in 1:(length(h2oPredsraw)*0.1)) {
          
          tempavg=mean(h2oPredsraw[W+(MMNum*0.1*length(h2oPredsraw))]) 
          h2otempAVG=c(h2otempAVG,tempavg)
          h2oTClass=c(h2oTClass,h2oClassraw[W])
        }
        h2oAVG[[z]]<-h2otempAVG
        h2oClass[[z]]<-h2oTClass
      }
      
      
      h2oValues<-data.frame()
      MatrixInfo<-NULL
      
      for(z in 1:5){
        
        h2oROC<-roc(response=as.factor(h2oClass[[z]]),predictor =h2oAVG[[z]],direction = "<")
        h2oValues<-rbind(h2oValues,data.frame(MMIDShort[x+1],auc(h2oROC)))
        
        Type=NULL
        CurrentInfo<-list(Posclass=which(h2oClass[[z]]==1),NegClass =which(h2oClass[[z]]==0),
                          PosPred =which(h2oAVG[[z]]>=0.5),NegPred=which(h2oAVG[[z]]<0.5))
        
        for(y in 1:length(h2oClass[[z]])) {
          #4=true pos, , 3=false pos, 2=true neg, 1= false neg
          Type=c(Type,(if((y %in% CurrentInfo$PosPred)==TRUE){if((y %in% CurrentInfo$Posclass)==TRUE){4}else{3}}else{if((y %in% CurrentInfo$NegClass)==TRUE){2}else{1}}))
          
        }
        MatrixInfo[z]=as.data.frame(Type)
      }
      for(p in 1:length(MatrixInfo)){
        Bh2oAllMatrixInfo[(5*x)+p]=(as.data.frame(MatrixInfo[[p]]))
        }
      B_h2odata<-rbind(B_h2odata,data.frame(paste('h2o',MMIDShort[x+1]),mean(h2oValues[,2]),sd(h2oValues[,2])))
    } 
    
    names(h2oValues)<-c('ID','AUC')
    names(B_h2odata)<-c('ID','Mean AUC', 'SD AUC')
  
    

    
    
    
  #orth 
    
    Bh2oOrthMatrixInfo=list()
    Orth_h2odata=data.frame()
      positionOh2oclass=c(221,223,225,227,229)
      positionOh2opred=c(222,224,226,228,230)
      Orthh2oValues<-data.frame()
      MatrixInfo<-NULL
      for(z in 1:5){
        Orthh2oROC<-roc(as.factor(FinalData[[positionOh2oclass[z]]]),predictor =FinalData[[positionOh2opred[z]]],direction = "<")
        Orthh2oValues<-rbind(Orthh2oValues,data.frame(MMIDShort[x+1],auc(Orthh2oROC)))
        
        Type=NULL
        CurrentInfo<-list(Posclass=which(FinalData[[positionOh2oclass[z]]]==1),NegClass =which(FinalData[[positionOh2oclass[z]]]==0),
                          PosPred =which(FinalData[[positionOh2opred[z]]]>=0.5),NegPred=which(FinalData[[positionOh2opred[z]]]<0.5))
        
        for(y in 1:length(FinalData[[positionOh2oclass[z]]])) {
          #4=true pos, , 3=false pos, 2=true neg, 1= false neg
          Type=c(Type,(if((y %in% CurrentInfo$PosPred)==TRUE){if((y %in% CurrentInfo$Posclass)==TRUE){4}else{3}}else{if((y %in% CurrentInfo$NegClass)==TRUE){2}else{1}}))
          
        }
        MatrixInfo[z]=as.data.frame(Type)
      }
      for(p in 1:length(MatrixInfo)){
        Bh2oOrthMatrixInfo[p]=(as.data.frame(MatrixInfo[[p]]))
      }
      Orth_h2odata<-rbind(Orth_h2odata,data.frame('h2oorth',mean(Orthh2oValues[,2]),sd(Orthh2oValues[,2])))
    
    names(Orthh2oValues)<-c('ID','AUC')
    names(Orth_h2odata)<-c('ID','Mean AUC', 'SD AUC')
    

#  
Orth_elmdata=data.frame()
BelmOrthMatrixInfo=list()
positionOelmclass=c(231,233,235,237,239)
positionOelmpred=c(232,234,236,238,240)
OrthelmValues<-data.frame()
MatrixInfo<-NULL
for(z in 1:5){
  OrthelmROC<-roc(as.factor(FinalData[[positionOelmclass[z]]]),FinalData[[positionOelmpred[z]]],direction = "<")
  OrthelmValues<-rbind(OrthelmValues,data.frame(MMIDShort[x+1],auc(OrthelmROC)))
  
  Type=NULL
  CurrentInfo<-list(Posclass=which(FinalData[[positionOelmclass[z]]]==1),NegClass =which(FinalData[[positionOelmclass[z]]]==0),
                    PosPred =which(FinalData[[positionOelmpred[z]]]>=0.5),NegPred=which(FinalData[[positionOelmpred[z]]]<0.5))
  
  for(y in 1:length(FinalData[[positionOelmclass[z]]])) {
    Type=c(Type,(if((y %in% CurrentInfo$PosPred)==TRUE){if((y %in% CurrentInfo$Posclass)==TRUE){4}else{3}}else{if((y %in% CurrentInfo$NegClass)==TRUE){2}else{1}}))
    
  }
  MatrixInfo[z]=as.data.frame(Type)
}
for(p in 1:length(MatrixInfo)){
  BelmOrthMatrixInfo[p]=(as.data.frame(MatrixInfo[[p]]))
}
Orth_elmdata<-rbind(Orth_elmdata,data.frame('elmorth',mean(OrthelmValues[,2]),sd(OrthelmValues[,2])))

names(OrthelmValues)<-c('ID','AUC')
names(Orth_elmdata)<-c('ID','Mean AUC', 'SD AUC')

#AllMMelm

######average the elm MM preds####
x=10
positionelmAMMclass=c(1,23,45,67,89)+(2*x)
positionelmAMMpred=c(2,24,46,68,90)+(2*x)
MatrixInfo<-NULL
AllMMelmAVG<-list()
AllMMelmClass<-list()

MMNum=c(0,1,2,3,4,5,6,7,8,9)

for(z in 1:5){
  AllMMelmtempAVG=vector()
  AllMMelmTClass=vector()
  elmPredsraw=as.vector(FinalData[[positionelmAMMpred[z]]])
  elmClassraw=as.vector(FinalData[[positionelmAMMclass[z]]])
  for (W in 1:(length(elmPredsraw)*0.1)) {
    
    tempavg=mean(elmPredsraw[W+(MMNum*0.1*length(elmPredsraw))]) 
    AllMMelmtempAVG=c(AllMMelmtempAVG,tempavg)
    AllMMelmTClass=c(AllMMelmTClass,elmClassraw[W])
  }
  AllMMelmAVG[[z]]<-AllMMelmtempAVG
  AllMMelmClass[[z]]<-AllMMelmTClass
  }


AllMMelmMatrixInfo=list()
AllMMelmdata=data.frame()


allelmValues<-data.frame()
MatrixInfo<-NULL
for(z in 1:5){
  
   allelmROC<-roc(response=as.factor(AllMMelmClass[[z]]),predictor =AllMMelmAVG[[z]],direction = "<")
   allelmValues<-rbind(allelmValues,data.frame(MMIDShort[x+1],auc(allelmROC)))
  
  Type=NULL
  CurrentInfo<-list(Posclass=which(AllMMelmClass[[z]]==1),NegClass =which(AllMMelmClass[[z]]==0),
                    PosPred =which(AllMMelmAVG[[z]]>=0.5),NegPred=which(AllMMelmAVG[[z]]<0.5))
  
  for(y in 1:length(AllMMelmClass[[z]])) {
    #4=true pos, , 3=false pos, 2=true neg, 1= false neg
    Type=c(Type,(if((y %in% CurrentInfo$PosPred)==TRUE){if((y %in% CurrentInfo$Posclass)==TRUE){4}else{3}}else{if((y %in% CurrentInfo$NegClass)==TRUE){2}else{1}}))
    
  }
  MatrixInfo[z]=as.data.frame(Type)
}
for(p in 1:length(MatrixInfo)){
  AllMMelmMatrixInfo[p]=(as.data.frame(MatrixInfo[[p]]))
}
AllMMelmdata<-rbind(AllMMelmdata,data.frame(paste('elmall',MMIDShort[x+1]),mean(allelmValues[,2]),sd(allelmValues[,2])))


names(allelmValues)<-c('ID','AUC')
names(AllMMelmdata)<-c('ID','Mean AUC', 'SD AUC')



######average the h2o MM preds####
x=10
positionh2oAMMclass=c(111,133,155,177,199)+(2*x)
positionh2oAMMpred=c(112,134,156,178,200)+(2*x)


MatrixInfo<-NULL
AllMMh2oAVG<-list()
AllMMh2oClass<-list()

MMNum=c(0,1,2,3,4,5,6,7,8,9)

for(z in 1:5){
  AllMMh2otempAVG=vector()
  AllMMh2oTClass=vector()
  h2oPredsraw=as.vector(FinalData[[positionh2oAMMpred[z]]])
  h2oClassraw=as.vector(FinalData[[positionh2oAMMclass[z]]])
  for (W in 1:(length(h2oPredsraw)*0.1)) {
    
    tempavg=mean(h2oPredsraw[W+(MMNum*0.1*length(h2oPredsraw))]) 
    AllMMh2otempAVG=c(AllMMh2otempAVG,tempavg)
    AllMMh2oTClass=c(AllMMh2oTClass,h2oClassraw[W])
  }
  AllMMh2oAVG[[z]]<-AllMMh2otempAVG
  AllMMh2oClass[[z]]<-AllMMh2oTClass
}



AllMMh2oMatrixInfo=list()
AllMMh2odata=data.frame()


allh2oValues<-data.frame()
MatrixInfo<-NULL
for(z in 1:5){
  
  allh2oROC<-roc(response=as.factor(AllMMh2oClass[[z]]),predictor =AllMMh2oAVG[[z]],direction = "<")
  allh2oValues<-rbind(allh2oValues,data.frame(MMIDShort[x+1],auc(allh2oROC)))
  
  Type=NULL
  CurrentInfo<-list(Posclass=which(AllMMh2oClass[[z]]==1),NegClass =which(AllMMh2oClass[[z]]==0),
                    PosPred =which(AllMMh2oAVG[[z]]>=0.5),NegPred=which(AllMMh2oAVG[[z]]<0.5))
  
  for(y in 1:length(AllMMh2oClass[[z]])) {
    #4=true pos, , 3=false pos, 2=true neg, 1= false neg
    Type=c(Type,(if((y %in% CurrentInfo$PosPred)==TRUE){if((y %in% CurrentInfo$Posclass)==TRUE){4}else{3}}else{if((y %in% CurrentInfo$NegClass)==TRUE){2}else{1}}))
    
  }
  MatrixInfo[z]=as.data.frame(Type)
}
for(p in 1:length(MatrixInfo)){
  AllMMh2oMatrixInfo[p]=(as.data.frame(MatrixInfo[[p]]))
}
AllMMh2odata<-rbind(AllMMh2odata,data.frame(paste('h2oall',MMIDShort[x+1]),mean(allh2oValues[,2]),sd(allh2oValues[,2])))


names(allh2oValues)<-c('ID','AUC')
names(AllMMh2odata)<-c('ID','Mean AUC', 'SD AUC')
#####Summary####





outputsummary<-rbind(B_elmdata,AllMMelmdata,B_h2odata,AllMMh2odata,Orth_elmdata,Orth_h2odata)


write.csv( outputsummary,paste0(Output_Dir,"\\Output summary .csv"))
write.csv( TestSequencesKey,paste0(Output_Dir,"\\TestSequencesKey .csv"))
####Venn diagram####

Tempelm1=vector()
Tempelm2=vector()

   for (w in 1:5) {
     Tempelm1=as.vector(AllMMelmMatrixInfo[[w]]) 
     Tempelm2=c(Tempelm2,Tempelm1)
     }

Temph2o1=vector()
Temph2o2=vector()

for (w in 1:5) {
  Temph2o1=as.vector(AllMMh2oMatrixInfo[[w]]) 
  Temph2o2=c(Temph2o2,Temph2o1)
}

  BestELM=which(B_elmdata[,2]==max(B_elmdata[,2]))
  
  Tempelmmm1=vector()
  Tempelmmm2=vector()
  
  for (w in 1:5) {
    Tempelmmm1=as.vector(BelmAllMatrixInfo[[w+(5*(BestELM-1))]]) 
    Tempelmmm2=c(Tempelmmm2,Tempelmmm1)
  }
    
  
  Besth2o=which(B_h2odata[,2]==max(B_h2odata[,2]))
  
  Temph2omm1=vector()
  Temph2omm2=vector()
  
  for (w in 1:5) {
    Temph2omm1=as.vector(Bh2oAllMatrixInfo[[w+(5*(Besth2o-1))]]) 
    Temph2omm2=c(Temph2omm2,Temph2omm1)
  }
  
  
  
  
  
  Temporthelm1=vector()
  Temporthelm2=vector()
  
  for (w in 1:5) {
    Temporthelm1=as.vector(BelmOrthMatrixInfo[[w]]) 
    Temporthelm2=c(Temporthelm2,Temporthelm1)
  }
  
  Temporthh2o1=vector()
  Temporthh2o2=vector()
  
  for (w in 1:5) {
    Temporthh2o1=as.vector(Bh2oOrthMatrixInfo[[w]]) 
    Temporthh2o2=c(Temporthh2o2,Temporthh2o1)
  }
  
  
  
  

AllELM<-data.frame("Sequences" =TestSequencesKey[,1],"AllELM"=Tempelm2)
AllH2O<-data.frame("Sequences" =TestSequencesKey[,1],"AllH2O"=Temph2o2)
SpecificMMelm<-data.frame(TestSequencesKey[,1],Tempelmmm2)
names(SpecificMMelm)<-c("Sequences",MMID[BestELM])
SpecificMMh2o<-data.frame(TestSequencesKey[,1],Temph2omm2)
names(SpecificMMh2o)<-c("Sequences",MMID[Besth2o])
OrthELM<-data.frame("Sequences" =TestSequencesKey[,1],"OrthELM"=Temporthelm2)
OrthH2O<-data.frame("Sequences" =TestSequencesKey[,1],"OrthH2O"=Temporthh2o2)

Venn.data<-list("Bf+elm"=AllELM[which(AllELM[,2]==3),1],"DBNN"=AllH2O[which(AllH2O[,2]==3),1],
                SpecificMMelm[which(SpecificMMelm[,2]==3),1],
                SpecificMMh2o[which(SpecificMMh2o[,2]==3),1])
names(Venn.data)<-c("Bf+elm","DBNN",paste0("elm:",MMID[BestELM]),paste0("H2O:",MMID[Besth2o]))

AllMethod=vector()
ELMmethod=vector()
H2Omethod=vector()
DeepMethods=vector()


AllMethod=as.character(Venn.data[[1]][which(Venn.data[[1]][which(Venn.data[[1]][which(Venn.data[[1]] %in% Venn.data[[2]])]%in% Venn.data[[3]])]%in% Venn.data[[3]])])
  
  ELMmethod=as.character(Venn.data[[1]][which((Venn.data[[1]] %in% Venn.data[[3]])==TRUE)])
  DeepMethods=as.character(Venn.data[[1]][which((Venn.data[[1]] %in% Venn.data[[2]])==TRUE)])
  H2Omethod=as.character(Venn.data[[2]][which((Venn.data[[2]] %in% Venn.data[[4]])==TRUE)])


Venn_data_final=list(date(),"Bf+elm",as.character(Venn.data[[1]]),"DBNN",as.character(Venn.data[[2]]),paste0("elm:",MMID[BestELM]),as.character(Venn.data[[3]]),paste0("H2O:",MMID[Besth2o]),as.character(Venn.data[[4]]),c(" AllMethod  ",AllMethod),c(" ELMMethod  ",ELMmethod),c(" H2OMethod  ",H2Omethod),c(" DeepMethod  ",DeepMethods))
#names(Venn_data_final)<-c("Bf+elm","DBNN",paste0("elm:",MMID[BestELM]),paste0("H2O:",MMID[Besth2o]),"AllMethod","ELMMethod","H2OMethod","DeepMethod")
date()
lapply(Venn_data_final, write, "Venn Data.txt", append=TRUE, ncolumns=1000)


venn.diagram(Venn.data, fill = c("red", "green","blue","yellow"), alpha = c(0.2, 0.2,0.2,0.2), cex = 1,cat.fontface = 2,
                      lty =2, filename = "FalsePosVenn.tiff")








