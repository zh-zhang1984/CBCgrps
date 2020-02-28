twogrps <-
  function(df,gvar,varlist = NULL,
           p.rd=3, 
           skewvar=NULL,
           norm.rd=2,
           sk.rd=2,
           tabNA="no",#need to replace NaN with NA for all factors
           cat.rd=0, pnormtest=0.05,
           maxfactorlevels=30,
           minfactorlevels=10,
           sim = FALSE,#to use simulated p value
           workspace=2e5,ShowStatistic = F,ExtractP = 0.05){
    ##group varibale must be a factor
    df[,gvar]<-as.factor(df[,gvar])
    if(length(table(df[,gvar]))>2){
      stop("The gvar contains more than two levels, please  recheck or consider using multigrps() function")
    }
    #NaN is forced to be NA, NaN can cause problem
    df<-replace(df,is.na(df),NA)
    g1<-levels(df[,gvar])[1]
    g2<-levels(df[,gvar])[2]
    if(is.null(varlist)){
      varlist<-names(df)[!names(df)%in%gvar]
    }else if(sum(!(varlist%in%names(df)))>0){
      stop("varlist contains variables not in the data frame")
    } else {
      varlist = varlist
    }
    if(sum(!(skewvar%in%varlist))>0){
      stop("skewvar contains variables not in the data frame or varlist")
    }
    Table <- NULL;
    VarExtract <- NULL;#extract variable names with p value less than ExtractP
    #loop over variables
    for (var in varlist){
      if((class(df[,var])=="factor"|class(df[,var])=="character")&length(levels(factor(df[,var]))) > maxfactorlevels){
        print(paste("The factor/character variable", var,
                    "contains more than",
                    maxfactorlevels,"levels,","check the class of", var,
                    "or reset the maxfactorlevels",sep=' '))
        next 
      }else{	
        if(class(df[,var]) == "factor"|class(df[,var])=="character"|length(levels(factor(df[,var]))) <= minfactorlevels){
          if(var %in% skewvar){
            stop("skewvar contains categorical variables")
          }
          if(tabNA =="no"){
            df[,var] <- factor(df[,var])
          }else{
            df[,var]<-factor(df[,var],exclude = NULL)
          }
          tableTol<-table(df[,var],useNA=tabNA)
          per<-prop.table(tableTol)
          table.sub<-table(df[,var],df[,gvar],useNA=tabNA)
          per.sub<-prop.table(table.sub,2)
          if(nrow(table.sub)==1){
            p =1;statistic=NULL
          } else{
            p<-tryCatch({#using fisher's test when scarce data
              chisq.test(table.sub)$p.value
            }, warning = function(w) {
              fisher.test(table.sub,
                          workspace = workspace,
                          simulate.p.value = sim)$p.value
            })
            statistic<-tryCatch({#using fisher's test when scarce data
              chisq.test(table.sub)$statistic
            }, warning = function(w) {
              NULL
            })
          }
          names(statistic) <- NULL
          table1 <- data.frame("Variable"=paste("  ",levels(df[,var]),sep = ""),
                               paste(as.data.frame(tableTol)[,"Freq"]," (",
                                     round(as.data.frame(per)[,"Freq"]*100,cat.rd),
                                     ")",sep = ""),
                               paste(as.data.frame.matrix(table.sub)[,g1]," (",
                                     round(as.data.frame.matrix(per.sub)[,g1]*100,cat.rd),
                                     ")",sep = ""),
                               paste(as.data.frame.matrix(table.sub)[,g2]," (",
                                     round(as.data.frame.matrix(per.sub)[,g2]*100,cat.rd),
                                     ")",sep = ""),
                               p = "",
                               statistic = "",stringsAsFactors = F)
          newline <- c(paste(var,", n (%)",sep = ""),
                       rep("",3),
                       ifelse(p<1*10^(-p.rd),
                              paste("< ",1*10^(-p.rd),sep = ""),
                              round(p,p.rd)),
                       ifelse(is.null(statistic),"Fisher",
                              round(statistic,3)))
          table1 <- rbind(newline,table1)
          colnames(table1)<-c("Variables",
                              paste("Total (n = ",nrow(df),")",sep = ""),
                              paste(g1," (n = ",nrow(df[df[,gvar]==g1,]),")",sep = ""),
                              paste(g2," (n = ",nrow(df[df[,gvar]==g2,]),")",sep = ""),"p",
                              "statistic")
          rownames(table1) <- NULL
          Table<-rbind.data.frame(Table,table1,stringsAsFactors = F)
          if(p < ExtractP){
            VarExtract <- c(VarExtract,var)
          }
        }else{
          if((ad.test(df[,var])$p.value>=pnormtest&is.null(skewvar))|(!(var%in%skewvar)&!is.null(skewvar))){
            mean<-round(mean(df[,var],na.rm=T),norm.rd)
            sd<-round(sd(df[,var],na.rm=T),norm.rd)
            mean.1<-round(mean(df[df[,gvar]==g1,var],na.rm=T),norm.rd)
            sd.1<-round(sd(df[df[,gvar]==g1,var],na.rm=T),norm.rd)
            mean.2<-round(mean(df[df[,gvar]==g2,var],na.rm=T),norm.rd)
            sd.2<-round(sd(df[df[,gvar]==g2,var],na.rm=T),norm.rd)
            p<-t.test(df[,var]~df[,gvar])$p.value
            statistic <- t.test(df[,var]~df[,gvar])$statistic
            table1 <- data.frame("Variable"= paste(var,", Mean"," \U00B1 ","SD",sep = ""),
                                 paste(mean," \U00B1 ",sd,sep=""),
                                 paste(mean.1," \U00B1 ",sd.1,sep=""),
                                 paste(mean.2," \U00B1 ",sd.2,sep=""),
                                 p=ifelse(p<1*10^(-p.rd),
                                          paste("< ",1*10^(-p.rd),sep = ""),
                                          round(p,p.rd)),
                                 statistic = round(statistic,3),stringsAsFactors = F)
            colnames(table1)<-c("Variables",
                                paste("Total (n = ",nrow(df),")",sep = ""),
                                paste(g1," (n = ",nrow(df[df[,gvar]==g1,]),")",sep = ""),
                                paste(g2," (n = ",nrow(df[df[,gvar]==g2,]),")",sep = ""),"p",
                                "statistic")
            rownames(table1) <- NULL
            Table <- rbind.data.frame(Table,table1,stringsAsFactors = F)
            if(p < ExtractP){
              VarExtract <- c(VarExtract,var)
            }
          }else{
            median<-as.numeric(summary(df[,var])[3])
            IQR1<-as.numeric(summary(df[,var])[2])
            IQR3<-as.numeric(summary(df[,var])[5])
            median.1<-as.numeric(summary(df[df[,gvar]==g1,var])[3])
            IQR1.1<-as.numeric(summary(df[df[,gvar]==g1,var])[2])
            IQR3.1<-as.numeric(summary(df[df[,gvar]==g1,var])[5])
            median.2<-as.numeric(summary(df[df[,gvar]==g2,var])[3])
            IQR1.2<-as.numeric(summary(df[df[,gvar]==g2,var])[2])
            IQR3.2<-as.numeric(summary(df[df[,gvar]==g2,var])[5])
            p<-wilcox.test(df[,var]~df[,gvar])$p.value
            statistic <- wilcox.test(df[,var]~df[,gvar])$statistic
            table1<-data.frame("Variable"= paste(var,", Median"," (IQR)",sep = ""),
                               paste(round(median,sk.rd)," (",round(IQR1,sk.rd),", ",round(IQR3,sk.rd),")",sep = ""),
                               paste(round(median.1,sk.rd)," (",round(IQR1.1,sk.rd),", ",round(IQR3.1,sk.rd),")",sep = ""),
                               paste(round(median.2,sk.rd)," (",round(IQR1.2,sk.rd),", ",round(IQR3.2,sk.rd),")",sep = ""),
                               p=ifelse(p<1*10^(-p.rd),
                                        paste("< ",1*10^(-p.rd),sep = ""),
                                        round(p,p.rd)),
                               statistic = round(statistic,3),stringsAsFactors = F)
            colnames(table1)<-c("Variables",
                                paste("Total (n = ",nrow(df),")",sep = ""),
                                paste(g1," (n = ",nrow(df[df[,gvar]==g1,]),")",sep = ""),
                                paste(g2," (n = ",nrow(df[df[,gvar]==g2,]),")",sep = ""),"p",
                                "statistic")
            rownames(table1) <- NULL
            Table <- rbind.data.frame(Table,table1,stringsAsFactors = F)
            if(p < ExtractP){
              VarExtract <- c(VarExtract,var)
            }
          }
        }
      }
    }
    if(!ShowStatistic){
      Table$statistic <- NULL
    }
    Table <- rbind(colnames(Table),Table)
    colnames(Table) <- NULL
    return(list(Table=Table,VarExtract = VarExtract))       	      	
    #the end of the function      	
  }
