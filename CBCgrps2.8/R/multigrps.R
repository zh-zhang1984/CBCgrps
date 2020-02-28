multigrps <-
  function(df,gvar,
           p.rd=3,varlist = NULL,
           skewvar=NULL,
           norm.rd=2,
           sk.rd=2,
           tabNA="no",#need to replace NaN with NA for all factors
           cat.rd=0,pnormtest=0.05,
           maxfactorlevels=30,
           minfactorlevels=10,
           sim=FALSE,
           workspace=2e5,ShowStatistic = F){
    ##group varibale must be a factor
    df[,gvar]<-as.factor(df[,gvar])
    #NaN is forced to be NA, NaN can cause problem
    df<-replace(df,is.na(df),NA)
    for(i in 1:length(levels(df[,gvar]))){
      assign(paste("g", i, sep = ""), levels(df[,gvar])[i])    
    }
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
    Table <- NULL
    #loop over variables
    for (var in varlist){
      if((class(df[,var])=="factor"|class(df[,var])=="character")&length(levels(factor(df[,var])))>maxfactorlevels){
        print(paste("the factor variable",var,
                    "contains more than",
                    maxfactorlevels,"levels","check the class of", var,
                    "or reset the maxfactorlevels",sep=' '))
        next 
      }else{
        if(class(df[,var])=="factor"|length(levels(factor(df[,var])))<=minfactorlevels){
          if(var%in%skewvar){
            stop("skewvar contains categorical variables")
          }
          if(tabNA=="no"){
            df[,var]<-factor(df[,var])
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
          tabGrp <- NULL
          nameGrp <- NULL
          for (varGrp in levels(df[,gvar])) {
            tabGrp1 <- paste(table.sub[,varGrp]," (",
                             round(per.sub[,varGrp]*100,cat.rd),
                             ")",sep = "")
            nameGrp1 <- paste(varGrp," (","n = ",table(df[,gvar])[varGrp],")",sep = "")
            tabGrp <- cbind(tabGrp,tabGrp1)
            nameGrp <- c(nameGrp,nameGrp1)
          }
          colnames(tabGrp) <- levels(df[,gvar])
          rownames(tabGrp) <- rownames(table.sub)
          table1 <- data.frame("Variables" = paste("  ",levels(df[,var]),sep = ""),
                               paste(as.data.frame(tableTol)[,"Freq"]," (",
                                     round(as.data.frame(per)[,"Freq"]*100,cat.rd),
                                     ")",sep = ""),stringsAsFactors = F)
          table1 <- cbind(table1,tabGrp,p=" ",statistic=" ")
          table1 <- apply(table1, 2,as.character)
          newline <- c(paste(var,", n (%)",sep = ""),
                       rep("",length(levels(df[,gvar]))+1),
                       ifelse(p<1*10^(-p.rd),
                              paste("< ",1*10^(-p.rd),sep = ""),
                              round(p,p.rd)),
                       ifelse(is.null(statistic),"Fisher",
                              round(statistic,3)))
          table1 <- rbind(newline,table1)
          colnames(table1)<-c("Variables",
                              paste("Total (n = ",nrow(df),")",sep = ""),
                              nameGrp,"p",
                              "statistic")
          rownames(table1) <- NULL
          Table<-rbind.data.frame(Table,table1,stringsAsFactors = F)
        }else{
          if((ad.test(df[,var])$p.value>=pnormtest&is.null(skewvar))|(!(var%in%skewvar)&!is.null(skewvar))){
            tabGrp<-NULL
            nameGrp <- NULL
            for(varGrp in levels(df[,gvar])){
              tabGrp1 <- paste(round(mean(df[df[,gvar]==varGrp,var],na.rm=T),norm.rd),
                               " \U00B1 ",
                               round(sd(df[df[,gvar]==varGrp,var],na.rm=T),norm.rd),
                               sep = "")
              nameGrp1 <- paste(varGrp," (","n = ",table(df[,gvar])[varGrp],")",sep = "")
              tabGrp <- cbind(tabGrp,tabGrp1)
              nameGrp <- c(nameGrp,nameGrp1)
            }	 
            colnames(tabGrp) <- levels(df[,gvar])
            rownames(tabGrp) <- var
            p<-summary(aov(df[,var]~df[,gvar]))[[1]][1,"Pr(>F)"]
            statistic <- summary(aov(df[,var]~df[,gvar]))[[1]][1,"F value"]
            table1 <- data.frame(paste(var,", Mean"," \U00B1 ","SD",sep = ""),
                                 paste(round(mean(df[,var],na.rm=T),norm.rd),
                                       " \U00B1 ",
                                       round(sd(df[,var],na.rm=T),norm.rd),sep = ""),
                                 tabGrp,
                                 ifelse(p<1*10^(-p.rd),
                                        paste("< ",1*10^(-p.rd),sep = ""),
                                        round(p,p.rd)),
                                 round(statistic,3),stringsAsFactors = F)
            colnames(table1)<-c("Variables",
                                paste("Total (n = ",nrow(df),")",sep = ""),
                                nameGrp,"p",
                                "statistic")
            rownames(table1) <- NULL
            Table<-rbind.data.frame(Table,table1,stringsAsFactors = F)
          }else{
            median<-as.numeric(summary(df[,var])[3])
            IQR1<-as.numeric(summary(df[,var])[2])
            IQR3<-as.numeric(summary(df[,var])[5])
            tabGrp<-NULL
            nameGrp <- NULL
            for(varGrp in levels(df[,gvar])){
              tabGrp1 <- paste(round(summary(df[df[,gvar]==varGrp,var])[3],sk.rd),
                               " (",
                               round(summary(df[df[,gvar]==varGrp,var])[2],sk.rd),", ",
                               round(summary(df[df[,gvar]==varGrp,var])[5],sk.rd),")",
                               sep = "")
              nameGrp1 <- paste(varGrp," (","n = ",table(df[,gvar])[varGrp],")",sep = "")
              tabGrp <- cbind(tabGrp,tabGrp1)
              nameGrp <- c(nameGrp,nameGrp1)
            }		
            p<-kruskal.test(df[,var]~df[,gvar])$p.value
            statistic <- kruskal.test(df[,var]~df[,gvar])$statistic
            table1 <- data.frame(paste(var,", Median"," (IQR)",sep = ""),
                                 paste(round(median,sk.rd)," (",round(IQR1,sk.rd),", ",
                                       round(IQR3,sk.rd),")",sep = ""),
                                 tabGrp,
                                 ifelse(p<1*10^(-p.rd),
                                        paste("< ",1*10^(-p.rd),sep = ""),
                                        round(p,p.rd)),
                                 round(statistic,3),stringsAsFactors = F)
            colnames(table1)<-c("Variables",
                                paste("Total (n = ",nrow(df),")",sep = ""),
                                nameGrp,"p",
                                "statistic")
            rownames(table1) <- NULL
            Table<-rbind.data.frame(Table,table1,stringsAsFactors = F)
          }}
      }
    }
    if(!ShowStatistic){
     Table <- Table[,!colnames(Table) %in% "statistic"]
    }
    Table <- rbind(colnames(Table),Table)
    colnames(Table) <- NULL
    return(Table)       	      	     	
  }
