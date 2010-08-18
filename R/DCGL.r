.packageName <- "DCGL"
  #####################################################################################################
  ## a gene filtering function: Genes which have a Between-Experiment Mean  Expression Signal (BEMES) lower than the median of BEMES's 
  ## of all genes will be filtered out. 
  ## 'exprs' is the expression data for two conditions.
  # output: A data frame or matrix with a reduced number of rows.
  #####################################################################################################
"expressionBasedfilter" <- function(exprs) {
 	all.mean <- apply(exprs,1,mean)
  	median <- median(all.mean)
  	exprs.filter <- exprs[all.mean>median,]
  	return(exprs.filter)
}

  #################################################################################################### 
  ## a gene filtering function: Those genes not significantly more variable than the median gene are filtered out.
  ## 'exprs' is the expression data for two conditions.
  ## 'p' is the probability  cut-off of Chi Squared distribution.
  # output: A data frame or matrix with a reduced number of rows.
  ####################################################################################################
"varianceBasedfilter" <- function(exprs,p) {
 	n <- ncol(exprs)
 	Var.i <- apply(exprs,1,var)
 	Var.median <- median(Var.i)
 	quantity <- (n-1)*Var.i/Var.median
 	degree <- ncol(exprs)-1
 	prob <- pchisq(quantity, degree, ncp=0, lower.tail = F, log.p = FALSE)
 	exprs.filter <- exprs[prob<=p,]
 	return(exprs.filter)
}

  ###############################################################################################
  ## a link filtering function: user can set the threshold of correlation based on the relationship (a curve) between correlation and clustering coefficient.
  ## exprs: GEM, each column is a condition. Affy-like intensity (log 2 transformed) assumed.
  ###############################################################################################
"systematicLinkfilter" <-
function(exprs) {
N<- nrow(exprs);
exprs=as.matrix(exprs);
r=cor(t(exprs),method="pearson",use="pairwise.complete.obs");
CC0<- vector();
CCmean<- vector();
for(l in 1:100) {

Rth<- 0.01*l;
K <- rep(0, dim(r)[1]);
names(K) <- rownames(r);
CC <- K;
for( i in 1:dim(r)[1]){
  neighbor <- rownames(r)[i];
  for( j in 1:dim(r)[2] ){
    if( r[i,j] >= Rth ){
      K[i] = K[i]+1;
      neighbor <- c(neighbor, colnames(r)[j]);
    }
  }
  if( length(neighbor) > 2 ){
    for( m in 2:length(neighbor)){
      for( n in 2:length(neighbor) ){
        if( r[neighbor[m], neighbor[n]] >= Rth){
          CC[i] = CC[i] + 1;
        }
      }
    }
    CC[i] = CC[i]/(K[i]*(K[i]-1));
  }
 }
 K2<- 0;
 KK<- 0;
 CCC<- 0;
 for(j in 1:100) {
 K2<- K[j]*K[j]+K2;
 KK<- KK+K[j]; 
 CCC<- CC[i]+CCC;
 }  
CCCmean<- CCC/N;
K2mean<- K2/N;
Kmean<- KK/N;
C0<- (K2mean-Kmean)^2/(Kmean^3)/N;
CC0<- rbind(CC0,C0);
CCmean<- rbind(CCmean,CCCmean);
} 
CC_C0<- CCmean-CC0;    
cor <- seq(0.01,1,0.01)
#plot(x,CC_C0,xlab="correlation threshold",ylab="C-C0");
#lines(x,CC_C0);
C_r <- cbind(cor,CC_C0)
colnames(C_r)=c("cor","C_C0")
rownames(C_r)=c(1:nrow(C_r))
return(C_r)
}

  ######################################################################################
  ## a link filtering function: Gene links with q-values of coexpression values in either of two conditions higher than the cutoff (qth) are reserved.
  ## exprs.1 a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.  
  ## exprs.2 a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.  
  ## qth the cutoff of q-value; must be within [0,1].
  # output: A list with two components of lists: one lists the two rths (thresholds of correlation values) corresponding to the two conditions, the other lists the two matrices of filtered pairwise correlation values.
  ######################################################################################
"qLinkfilter" <-function(exprs.1,exprs.2,qth) {
  	# m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
  	n.1 <- ncol(exprs.1)
  	n.2 <- ncol(exprs.2)
  
  	degree.1 <- n.1-2
  	degree.2 <- n.2-2
  
  	genes <- rownames(exprs.1)
  	exprs.1 <- as.matrix(exprs.1)
  	exprs.2 <- as.matrix(exprs.2)
  	cor.pairwise.1 <- cor(t(exprs.1),method="pearson",use="pairwise.complete.obs")
  	cor.pairwise.2 <- cor(t(exprs.2),method="pearson",use="pairwise.complete.obs")
  	t.1 <- cor.pairwise.1*sqrt(degree.1)/sqrt(1-cor.pairwise.1*cor.pairwise.1)
  	t.2 <- cor.pairwise.2*sqrt(degree.2)/sqrt(1-cor.pairwise.2*cor.pairwise.2)
  
  	p0.1 <- 2*pt(-abs(t.1), degree.1, lower.tail = TRUE, log.p = FALSE)
  	p0.2 <- 2*pt(-abs(t.2), degree.2, lower.tail = TRUE, log.p = FALSE)
  	diag(p0.1) <- NA
  	diag(p0.2) <- NA

  	q.1 <- p.adjust(p0.1,method="BH")
  	q.2 <- p.adjust(p0.2,method="BH")
  	diag(q.1) <- 1
  	diag(q.2) <- 1
  	
  	rth.1 <- abs(cor.pairwise.1[which.min(abs(q.1-qth))])  
   	rth.2 <- abs(cor.pairwise.2[which.min(abs(q.2-qth))])
  	cor.pairwise.1 [q.1>=qth & q.2>=qth] <- 0
  	cor.pairwise.2 [q.1>=qth & q.2>=qth] <- 0
  
  	cor.filtered <- list(rth = list(rth.1 = rth.1, rth.2 = rth.2), cor.filtered = list(cor.filtered.1 = cor.pairwise.1, cor.filtered.2 = cor.pairwise.2))
		return(cor.filtered)
}

  ###############################################################################################################
  ## Select part of the correlation pairs, the max(abs(cor.1),abs(cor.2))
  ## exprs.1 a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.  
  ## exprs.2 a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.  
  ## percent percent of links to be reserved. 
  # output: A list with two components of lists: one lists the rth (thresholds of correlation values) for both conditions, the other lists the two matrices of filtered pairwise correlation values.
  ###############################################################################################################
"percentLinkfilter" <-function(exprs.1,exprs.2,percent) {
  	# m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
  	n.1 <- ncol(exprs.1)
  	n.2 <- ncol(exprs.2)
  
  	degree.1 <- n.1-2
  	degree.2 <- n.2-2
  
  	genes <- rownames(exprs.1)
  	exprs.1 <- as.matrix(exprs.1)
  	exprs.2 <- as.matrix(exprs.2)
  	cor.pairwise.1 <- cor(t(exprs.1),method="pearson",use="pairwise.complete.obs")
  	cor.pairwise.2 <- cor(t(exprs.2),method="pearson",use="pairwise.complete.obs")
  	diag(cor.pairwise.1) <- 0
  	diag(cor.pairwise.1) <- 0
  	cor.1.2.vector <- cbind(cor.pairwise.1[lower.tri(cor.pairwise.1, diag = FALSE)],cor.pairwise.2[lower.tri(cor.pairwise.2, diag = FALSE)])
  	cor.1.2.max <- apply(abs(cor.1.2.vector),1,max)  #max of the two cor
  	cor.1.2.max.sort <- sort(cor.1.2.max,decreasing = TRUE);
  	# cor.1.2.max.sort <- as.matrix(cor.1.2.max.sort)
  	Rth <- cor.1.2.max.sort[as.integer(length(cor.1.2.max.sort)*percent)];
  
  	cor.filtered.1 <- cor.pairwise.1
  	cor.filtered.2 <- cor.pairwise.2
  	cor.filtered.1[abs(cor.pairwise.1)<=Rth & abs(cor.pairwise.2)<=Rth] <- 0
  	cor.filtered.2[abs(cor.pairwise.1)<=Rth & abs(cor.pairwise.2)<=Rth] <- 0
  	cor.filtered <- list(rth=list(rth.1=Rth,rth.2=Rth),cor.filtered=list(cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2))
		return(cor.filtered)
}

#################################################################################################
# the LFC model.
#################################################################################################
"LFC" <- function(exprs,nbins=20,p=0.1,sign,figname = 'LFC.jpeg') {
 if(sign=='same'){
		exprs.min <- apply(abs(exprs),1,min)
  	exprs.max <- apply(abs(exprs),1,max)
  	exprs.diff <- exprs.max/exprs.min
  	delinks=list(exprs.max=exprs.max,exprs.diff=exprs.diff)
    }else {
	 			exprs.min <- apply(abs(exprs),1,min)
  			exprs.max <- apply(abs(exprs),1,max)
  			exprs.diff <- exprs.max/exprs.min
  			a <- exprs.max
	 			delinks=list(exprs.max=exprs.diff,exprs.diff=a)
	 	}
  exprs.max=delinks$exprs.max
  exprs.diff=delinks$exprs.diff
  n.tol = length(exprs.max)
	num.bin = ceiling(n.tol/nbins)
	exprs.max.sorted = sort(exprs.max)
	steps = min(exprs.max)
	for (i in 1:nbins) {
		if (i == nbins) {
			steps = c(steps, exprs.max.sorted[n.tol]+1)
		} else {
			steps = c(steps, exprs.max.sorted[i*num.bin])
		}	
	}
	exprs.bin<-rep(1,n.tol);
	for (i in 1:nbins) {
		exprs.bin[exprs.max>=steps[i] & exprs.max<steps[i+1]]<-i
	}
	
	steps.x<-(steps[1:(length(steps)-1)]+steps[2:length(steps)])/2
	exprs.x.sub<-exprs.y.sub<-NULL
	for (i in 1:nbins) {
		if (sum(exprs.bin==i)>0) {
			exprs.diff.sub<-exprs.diff[exprs.bin==i]
			diff.rank<-sort(exprs.diff.sub,decreasing=T,index.return=T)$ix
			exprs.x.sub<-c(exprs.x.sub,exprs.max[exprs.bin==i][diff.rank[ceiling(length(exprs.diff.sub)*p)]])
			exprs.y.sub<-c(exprs.y.sub,exprs.diff[exprs.bin==i][diff.rank[ceiling(length(exprs.diff.sub)*p)]])
			
		}
	}
  jpeg(figname)
  		if(sign=='same'){
  				plot((exprs.max),log2(exprs.diff),pch=16,col='darkgray',xlab='max(corr)',ylab='log2(corr ratio)',main='LFC')
  				points((exprs.x.sub),log2(exprs.y.sub),pch=16)
  		}
  	else {
  				plot(log2(exprs.max),(exprs.diff),pch=16,col='darkgray',xlab='log2(corr ratio)',ylab='max(corr)',main='LFC')
  				points(log2(exprs.x.sub),(exprs.y.sub),pch=16)
  		}
	
	mm<-glm(y~x,data=data.frame(y=exprs.y.sub,x=1/exprs.x.sub))
	x<-seq(min(exprs.max),max(exprs.max),length.out=500)
	y<-predict(mm,newdata=data.frame(x=1/x))
	exprs.diff.threshold<-predict(mm,newdata=data.frame(x=1/exprs.max))
	delink<- (exprs.diff>exprs.diff.threshold)
	
		if(sign=='same'){
			points((exprs.max[delink]),log2(exprs.diff[delink]),pch=16,col='lightblue')
			lines((x[!is.na(y)]),log2(y[!is.na(y)]),col='red')
			points((exprs.x.sub),log2(exprs.y.sub),pch=16)
		} else {
			points(log2(exprs.max[delink]),(exprs.diff[delink]),pch=16,col='lightblue')
			exprs.max <- sort(exprs.max);
			exprs.diff.threshold <- sort(exprs.diff.threshold);
			lines(log2(exprs.max),(exprs.diff.threshold),col='red')
			points(log2(exprs.x.sub),(exprs.y.sub),pch=16)
		}
	dev.off()
	delink
  }

#####################################################################################
## 'exprs.1' a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.  
# 'exprs.2' a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.  
# 'N' ramdom sampling times to form the NULL distribution. If N>0 ramdom sampling is used to form the NULL distribution. 
#####################################################################################

"DCp" <- function(exprs.1,exprs.2,method=Linkfilter.methods,cutoff,N=0) {  
  	#Linkfilter.method=c('rth','qth','percent')
  	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
  	genes = rownames(exprs.1)
  	if (nrow(exprs.1)!=nrow(exprs.2)) stop("invalid row lengths") else {
	 cor.filtered = switch(method,
		rth = {
  			genes <- rownames(exprs.1)
  			exprs.1 <- as.matrix(exprs.1)
  			exprs.2 <- as.matrix(exprs.2)
  			cor.pairwise.1 <- cor(t(exprs.1),method="pearson",use="pairwise.complete.obs")
  			cor.pairwise.2 <- cor(t(exprs.2),method="pearson",use="pairwise.complete.obs")
  			diag(cor.pairwise.1) <- 0
  			diag(cor.pairwise.2) <- 0
    		cor.filtered.1 <- cor.pairwise.1
  			cor.filtered.2 <- cor.pairwise.2
  			cor.filtered.1[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
  			cor.filtered.2[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
   			list(rth=list(rth.1=cutoff, rth.2=cutoff),cor.filtered=list(cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2))
		
		},
	 	qth =  	qLinkfilter(exprs.1,exprs.2,cutoff),
		percent = percentLinkfilter(exprs.1,exprs.2,cutoff) )

	
	rth = cor.filtered$rth
	rth.1 = rth$rth.1
	rth.2 = rth$rth.2
	cor.filtered.1 = cor.filtered$cor.filtered$cor.filtered.1
	cor.filtered.2 = cor.filtered$cor.filtered$cor.filtered.2
	rownames(cor.filtered.1)=rownames(cor.filtered.2)=colnames(cor.filtered.1)=colnames(cor.filtered.2)<-rownames(exprs.1)
	rm(cor.filtered)  
  #####################################################################################################################################################
  ## For one gene, there are n-1 correlation value pairs also q value pairs(From the two conditions). For a correlation pair for this gene,
  ## if one of the q values is less than the threshold, the pair will be retained. Then there are 'number.i.uniq' unique pairs retained that is there 
  ## are two vectors of correlation values. 
  ## Then a length normalized Euclidean Distance for the two vectors will be calculated (LNED). 
  ##################################################################################################################################################### 
  	squares = (cor.filtered.1-cor.filtered.2)^2
  	number_uniq = apply(cor.filtered.1!=0 | cor.filtered.2!=0,1,sum)
  	ss = apply(squares,1,sum)
  	LNED.result = as.vector(matrix(NA,m,1))
  	LNED.result[number_uniq!=0] = sqrt(ss[number_uniq!=0])/sqrt(number_uniq[number_uniq!=0])
  	

 
  
 ########################################################################################################################
 ## Disturb the sample labels for the two conditions and re-assign the samples to two datasets,then calculate the 'dC0' for 
 ## N times and then pool all the dC0 together to construct a 'NULL' distribution.
 #########################################################################################################################
		if(N>0){
		dC0 <- vector()
		exprs <- cbind(exprs.1,exprs.2)
		n.1 = ncol(exprs.1)
		n.2 = ncol(exprs.2)
		for(j in 1:N) {
			seq <- sample(n.1+n.2)
			exprs_1 <- exprs[,seq[1:n.1]];
			exprs_2 <- exprs[,seq[(n.1+1):(n.1+n.2)]];
			cor_filtered = switch(method,
		rth = {
  			genes <- rownames(exprs.1)
  			exprs.1 <- as.matrix(exprs.1)
  			exprs.2 <- as.matrix(exprs.2)
  			cor.pairwise.1 <- cor(t(exprs.1),method="pearson",use="pairwise.complete.obs")
  			cor.pairwise.2 <- cor(t(exprs.2),method="pearson",use="pairwise.complete.obs")
  			diag(cor.pairwise.1) <- 0
  			diag(cor.pairwise.2) <- 0
    		cor.filtered.1 <- cor.pairwise.1
  			cor.filtered.2 <- cor.pairwise.2
  			cor.filtered.1[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
  			cor.filtered.2[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
   			list(rth=list(rth.1=cutoff, rth.2=cutoff),cor.filtered=list(cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2))
		
		},
	 	qth =  	qLinkfilter(exprs.1,exprs.2,cutoff),
		percent = percentLinkfilter(exprs.1,exprs.2,cutoff) )
  			cor_filtered_1 = cor_filtered$cor.filtered$cor.filtered.1
  			cor_filtered_2 = cor_filtered$cor.filtered$cor.filtered.2
  			rm(cor_filtered)
  			squares = (cor_filtered_1-cor_filtered_2)^2
  			number_uniq_0 = apply(cor_filtered_1!=0 | cor_filtered_2!=0,1,sum)
  			ss = apply(squares,1,sum)
  			LNED.result.0 = numeric(number_uniq_0)
  			LNED.result.0[number_uniq_0!=0] = sqrt(ss[number_uniq_0!=0])/sqrt(number_uniq[number_uniq_0!=0])

  			dC0 = c(dC0,LNED.result.0)
   		}


		p.value <- vector()
		for(k in 1:m){
			p <- sum(LNED.result[k]<dC0)/(N*m);
			p.value <- c(p.value,p);
		}
  		FWER <- p.adjust(p.value,method="BH")
  		Result <- cbind(LNED.result,number_uniq,p.value,FWER);
  		row.names(Result) <- genes;
  		colnames(Result) <- paste(c("dC","length","p.value","FWER"));
  		# dec.idx <-sort(as.numeric(Result[,4]),method = "quick", index.return=TRUE, decreasing=FALSE)$ix
  		# Result <- Result[dec.idx,]
	}
	else { 
	##############################################
  ## Don't use the permutation step.
  ###############################################
  
  		Result<- cbind(LNED.result,number_uniq);
  		row.names(Result) <- genes;
  		colnames(Result) <- c("dC","length");
  		# dec.idx <-sort(as.numeric(Result[,1]),method = "quick", index.return=TRUE,decreasing=TRUE)$ix
  		# Result <- as.matrix(Result[dec.idx,])
  	} 
  	return(Result)
	}
}

"DCe" <-
function(exprs.1,exprs.2,method=Linkfilter.methods,cutoff,nbins=20,p=0.1,figname = c('LFC.s.jpeg','LFC.d.jpeg')) {
 	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
	genes = rownames(exprs.1)
  	if (nrow(exprs.1)!=nrow(exprs.2)) stop("invalid row lengths") else {
		cor.filtered = switch(method,
		rth = {
  			genes <- rownames(exprs.1)
  			exprs.1 <- as.matrix(exprs.1)
  			exprs.2 <- as.matrix(exprs.2)
  			cor.pairwise.1 <- cor(t(exprs.1),method="pearson",use="pairwise.complete.obs")
  			cor.pairwise.2 <- cor(t(exprs.2),method="pearson",use="pairwise.complete.obs")
  			diag(cor.pairwise.1) <- 0
  			diag(cor.pairwise.2) <- 0
    		cor.filtered.1 <- cor.pairwise.1
  			cor.filtered.2 <- cor.pairwise.2
  			cor.filtered.1[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
  			cor.filtered.2[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
   			list(rth=list(rth.1=cutoff, rth.2=cutoff),cor.filtered=list(cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2))
		
		},
	 	qth =  	qLinkfilter(exprs.1,exprs.2,cutoff),
		percent = percentLinkfilter(exprs.1,exprs.2,cutoff) )

	rth = cor.filtered$rth
	cor.filtered.1 = cor.filtered$cor.filtered$cor.filtered.1
	cor.filtered.2 = cor.filtered$cor.filtered$cor.filtered.2
	cor.filtered.1[lower.tri(cor.filtered.1, diag = TRUE)] = 0
	cor.filtered.2[lower.tri(cor.filtered.2, diag = TRUE)] = 0
	rownames(cor.filtered.1)=rownames(cor.filtered.2)=colnames(cor.filtered.1)=colnames(cor.filtered.2)<-rownames(exprs.1)
	rm(cor.filtered)
	
  	#############################################################
  	## decide three sets of correlation pairs and organize them into two-columned matrices.
  	#############################################################  	
  	
  	idx.same = (cor.filtered.1*cor.filtered.2)>0
  	idx.diff = (cor.filtered.1*cor.filtered.2)<0
  	idx.switched = (cor.filtered.1*cor.filtered.2<0)& ( abs(cor.filtered.1)>=rth$rth.1 & abs(cor.filtered.2)>=rth$rth.2 )
  	
  	cor.same = cbind(cor.filtered.1[idx.same],cor.filtered.2[idx.same])
  	cor.switched = cbind(cor.filtered.1[idx.switched],cor.filtered.2[idx.switched])  	
  	cor.diff = cbind(cor.filtered.1[idx.diff & (!idx.switched)],cor.filtered.2[idx.diff & (!idx.switched)])
 	n.switchedDCL = nrow(cor.switched)
 	
 	name.rows = matrix(rep(genes,m),m,m,byrow=T)
  	name.columns = matrix(rep(genes,m),m,m)
  	
	name.same = cbind(name.rows[idx.same],name.columns[idx.same])
	name.switched = cbind(name.rows[idx.switched],name.columns[idx.switched])
	name.diff= cbind(name.rows[idx.diff & (!idx.switched)],name.columns[idx.diff & (!idx.switched)])

	name.all = rbind(name.same,name.switched,name.diff)
  #############################################################
  ## Determine DCLs from same sign correlation pairs
  #############################################################
  	if(nrow(cor.same)>1){
		de.s = LFC(cor.same,nbins,p,sign="same",figname = figname[1])
		DCL.same = cor.same[de.s,]
		name.same = name.same[de.s,]
		n.sameDCL = nrow(DCL.same)
		DCL.same <- data.frame(name.same,DCL.same);
		colnames(DCL.same) <- c("Gene.1","Gene.2","cor.1","cor.2")
 	} else stop("only one or no same-signed pair in all!")

  #############################################################
  ## Determine DCLs from different sign correlation pairs
  #############################################################	
  	if(nrow(cor.diff)>1){
		de.d = LFC(cor.diff,nbins,p,sign="diff",figname = figname[2])
		DCL.diff = cor.diff[de.d,]
		name.diff = name.diff[de.d,]
  		n.diffDCL = nrow(DCL.diff)
  		DCL.diff <- data.frame(name.diff,DCL.diff);
		colnames(DCL.diff) <- c("Gene.1","Gene.2","cor.1","cor.2")
	} else stop("only one or no differently-signed pair in all!")

	
################################################################################################
## Determine Switched DCLs if they exist 
################################################################################################
	
	pairs=rbind(name.same,name.diff,name.switched);
	if(n.switchedDCL>0) {
		DCL.switched <- data.frame(name.switched,cor.switched);
		colnames(DCL.switched) <- c("Gene.1","Gene.2","cor.1","cor.2")
		cor.max <- apply(abs(cor.switched),1,max)
		middle <-sort(cor.max,method = "quick", index.return=TRUE,decreasing=TRUE)$ix ########
		DCL.switched<- DCL.switched[middle,]
	}

	library(igraph);
####################################
## All links
####################################
	g.all <- graph.data.frame(name.all);
	gene.all <- as.matrix(V(g.all)$name);
	de.all <- degree(g.all);  
#####################################
## DCLs
#####################################
	g <- graph.data.frame(pairs);
	gene.1 <- as.matrix(V(g)$name);
	de <- degree(g);  
######################################
##DCLs of same sign
######################################
	g.same <- graph.data.frame(name.same);
	g.same.name <- as.matrix(V(g.same)$name);
	degree.same <- as.matrix(degree(g.same));

########################################
## DCLs of different sign
########################################

	g.diff <- graph.data.frame(name.diff);
	g.diff.name <- as.matrix(V(g.diff)$name);
	degree.diff <- as.matrix(degree(g.diff));

#######################################
## DCLs of switched correlation
#######################################
	if(n.switchedDCL>0) {
		g.switch <- graph.data.frame(name.switched);
		g.switch.name <- as.matrix(V(g.switch)$name);
		degree.switch <- as.matrix(degree(g.switch));
	} else { degree.switch = matrix(0,1,1)
	DCL.switched = matrix("NULL",1,1)
	}

#######################################
## Numbers for DCLs of different type. 
#######################################

	degree.bind <- matrix(0,m,5)
	row.names(degree.bind) <- genes
	colnames(degree.bind) <- c("All.links","DC.links","DCL.same","DCL.diff","DCL.switched")

	degree.bind[gene.all,1]=de.all
	degree.bind[gene.1,2]=de
	degree.bind[g.same.name,3]=degree.same
	degree.bind[g.diff.name,4]=degree.diff
	if(n.switchedDCL>0) {
	degree.bind[g.switch.name,5]=degree.switch
	}


########################################################
## DCGs Identification
########################################################

 	prob <- nrow(pairs)/nrow(name.all)
 	p.value <- vector()
 	for(i in 1:m) {
 		x <- pbinom(degree.bind[i,2]-1, degree.bind[i,1], prob, lower.tail = F, log.p = FALSE);
 		p.value <- rbind(p.value,x);
	}
 	q.value <- p.adjust(p.value,method="bonferroni");
 
 	degree.bind <- cbind(degree.bind,p.value,q.value)
 	colnames(degree.bind) <- c("All.links","DC.links","DCL_same","DCL_diff","DCL_switch","p","q")

 	middle <-sort(as.numeric(degree.bind[,6]),method = "quick", decreasing=FALSE,index.return=TRUE)$ix 
 	DCGs <- degree.bind[middle,]
 
 #########################################################
 
 	Result <- list(DCGs=DCGs,DCL.same=DCL.same,DCL.diff=DCL.diff,DCL.switched=DCL.switched)
	return(Result)
	}
}

##############################################################################################################
## ASC
##############################################################################################################
"ASC" <-
function(exprs.1,exprs.2,method=Linkfilter.methods,cutoff) {
 	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
	genes = rownames(exprs.1)
  	if (nrow(exprs.1)!=nrow(exprs.2)) stop("invalid row lengths") else {
  		cor.filtered = switch(method,
		rth = {
  			genes <- rownames(exprs.1)
  			exprs.1 <- as.matrix(exprs.1)
  			exprs.2 <- as.matrix(exprs.2)
  			cor.pairwise.1 <- cor(t(exprs.1),method="pearson",use="pairwise.complete.obs")
  			cor.pairwise.2 <- cor(t(exprs.2),method="pearson",use="pairwise.complete.obs")
  			diag(cor.pairwise.1) <- 0
  			diag(cor.pairwise.2) <- 0
    		cor.filtered.1 <- cor.pairwise.1
  			cor.filtered.2 <- cor.pairwise.2
  			cor.filtered.1[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
  			cor.filtered.2[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
   			list(rth=list(rth.1=cutoff, rth.2=cutoff),cor.filtered=list(cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2))
		
		},
	 	qth =  	qLinkfilter(exprs.1,exprs.2,cutoff),
		percent = percentLinkfilter(exprs.1,exprs.2,cutoff) )


		rth = cor.filtered$rth
		rth.1 = rth$rth.1
		rth.2 = rth$rth.2
		cor.filtered.1 = cor.filtered$cor.filtered$cor.filtered.1
		cor.filtered.2 = cor.filtered$cor.filtered$cor.filtered.2
		rownames(cor.filtered.1)=rownames(cor.filtered.2)=colnames(cor.filtered.1)=colnames(cor.filtered.2)<-rownames(exprs.1)
		rm(cor.filtered)  
  		degree.1 = apply(abs(cor.filtered.1)>=rth.1 & abs(cor.filtered.2)<rth.2,1,sum)
  		degree.2 = apply(abs(cor.filtered.1)<rth.1 & abs(cor.filtered.2)>=rth.2,1,sum)
  		ASC = as.vector(matrix(NA,m,1))
  		#rownames(ASC) = genes
  		ASC = (degree.1+degree.2)/2
  	}
  return(ASC)
  }

##############################################################################################################
## LRC
##############################################################################################################
"LRC" <-
function(exprs.1,exprs.2,method=Linkfilter.methods,cutoff) {
 	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
	genes = rownames(exprs.1)
  	if (nrow(exprs.1)!=nrow(exprs.2)) stop("invalid row lengths") else {
  		cor.filtered = switch(method,
		rth = {
  			genes <- rownames(exprs.1)
  			exprs.1 <- as.matrix(exprs.1)
  			exprs.2 <- as.matrix(exprs.2)
  			cor.pairwise.1 <- cor(t(exprs.1),method="pearson",use="pairwise.complete.obs")
  			cor.pairwise.2 <- cor(t(exprs.2),method="pearson",use="pairwise.complete.obs")
  			diag(cor.pairwise.1) <- 0
  			diag(cor.pairwise.2) <- 0
    		cor.filtered.1 <- cor.pairwise.1
  			cor.filtered.2 <- cor.pairwise.2
  			cor.filtered.1[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
  			cor.filtered.2[abs(cor.pairwise.1)<=cutoff & abs(cor.pairwise.2)<=cutoff] <- 0
   			list(rth=list(rth.1=cutoff, rth.2=cutoff),cor.filtered=list(cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2))
		
		},
	 	qth =  	qLinkfilter(exprs.1,exprs.2,cutoff),
		percent = percentLinkfilter(exprs.1,exprs.2,cutoff) )

		rth = cor.filtered$rth
		rth.1 = rth$rth.1
		rth.2 = rth$rth.2
		cor.filtered.1 = cor.filtered$cor.filtered$cor.filtered.1
		cor.filtered.2 = cor.filtered$cor.filtered$cor.filtered.2
		rownames(cor.filtered.1)=rownames(cor.filtered.2)=colnames(cor.filtered.1)=colnames(cor.filtered.2)<-rownames(exprs.1)
		rm(cor.filtered)  
  		degree.1 = apply(abs(cor.filtered.1)>=rth.1,1,sum)
  		degree.2 = apply(abs(cor.filtered.2)>=rth.2,1,sum)
  		degree.1[degree.1==0] = 1
  		degree.2[degree.2==0] = 1
  		LRC = as.vector(matrix(NA,m,1))
  		LRC = abs(log10(degree.2/degree.1))
  	}
return(LRC)
}
  	
##############################################################################################################
## WGCNA
##############################################################################################################
"WGCNA" <-
function(exprs.1,exprs.2,power=12,variant='WGCNA') {  	
 	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
	genes = rownames(exprs.1)
  	if (nrow(exprs.1)!=nrow(exprs.2)) stop("invalid row lengths") else {
  		genes <- rownames(exprs.1)
  		exprs.1 <- as.matrix(exprs.1)
  		exprs.2 <- as.matrix(exprs.2)
  		cor.pairwise.1 <- cor(t(exprs.1),method="pearson",use="pairwise.complete.obs")
  		cor.pairwise.2 <- cor(t(exprs.2),method="pearson",use="pairwise.complete.obs")
  		diag(cor.pairwise.1) <- 0
  		diag(cor.pairwise.2) <- 0
  		cor.pairwise.1 <- ((1+cor.pairwise.1)/2)^power
  		cor.pairwise.2 <- ((1+cor.pairwise.2)/2)^power
  		WGCNA = switch(variant,
  			DCp={ 	 				
  				squares = (cor.pairwise.1-cor.pairwise.2)^2
  				ss = apply(squares,1,sum)
   				sqrt(ss/(m-1))
   			},
  			WGCNA={
  				connectivity.1 = apply(cor.pairwise.1,1,sum)
   				connectivity.2 = apply(cor.pairwise.2,1,sum)
  				max.1 = max(connectivity.1)
  				max.2 = max(connectivity.2)
  				abs(connectivity.1/max.1 - connectivity.2/max.2)	
   				}
   				)
 }
  	return(WGCNA) 	  
}
