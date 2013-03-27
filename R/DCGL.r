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
"qLinkfilter" <-function(exprs.1,exprs.2,cutoff=0.25,r.method=c('pearson','spearman')[1],q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr")[1]) {
        # m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
        degree.1 <- ncol(exprs.1)-2
        degree.2 <- ncol(exprs.2)-2

        genes <- rownames(exprs.1)
        exprs.1 <- as.matrix(exprs.1)
        exprs.2 <- as.matrix(exprs.2)
        cor.1 <- cor(t(exprs.1),method=r.method,use="pairwise.complete.obs")
        cor.2 <- cor(t(exprs.2),method=r.method,use="pairwise.complete.obs")
        cor.1 <- cor.1[lower.tri(cor.1,diag=F)]
        cor.2 <- cor.2[lower.tri(cor.2,diag=F)]

        rm(exprs.1); rm(exprs.2)

	
        t.1 <- cor.1*sqrt(degree.1)/sqrt(1-cor.1*cor.1)
        t.2 <- cor.2*sqrt(degree.2)/sqrt(1-cor.2*cor.2)

        p0.1 <- 2*pt(-abs(t.1), degree.1, lower.tail = TRUE, log.p = FALSE)
        p0.2 <- 2*pt(-abs(t.2), degree.2, lower.tail = TRUE, log.p = FALSE)
        #diag(p0.1) <- NA
        #diag(p0.2) <- NA
        rm(t.1); rm(t.2)

        q.1<- p.adjust(p0.1, method = q.method)
        q.2<- p.adjust(p0.2, method = q.method)
        #dim(q.1)<- dim(p0.1); 
        #dim(q.2)<- dim(p0.2)
        #diag(q.1) <- 1
        #diag(q.2) <- 1
        rm(p0.1); rm(p0.2)

        rth.1 <- abs(cor.1[which.min(abs(q.1-cutoff))])
        rth.2 <- abs(cor.2[which.min(abs(q.2-cutoff))])
        cor.1[q.1>=cutoff & q.2>=cutoff] <- cor.2[q.1>=cutoff & q.2>=cutoff] <- 0
        
        #cor.1 <- cor.2 <- diag(rep(0,length(genes)))
        #cor.1[lower.tri(cor.1,diag=F)] = cor.1; cor.1 = cor.1+t(cor.1)
        #cor.2[lower.tri(cor.2,diag=F)] = cor.2; cor.2 = cor.2+t(cor.2)
        #rownames(cor.1) <- rownames(cor.2) <- colnames(cor.1) <- colnames(cor.2) <- genes


	name.row <- matrix(rep(genes,length(genes)),length(genes),length(genes))
	name.col <- matrix(rep(genes,length(genes)),length(genes),length(genes),byrow=T)
	name.pairs <- matrix(paste(name.row,name.col,sep=','),length(genes),length(genes))
	rm(list=c('name.row','name.col'))
	name.pairs <- name.pairs[lower.tri(name.pairs,diag=F)]
	names(cor.1) <- names(cor.2) <- name.pairs


        cor.filtered <- list(rth.1 = rth.1, rth.2 = rth.2, cor.filtered.1 = cor.1, cor.filtered.2 = cor.2, genes=genes)
        return(cor.filtered)
}


"rLinkfilter" <- function(exprs.1,exprs.2,cutoff=0.8,r.method=c('pearson','spearman')[1]) {

 	genes <- rownames(exprs.1)
        exprs.1 <- as.matrix(exprs.1)
        exprs.2 <- as.matrix(exprs.2)
        cor.filtered.1 <- cor(t(exprs.1),method=r.method,use="pairwise.complete.obs")
        cor.filtered.2 <- cor(t(exprs.2),method=r.method,use="pairwise.complete.obs")
        cor.filtered.1 <- cor.filtered.1[lower.tri(cor.filtered.1,diag=F)]
        cor.filtered.2 <- cor.filtered.2[lower.tri(cor.filtered.2,diag=F)]
        rm(exprs.1); rm(exprs.2)

        cor.filtered.1[abs(cor.filtered.1)<=cutoff & abs(cor.filtered.2)<=cutoff] <- 0
        cor.filtered.2[abs(cor.filtered.1)<=cutoff & abs(cor.filtered.2)<=cutoff] <- 0

	name.row <- matrix(rep(genes,length(genes)),length(genes),length(genes))
	name.col <- matrix(rep(genes,length(genes)),length(genes),length(genes),byrow=T)
	name.pairs <- matrix(paste(name.row,name.col,sep=','),length(genes),length(genes))
	rm(list=c('name.row','name.col'))
	name.pairs <- name.pairs[lower.tri(name.pairs,diag=F)]
	names(cor.filtered.1) <- names(cor.filtered.2) <- name.pairs

        list(rth.1=cutoff, rth.2=cutoff,cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2, genes=genes)

}
  ###############################################################################################################
  ## Select part of the correlation pairs, the max(abs(cor.1),abs(cor.2))
  ## exprs.1 a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.  
  ## exprs.2 a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.  
  ## percent percent of links to be reserved. 
  # output: A list with two components of lists: one lists the rth (thresholds of correlation values) for both conditions, the other lists the two matrices of filtered pairwise correlation values.
  ###############################################################################################################
"percentLinkfilter" <-function(exprs.1,exprs.2,cutoff=0.25,r.method=c('pearson','spearman')[1]) {
  	# m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
  	n.1 <- ncol(exprs.1)
  	n.2 <- ncol(exprs.2)
  
  	#degree.1 <- n.1-2
  	#degree.2 <- n.2-2
  
  	genes <- rownames(exprs.1)
  	exprs.1 <- as.matrix(exprs.1)
  	exprs.2 <- as.matrix(exprs.2)
  	cor.filtered.1 <- cor(t(exprs.1),method=r.method,use="pairwise.complete.obs")
  	cor.filtered.2 <- cor(t(exprs.2),method=r.method,use="pairwise.complete.obs")
  	cor.filtered.1 <- cor.filtered.1[lower.tri(cor.filtered.1,diag=F)]
	cor.filtered.2 <- cor.filtered.2[lower.tri(cor.filtered.2,diag=F)]
	rm(exprs.1); rm(exprs.2)
	
	#diag(cor.filtered.1) <- 0
  	#diag(cor.filtered.2) <- 0
	cor.1.2.max <- pmax(abs(cor.filtered.1),abs(cor.filtered.2))
	#cor.1.2.max <- cor.1.2.max[lower.tri(cor.1.2.max,diag=F)]
  	#cor.1.2.vector <- cbind(cor.pairwise.1[lower.tri(cor.pairwise.1, diag = FALSE)],cor.pairwise.2[lower.tri(cor.pairwise.2, diag = FALSE)])
  	#cor.1.2.max <- apply(abs(cor.1.2.vector),1,max)  #max of the two cor
  	#cat('passed lower tri extraction \n')
	cor.1.2.max <- sort(cor.1.2.max,decreasing = TRUE);
	#cat('passed max corr operation \n')
  	# cor.1.2.max.sort <- as.matrix(cor.1.2.max.sort)
  	Rth <- cor.1.2.max[as.integer(length(cor.1.2.max)*cutoff)];
  	rm(cor.1.2.max)
  	#cor.filtered.1 <- cor.pairwise.1
  	#cor.filtered.2 <- cor.pairwise.2
	#rm(cor.pairwise.1); rm(cor.pairwise.2)
  	cor.filtered.1[abs(cor.filtered.1)<=Rth & abs(cor.filtered.2)<=Rth] <- 0
  	cor.filtered.2[abs(cor.filtered.1)<=Rth & abs(cor.filtered.2)<=Rth] <- 0
	#cat('passed 0 substitution \n')
  	
	#cor.2 <- cor.1 <- diag(rep(0,length(genes)))
	#cor.1[lower.tri(cor.1,diag=F)] = cor.filtered.1; cor.1 = cor.1+t(cor.1)
	#cor.2[lower.tri(cor.2,diag=F)] = cor.filtered.2; cor.2 = cor.2+t(cor.2)
	#rownames(cor.1) <- rownames(cor.2) <- colnames(cor.1) <- colnames(cor.2) <- genes
	#rm(cor.filtered.1); rm(cor.filtered.2)
	#cat('reform to two matrices\n')

	name.row <- matrix(rep(genes,length(genes)),length(genes),length(genes))
	name.col <- matrix(rep(genes,length(genes)),length(genes),length(genes),byrow=T)
	name.pairs <- matrix(paste(name.row,name.col,sep=','),length(genes),length(genes))
	name.pairs <- name.pairs[lower.tri(name.pairs,diag=F)]
	rm(list=c('name.row','name.col'))
	names(cor.filtered.1) <- names(cor.filtered.2) <- name.pairs

	cor.filtered <- list(rth.1=Rth,rth.2=Rth,cor.filtered.1 = cor.filtered.1, cor.filtered.2 = cor.filtered.2, genes=genes)
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
 }else {
	 		exprs.min <- apply(abs(exprs),1,min)
  			exprs.max <- apply(abs(exprs),1,max)
  			exprs.diff <- exprs.max/exprs.min
  			a <- exprs.max
	 		exprs.max <- exprs.diff
			exprs.diff <- a	
	 	}
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
	# For bin 1# to 2# to 3# ……, each time get the decisive point's x & y coordinates - the gradually elongating vectors exprs.x.sub & exprs.y.sub.
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
  		} else {
  				plot(log2(exprs.max),(exprs.diff),pch=16,col='darkgray',xlab='log2(corr ratio)',ylab='max(corr)',main='LFC')
  		}
	# important changement: steps.x substitutes for exprs.x.sub.
	mm<-glm(y~x,data=data.frame(y=exprs.y.sub,x=1/steps.x))
	exprs.diff.threshold<-predict(mm,newdata=data.frame(x=1/exprs.max))
	delink<- (exprs.diff>exprs.diff.threshold)
	
		if(sign=='same'){
			points((exprs.max[delink]),log2(exprs.diff[delink]),pch=16,col='lightblue')
			lines((exprs.max),log2(exprs.diff.threshold),col='red')
			points((steps.x),log2(exprs.y.sub),pch=16)
		} else {
			points(log2(exprs.max[delink]),(exprs.diff[delink]),pch=16,col='lightblue')
			exprs.max <- sort(exprs.max);
			exprs.diff.threshold <- sort(exprs.diff.threshold);
			lines(log2(exprs.max),(exprs.diff.threshold),col='red')
			points(log2(steps.x),(exprs.y.sub),pch=16)
		}
	dev.off()
	delink
  }

  #####################################################################################
  ## 'exprs.1' a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.  
  # 'exprs.2' a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.  
  # 'N' ramdom sampling times to form the NULL distribution. If N>0 ramdom sampling is used to form the NULL distribution. 
  #####################################################################################

"DCp" <- function(exprs.1,exprs.2,r.method=c('pearson','spearman')[1],link.method=c('qth','rth','percent')[1],cutoff=0.25,N=0,N.type=c('pooled','gene_by_gene')[1],q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr")[1]) {  
	cor.filtered.1<-cor.filtered.2<-NULL
	#Linkfilter.method=c('rth','qth','percent')
	if (nrow(exprs.1)!=nrow(exprs.2)) stop('the two expression matrices must have the same number of rows (genes).')
	if (length(rownames(exprs.1))==0 | length(rownames(exprs.2))==0) stop('the expression matrices must have row names specifying the gene names.')
	if ( min(ncol(exprs.1),ncol(exprs.2))<3 ){
		stop('each expression matrix must have at least three or more columns.')
	} else if (min(ncol(exprs.1),ncol(exprs.2))<5 ) {
		warning('the minimum number of columns is less than five and the result may not be reliable.')
	}

	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
  	if (m>5000) warning('the number of genes exceeds 5000 and the program may takes long time to run.')
	genes = rownames(exprs.1)
  	
	cor.filtered = switch(link.method,
		rth = rLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
	 	qth =  	qLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
		percent = percentLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff)
	)

	#attach(cor.filtered)
	rth.1 = cor.filtered$rth.1
	rth.2 = cor.filtered$rth.2
	cor.filtered.1 = cor.filtered$cor.filtered.1
	cor.filtered.2 = cor.filtered$cor.filtered.2

#	rth.1 = rth$rth.1
#	rth.2 = rth$rth.2
	
#	rownames(cor.filtered.1)=rownames(cor.filtered.2)=colnames(cor.filtered.1)=colnames(cor.filtered.2)<-genes
 
  #####################################################################################################################################################
  ## For one gene, there are n-1 correlation value pairs also q value pairs(From the two conditions). For a correlation pair for this gene,
  ## if one of the q values is less than the threshold, the pair will be retained. Then there are 'number.i.uniq' unique pairs retained that is there 
  ## are two vectors of correlation values. 
  ## Then a length normalized Euclidean Distance for the two vectors will be calculated (LNED). 
  ##################################################################################################################################################### 
	calc.dC <- function(cor.filtered.1,cor.filtered.2,genes) {
		nzero.vec <- (cor.filtered.1 != 0)|(cor.filtered.2 != 0 )
		nzero.sm <- diag(rep(0,length(genes)))
		nzero.sm[lower.tri(nzero.sm,diag=F)] <- nzero.vec; nzero.sm = nzero.sm+t(nzero.sm)
		number_uniq <- apply(nzero.sm,1,sum)
 
 		squares = (cor.filtered.1-cor.filtered.2)^2
  #	number_uniq = apply(cor.filtered.1!=0 | cor.filtered.2!=0,1,sum)
#  	ss = apply(squares,1,sum)
		squares.sm <- diag(rep(0,length(genes)))
		squares.sm[lower.tri(squares.sm,diag=F)] <- squares; squares.sm = squares.sm+t(squares.sm)
		ss = apply(squares.sm,1,sum)  
		LNED.result = as.vector(matrix(NA,length(genes),1))
  		LNED.result[number_uniq!=0] = sqrt(ss[number_uniq!=0])/sqrt(number_uniq[number_uniq!=0])
  		names(LNED.result) <- genes
		list(dC=LNED.result,length=number_uniq)
	}

	dC.length = calc.dC(cor.filtered.1,cor.filtered.2,genes)	
	dC = dC.length$dC
	number_uniq = dC.length$length
#	rm(list=c('nzero.vec','nzero.sm','squares','squares.sm','ss'))
	
	#detach(cor.filtered)

 
  
 ########################################################################################################################
 ## Disturb the sample labels for the two conditions and re-assign the samples to two datasets,then calculate the 'dC0' for 
 ## N times and then pool all the dC0 together to construct a 'NULL' distribution.
 #########################################################################################################################
	if(N>0){
		dC0 <- matrix(nrow=length(genes),ncol=N)
		rownames(dC0) <- genes
		exprs <- cbind(exprs.1,exprs.2)
		n.1 = ncol(exprs.1)
		n.2 = ncol(exprs.2)
		cat.j = 0
		for(j in 1:N) {
			if ( (j*100/N)%/%10>cat.j) {
				cat.j = cat.j+1
				cat(cat.j*10,'%','\n')
			}
			seq <- sample(n.1+n.2)
			exprs.1 <- exprs[,seq[1:n.1]];
			exprs.2 <- exprs[,seq[(n.1+1):(n.1+n.2)]];
			rownames(exprs.1) <- rownames(exprs.2) <- genes

			cor.filtered = switch(link.method,
                		rth = rLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
                		qth =   qLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
                		percent = percentLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff)
        		)
#			attach(cor.filtered)
#  			corfiltered_1 = cor_filtered$cor.filtered$cor.filtered.1
#  			cor_filtered_2 = cor_filtered$cor.filtered$cor.filtered.2
#  			rm(cor_filtered)

#  			squares = (cor_filtered_1-cor_filtered_2)^2
#  			number_uniq_0 = apply(cor_filtered_1!=0 | cor_filtered_2!=0,1,sum)
#  			ss = apply(squares,1,sum)
#  			LNED.result.0 = numeric(number_uniq_0)
#  			LNED.result.0[number_uniq_0!=0] = sqrt(ss[number_uniq_0!=0])/sqrt(number_uniq[number_uniq_0!=0])
			dC0.j <- calc.dC(cor.filtered$cor.filtered.1,cor.filtered$cor.filtered.2,cor.filtered$genes)
  			dC0[,j] = dC0.j$dC
			
#			detach(cor.filtered)
   		}

		p.value = switch(N.type,
			gene_by_gene = apply(cbind(dC0,dC),1,function(x) sum(x[1:(length(x)-1)]>x[length(x)],na.rm=T)/length(!is.na(x[1:(length(x)-1)]))),
			pooled = sapply(dC,function(x) sum(as.vector(dC0)>x,na.rm=T)/length(!is.na(as.vector(dC0))))
		)
#		p.value <- vector()
#		for(k in 1:m){
#			p <- sum(LNED.result[k]<dC0)/(N*m);
#			p.value <- c(p.value,p);
#		}
  		q.value <- p.adjust(p.value,method=q.method)
  		Result <- data.frame(dC=dC,links=number_uniq,p.value=p.value,q.value=q.value);
  		row.names(Result) <- genes;
  		#colnames(Result) <- paste(c("dC","length","p.value","q.value"));
  		# dec.idx <-sort(as.numeric(Result[,4]),method = "quick", index.return=TRUE, decreasing=FALSE)$ix
  		# Result <- Result[dec.idx,]
	}
	else { 

  
  		Result<- data.frame(dC=dC,links=number_uniq,p.value=rep(NA,length(dC)),q.value=rep(NA,length(dC)));
  		row.names(Result) <- genes;
  		#colnames(Result) <- c("dC","length");
  		# dec.idx <-sort(as.numeric(Result[,1]),method = "quick", index.return=TRUE,decreasing=TRUE)$ix
  		# Result <- as.matrix(Result[dec.idx,])
  	} 
  	return(Result)
	
}
  #####################################################################################
  ## 'exprs.1' a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.  
  # 'exprs.2' a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.  
  # 'link.method' a character string for link method.
  #####################################################################################

"DCe" <-
function(exprs.1,exprs.2,link.method=c('qth','rth','percent')[1],cutoff=0.25,r.method=c('pearson','spearman')[1],q.method=c("BH","holm", "hochberg", "hommel", "bonferroni", "BY","fdr")[1],nbins=20,p=0.1,figname = c('LFC.s.jpeg','LFC.d.jpeg')) {
	cor.filtered.1<-cor.filtered.2<-rth.1<-rth.2<-NULL
	if (nrow(exprs.1)!=nrow(exprs.2)) stop('the two expression matrices must have the same number of rows (genes).')
	if (length(rownames(exprs.1))==0 | length(rownames(exprs.2))==0) stop('the expression matrices must have row names specifying the gene names.')
	if ( min(ncol(exprs.1),ncol(exprs.2))<3 ){
		stop('each expression matrix must have at least three or more columns.')
	} else if (min(ncol(exprs.1),ncol(exprs.2))<5 ) {
		warning('the minimum number of columns is less than five and the result may not be reliable.')
	}

 	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
	genes = rownames(exprs.1)

	cor.filtered = switch(link.method,
		rth = rLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
	 	qth =  	qLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff),
		percent = percentLinkfilter(exprs.1,exprs.2,r.method=r.method,cutoff=cutoff)
	)

	#attach(cor.filtered)
	rth.1 = cor.filtered$rth.1
	rth.2 = cor.filtered$rth.2
	cor.filtered.1 = cor.filtered$cor.filtered.1
	cor.filtered.2 = cor.filtered$cor.filtered.2

	
  	#############################################################
  	## decide three sets of correlation pairs and organize them into two-columned matrices.
  	#############################################################  	
  	
  	idx.same = (cor.filtered.1*cor.filtered.2)>0
  	idx.diff = (cor.filtered.1*cor.filtered.2)<0
  	idx.switched = (cor.filtered.1*cor.filtered.2<0)& ( abs(cor.filtered.1)>=rth.1 & abs(cor.filtered.2)>=rth.2 )
  	
  	cor.same = cbind(cor.filtered.1[idx.same],cor.filtered.2[idx.same])
  	rownames(cor.same) <- names(cor.filtered.1)[idx.same]
	cor.switched = cbind(cor.filtered.1[idx.switched],cor.filtered.2[idx.switched])
	rownames(cor.switched) <- names(cor.filtered.1)[idx.switched]  	
  	cor.diff = cbind(cor.filtered.1[idx.diff & (!idx.switched)],cor.filtered.2[idx.diff & (!idx.switched)])
	rownames(cor.diff) <- names(cor.filtered.1)[idx.diff & (!idx.switched)]
 	#detach(cor.filtered)

	n.switchedDCL = nrow(cor.switched)
# use strsplit to get two-column edge specification. 	
# 	name.rows = matrix(rep(genes,m),m,m,byrow=T)
#  	name.columns = matrix(rep(genes,m),m,m)
  	
	if ( is.null( rownames(cor.same) ) ) {name.same = NULL}
	if ( !is.null( rownames(cor.same) ) ) {
		name.same = strsplit(rownames(cor.same),',')
		name.same = matrix(unlist(name.same),length(name.same),2,byrow=T)
	}
	if ( is.null( rownames(cor.switched) ) ) {name.switched = NULL}
	if ( !is.null( rownames(cor.switched) ) ) {
		name.switched = strsplit(rownames(cor.switched),',')
		name.switched = matrix(unlist(name.switched),length(name.switched),2,byrow=T)
	}
	if ( is.null( rownames(cor.diff) ) ) {name.diff = NULL}
	if ( !is.null( rownames(cor.diff) ) ) {
		name.diff = strsplit(rownames(cor.diff),',')
		name.diff = matrix(unlist(name.diff),length(name.diff),2,byrow=T)
	}
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
	p.value <- pbinom(degree.bind[,'DC.links']-1, degree.bind[,'All.links'], prob, lower.tail = F, log.p = FALSE);
 	q.value <- p.adjust(p.value,method=q.method);
 
 	degree.bind <- cbind(degree.bind,p.value,q.value)
 	colnames(degree.bind) <- c("All.links","DC.links","DCL_same","DCL_diff","DCL_switch","p","q")

 	middle <-sort(as.numeric(degree.bind[,'q']),method = "quick", decreasing=FALSE,index.return=TRUE)$ix 
 	DCGs <- degree.bind[middle,]
 
 #########################################################
 
	DCLs <- rbind(data.frame(DCL.same,type='same signed'),data.frame(DCL.diff,type='diff signed'))
	if (n.switchedDCL>0) DCLs <- rbind(DCLs,data.frame(DCL.switched,type='switched opposites'))
	DCLs <- data.frame(DCLs,cor.diff=FALSE)
	DCLs[,'cor.diff'] <- abs(DCLs[,'cor.1']-DCLs[,'cor.2'])
 	Result <- list(DCGs=DCGs,DCLs=DCLs)
	return(Result)

}


  ##########################################################DCsum##########################################
  ## DCsum 
  # Summary from DCp and DCe 
  # DCpcutoff is a numeric of 'q.value' for DCGs threshold from DCp.
  # DCecutoff is a numeric of 'q' for DCGs threshold from DCe.
  ##########################################################################################################
"DCsum" <- function(DCp.res,DCe.res,DCpcutoff=0.25,DCecutoff=0.25) {
	
	if (sum(is.na(DCp.res$q.value))==nrow(DCp.res)) {
		DCp.res <- DCp.res[order(-DCp.res[,'dC']),]
		DCp.DCG <- DCp.res[1:ceiling(nrow(DCp.res)*DCpcutoff),]
	} else {
		DCp.DCG <- DCp.res[DCp.res$q.value<DCpcutoff,]
	}
        
	DCp.DCG <- data.frame(DCG=rownames(DCp.DCG),DCp.DCG)
	
	DCe.DCG <- DCe.res$DCGs[DCe.res$DCGs[,'q']<DCecutoff,]
	DCe.DCG <- data.frame(DCG=rownames(DCe.DCG),DCe.DCG)
	

	DCG <- merge(DCp.DCG,DCe.DCG,by.x='DCG',by.y='DCG')
	if (nrow(DCG)<1) stop('there is no intersection between DCp DCGs and DCe DCGs')

		
	colnames(DCG) <- c('DCG','dC','All.links.DCp','DCp.p','DCp.q','All.links.DCe','DC.links','DCL.same','DCL.diff','DCL.switch','DCe.p','DCe.q')
	
	DCL0 <- DCe.res$DCLs
	DCL.rddt1 <- merge(DCL0,data.frame(DCG=DCG$DCG),by.y='DCG',by.x='Gene.1'); DCL.rddt1 <- data.frame(DCL.rddt1,DCG=DCL.rddt1[,'Gene.1'])
	DCL.rddt2 <- merge(DCL0,data.frame(DCG=DCG$DCG),by.y='DCG',by.x='Gene.2'); DCL.rddt2 <- data.frame(DCL.rddt2,DCG=DCL.rddt2[,'Gene.2'])
	DCL.rddt <- rbind(DCL.rddt1,DCL.rddt2)
	DCL.rddt.names <- paste(DCL.rddt$Gene.1,DCL.rddt$Gene.2,sep='; ')
	
	link.count <- as.data.frame(table(DCL.rddt.names))
	link.idx <- link.count[,2]
	names(link.idx) <- link.count[,1]
	DCL <- data.frame(unique(DCL.rddt[,1:(ncol(DCL.rddt)-1)]),DCG=NA)
	rownames(DCL.rddt) <- paste('pair',1:nrow(DCL.rddt))
	rownames(DCL.rddt)[DCL.rddt.names %in% names(link.idx[link.idx==1])] <- DCL.rddt.names[DCL.rddt.names %in% names(link.idx[link.idx==1])]  
	rownames(DCL) <- paste(DCL$Gene.1,DCL$Gene.2,sep='; ')
	DCL[names(link.idx)[link.idx==1],'DCG'] <- as.character(DCL.rddt[names(link.idx[link.idx==1]),'DCG'])
	if ( dim(DCL[names(link.idx)[link.idx==2],])[1]>0 ) {
		DCL[names(link.idx)[link.idx==2],'DCG'] <- names(link.idx)[link.idx==2]
	}
	
	#list(DCGs=DCG,DCLs=DCL,DCL.rddt=DCL.rddt,link.idx=link.idx)
	list(DCGs=DCG,DCLs=DCL)
}
##############################################################################################################
## ASC
##############################################################################################################
"ASC" <-
function(exprs.1,exprs.2,link.method=c('qth','rth','percent')[1],cutoff) {
 	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
	genes = rownames(exprs.1)
  	if (nrow(exprs.1)!=nrow(exprs.2)) stop("invalid row lengths") else {
  		cor.filtered = switch(link.method,
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


		#rth = cor.filtered$rth
		rth.1 = cor.filtered$rth.1
		rth.2 = cor.filtered$rth.2
		cor.filtered.1 = cor.filtered$cor.filtered.1
		cor.filtered.2 = cor.filtered$cor.filtered.2
		tri.cor.1 <- diag(rep(0,length(genes)))
		tri.cor.1[lower.tri(tri.cor.1,diag=F)] <- cor.filtered.1; who.cor.1 = tri.cor.1+t(tri.cor.1)
		tri.cor.2 <- diag(rep(0,length(genes)))
		tri.cor.2[lower.tri(tri.cor.2,diag=F)] <- cor.filtered.2; who.cor.2 = tri.cor.2+t(tri.cor.2)
		rownames(who.cor.1)=rownames(who.cor.2)=colnames(who.cor.1)=colnames(who.cor.2)<-rownames(exprs.1)
		rm(cor.filtered)  
  		degree.1 = apply(abs(who.cor.1)>=rth.1 & abs(who.cor.2)<rth.2,1,sum)
  		degree.2 = apply(abs(who.cor.1)<rth.1 & abs(who.cor.2)>=rth.2,1,sum)
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
function(exprs.1,exprs.2,link.method=c('qth','rth','percent')[1],cutoff) {
 	m <- nrow(exprs.1) # exprs.1, exprs.2 is the expression data for different conditions.
	genes = rownames(exprs.1)
  	if (nrow(exprs.1)!=nrow(exprs.2)) stop("invalid row lengths") else {
  		cor.filtered = switch(link.method,
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

		#rth = cor.filtered$rth
		rth.1 = cor.filtered$rth.1
		rth.2 = cor.filtered$rth.2
		cor.filtered.1 = cor.filtered$cor.filtered.1
		cor.filtered.2 = cor.filtered$cor.filtered.2
		tri.cor.1 <- diag(rep(0,length(genes)))
		tri.cor.1[lower.tri(tri.cor.1,diag=F)] <- cor.filtered.1; who.cor.1 = tri.cor.1+t(tri.cor.1)
		tri.cor.2 <- diag(rep(0,length(genes)))
		tri.cor.2[lower.tri(tri.cor.2,diag=F)] <- cor.filtered.2; who.cor.2 = tri.cor.2+t(tri.cor.2)
		rownames(who.cor.1)=rownames(who.cor.2)=colnames(who.cor.1)=colnames(who.cor.2)<-rownames(exprs.1)
		rm(cor.filtered)  
  		degree.1 = apply(abs(who.cor.1)>=rth.1,1,sum)
  		degree.2 = apply(abs(who.cor.2)>=rth.2,1,sum)
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

  ########################################################DRsort##########################################################
  ## DRsort To identify DRGs (Differentially-regulated Genes) and DRLs (Differentially-regulated Links)
  # DCGs: a data frame or matrix for DCGs list.
  # DCLs: a data frame or matrix for DCLs list.
  # tf2target: a data frame or matrix for regulator-to-target interaction pairs.
  ######################################################################################################################
"DRsort" <- function(DCGs, DCLs, tf2target, exprs) {

	######有表达谱的TF#############
	TFinexpression<-merge(tf2target,exprs,by.x='TF',by.y='row.names')
	TFinexpression<-unique(as.character(TFinexpression[,'TF']))

	tf2target <- matrix(as.character(unlist(tf2target)),ncol = ncol(tf2target))
	tf2target <- data.frame(tf2target)
	colnames(tf2target) <- c('TF','target')

	#######################################     DCG: add one field 'DCGisTF' to indicate whether the DCG is a TF or not        ##################################
	#cat('DCG: Add one field \'DCGisTF\' to indicate whether the DCG is a TF or not\n')
	TFs <- unique(tf2target$TF)
	DCGs <- data.frame(DCGisTF=FALSE,DCGs)
	TF.idx <- as.character(DCGs$DCG) %in% TFs
	DCGs$DCGisTF[TF.idx] <- TRUE
	DCGs <- DCGs[,c('DCG','DCGisTF','dC','DCp.p','All.links.DCe','DC.links','DCL.same','DCL.diff','DCL.switch')]
	
	#DCGisTF<-merge(tf2target,DCG,by.x="TF",by.y="DCG",all.y=T)
	#colnames(DCGisTF)[1:2] <- c('DCG','target')	
	#DCGisTF <- DCGisTF[,c('DCG','target','geneid','keggid','keggname','dC','DCp.p','All.links.DCe','DC.links','DCL.same','DCL.diff','DCL.switch')]
	#DCGisTF[!is.na(DCTisTF$target),'target'] <- 'True'
	#DCGisTF[is.na(DCTisTF$target),'target'] <- 'False'
	#colnames(DCGisTF)[2] <- 'DCGisTF'
	#DCGisTF <- unique(DCGisTF)

	#DCG <- DCGisTF
	#######################################   DCG: add one field 'Upstream_TFofDCG' to indicate the TF(s) of the DCG    #############################
	#cat('DCG: Add one field \'Upstream_TFofDCG\' to indicate the TF(s) of the DCG\n')
	DCGistarget<-merge(tf2target,DCGs,by.x="target",by.y="DCG",all.y=T)
	colnames(DCGistarget)[1:2] <- c('DCG','Upstream_TFofDCG')
	DCGistarget <- DCGistarget[,c('DCG','Upstream_TFofDCG','DCGisTF','dC','DCp.p','All.links.DCe','DC.links','DCL.same','DCL.diff','DCL.switch')]
	DCG2TF <- data.frame(DCG=DCGistarget$DCG,TF=DCGistarget$Upstream_TFofDCG)    ######### this object is outputted for easing TED calculation.
	
	
	
	DCGs.TF <- as.character(unlist(by(DCG2TF$TF,factor(DCG2TF$DCG,ordered=T),paste,collapse=';')))	
	DCG2TF.byDCG <- data.frame(DCG=levels(factor(DCG2TF$DCG,ordered=T)),TF=DCGs.TF)
	DCGs <- merge(DCG2TF.byDCG,DCGs,by.x='DCG',by.y='DCG',all.y=T)
	colnames(DCGs) <- c('DCG','Upstream_TFofDCG','DCGisTF','dC','DCp.p','All.links.DCe','DC.links','DCL.same','DCL.diff','DCL.switch')	

	#######################################  DCL: Add one field 'TF' to indicate the internal TF of the DCL pair   ############################
	#colnames(DCL) c('Gene.1','Gene.2',,'cor.1','cor.2','type','DCG','cor.diff')
	#cat('DCL: Add one field \'TF\' to indicate the internal TF of the DCL pair\n')
	tf.target <- paste(tf2target[,'TF'],tf2target[,'target'],sep='; ')
	DCL.pair1 <- paste(DCLs$Gene.1,DCLs$Gene.2, sep='; ')
	TF1.idx <- DCL.pair1 %in% tf.target
	
	DCL.pair2 <- paste(DCLs$Gene.2,DCLs$Gene.1, sep='; ')
	TF2.idx <- DCL.pair2 %in% tf.target
	
	DCLs <- data.frame(TF=NA,DCLs)
	DCLs[TF1.idx,'TF'] <- as.character(DCLs[TF1.idx,'Gene.1'])
	DCLs[TF2.idx,'TF'] <- as.character(DCLs[TF2.idx,'Gene.2'])
	DCLs[TF1.idx & TF2.idx, 'TF'] <- paste(as.character(DCLs[TF1.idx & TF2.idx,'Gene.1']),as.character(DCLs[TF1.idx & TF2.idx,'Gene.2']),sep='; ')
	
	colnames(DCLs) <- c('internal.TF','Gene.1','Gene.2','cor.1','cor.2','type','cor.diff','DCG')

	#DCL1 <- DCL
	#######################################  DCL: Add one field 'common.TF' to indicate the common TF(s) of the DCL pair        ##############
	#colnames(DCL) c('internal.TF','Gene.1','Gene.2','cor.1','cor.2','type','DCG','cor.diff')
	#cat('DCL: Add one field \'common.TF\' to indicate the common TF(s) of the DCL pair \n')
	tfbridgedDCL<-merge(DCLs,tf2target,by.x="Gene.1",by.y="target") ####nrow of tfbridgedDCL is less than nrow of DCL###
	colnames(tfbridgedDCL)<- c('Gene.1','internal.TF','Gene.2','cor.1','cor.2','type','cor.diff','DCG','TF.of.g1')   ## add the 'TF' that regulates 'Gene.1'
	
	tfbridgedDCL <- merge(tf2target,tfbridgedDCL,by.x=c('TF','target'),by.y=c('TF.of.g1','Gene.2'))
	colnames(tfbridgedDCL) <- c('common.TF','Gene.2','Gene.1','internal.TF', 'cor.1','cor.2','type','cor.diff','GeneisDCG')  ## extract the rows in which 'TF' that regulates Gene.1 and Gene.2 both.
	tfbridgedDCL <- unique(tfbridgedDCL)
	
	tfbridgedDCL<-data.frame(common.TFisDCG=FALSE,tfbridgedDCL)
	common.TFDCG.idx <- as.character(tfbridgedDCL$common.TF) %in% as.character(DCGs[DCGs[,'DCGisTF'],'DCG'])
	tfbridgedDCL$common.TFisDCG[common.TFDCG.idx]<-TRUE
#	tfbridgedDCL<-tfbridgedDCL[,c('common.TF','common.TFisDCG','Gene.2','Gene.1','internal.TF', 'geneid.1','keggid.1','keggname.1','geneid.2','keggid.2','keggname.2','cor.1','cor.2','type','GeneisDCG')]
	
	tfbridgedDCL<-data.frame(common.TFinexprs=FALSE,tfbridgedDCL)
	common.TF.exprsid <- as.character(tfbridgedDCL$common.TF) %in% TFinexpression
	tfbridgedDCL$common.TFinexprs[common.TF.exprsid]<-TRUE
	tfbridgedDCL<-tfbridgedDCL[,c('common.TF','Gene.1','Gene.2','internal.TF','common.TFisDCG','common.TFinexprs', 'cor.1','cor.2','type','cor.diff','GeneisDCG')]

	commonTF <- as.character(unlist( by(tfbridgedDCL$common.TF,factor(as.character(paste(tfbridgedDCL$Gene.1,tfbridgedDCL$Gene.2,sep='; ')),ordered=T),paste,collapse='; ')  ))        
	DCLsst.pairID <- levels(factor(as.character(paste(tfbridgedDCL$Gene.1,tfbridgedDCL$Gene.2,sep='; ')),ordered=T))
	commonTF <- data.frame(pairID=DCLsst.pairID,common.TF=commonTF)
	DCL.pairID <- paste(DCLs$Gene.1,DCLs$Gene.2,sep='; ')
	DCLs <- data.frame(pairID=DCL.pairID,DCLs)
	
	DCLs <- merge(commonTF,DCLs,by.x='pairID',by.y='pairID',all.y=T)
	
	DRGs<-DCGs[DCGs[,'DCGisTF'] | !is.na(DCGs[,'Upstream_TFofDCG']),]
	DRLs<-DCLs[!is.na(DCLs[,'common.TF']) | !is.na(DCLs[,'internal.TF']),]	
		
	
	list(DCGs=DCGs, DCLs=DCLs, DRGs=DRGs, DRLs=DRLs, DCG2TF=DCG2TF, TF_bridged_DCL=tfbridgedDCL)
}


  ###############################################visualization of TF2target_DCL#########################################################
  ## visualization of DCL is TF2target relation
  # TF2target_DCL: Tf2target identified from DCLs
  #######################################################################################################################################
 visTF2target_DCL<-function(TF2target_DCL.int,vsize=vsize,asize=asize,lcex=lcex,ewidth=ewidth,figname){
	if (!require(igraph)) stop('package igraph is required!')
	if ( is.null(TF2target_DCL.int)) stop ('There are no interesting genes in TF2target_DCL.\n')
	if ( nrow(TF2target_DCL.int)==0 ) stop('There are no TF2target_DCL.\n')
	# colnames of DCListf2target c('TF','keggname.TF','Gene','keggname.Gene','DCG')
	#if (nrow(DCListf2target)==0) stop ('there are no tf2target interaction pairs.\n')
	if (nrow(TF2target_DCL.int)>0) {
		DCL <- TF2target_DCL.int
		DCG <- unique(unlist(strsplit(DCL$DCG,'; '))) #provided that nested DCGs are separated as "A; B"
		TF <- unique(DCL$TF)
		node <- unique(c(as.character(DCL$TF),as.character(DCL$Gene)))
		node.classes <- data.frame(TF= node %in% TF, DCG= node %in% DCG)
		rownames(node.classes) <- as.character(node)
		
		relation<-data.frame(DCL$TF,DCL$Gene)
		if (nrow(relation)>1000) {
			warning ("the edge number of TF2target_DCL(>1000) is too large to display clearly and the program maybe need long time to run.\n")
		}
	
		g <- graph.data.frame(relation)
		
		#browser()
		if (length(V(g))!=nrow(node.classes)) {
			stop('oops - why rows of node classes not equal to nodes of graph?\n')
		} else {
			node.classes = node.classes[V(g)$name,]
			V(g)$shape <- 'circle'; V(g)$shape[node.classes$TF] <- 'square'
			V(g)$color <- 'skyblue'; V(g)$color[node.classes$DCG] <- 'pink'
			E(g)$color <- "black"
			E(g)$width <- ewidth
			E(g)$arrow.size <- asize
			V(g)$size <- vsize
			V(g)$label.color <- "black";
			V(g)$label.cex <- lcex; 
			V(g)$label <- V(g)$name
			pdf(figname[1])
			plot(g,layout=layout.fruchterman.reingold)
			dev.off()
			cat("The graph of", figname[1], "has been completed and saved in your working directory.\n")
		}
	}
}


  ###############################################visualization of TF_bridged_DCL#########################################################
  ## visualization of TF_bridged_DCL
  # TF_bridged_DCL: both genes in a DCL regulated by a common regulator.
  #####################################################################################################################################
 visTF_bridged_DCL <- function(TF2target_DCL,TF_bridged_DCL,vsize=vsize,asize=asize,lcex=lcex,ewidth=ewidth,figname){
	if (!require(igraph)) stop('package igraph is required!\n')
	if ( is.null(TF_bridged_DCL) ) stop('There are no interesting genes in TF_bridged_DCL.\n')
	if (nrow(TF_bridged_DCL)==0) stop('There are no TF_bridged_DCL.\n')
	if (length(TF_bridged_DCL[is.na(TF_bridged_DCL[,'internal.TF']),'Gene.1']) >0 ){
		relation.DCL <- data.frame(node1=TF_bridged_DCL[is.na(TF_bridged_DCL[,'internal.TF']),'Gene.1'],node2=TF_bridged_DCL[is.na(TF_bridged_DCL[,'internal.TF']),'Gene.2'],regulation=FALSE,DCL=TRUE)
	} else {relation.DCL=NULL}
	if ( length( TF_bridged_DCL[!is.na(TF_bridged_DCL[,'internal.TF']) & TF_bridged_DCL[,'internal.TF']==TF_bridged_DCL[,'Gene.1'],'Gene.1'] ) >0 ) {
		relation.DCL.oneReg.1 <- data.frame(node1=TF_bridged_DCL[!is.na(TF_bridged_DCL[,'internal.TF']) & TF_bridged_DCL[,'internal.TF']==TF_bridged_DCL[,'Gene.1'],'Gene.1'], node2=TF_bridged_DCL[!is.na(TF_bridged_DCL[,'internal.TF']) & TF_bridged_DCL[,'internal.TF']==TF_bridged_DCL[,'Gene.1'],'Gene.2'], regulation=TRUE, DCL=TRUE)
	} else { relation.DCL.oneReg.1=NULL}
	if ( length( TF_bridged_DCL[!is.na(TF_bridged_DCL[,'internal.TF']) & TF_bridged_DCL[,'internal.TF']==TF_bridged_DCL[,'Gene.2'],'Gene.2'] ) >0 ){
		relation.DCL.oneReg.2 <- data.frame(node1=TF_bridged_DCL[!is.na(TF_bridged_DCL[,'internal.TF']) & TF_bridged_DCL[,'internal.TF']==TF_bridged_DCL[,'Gene.2'],'Gene.2'], node2=TF_bridged_DCL[!is.na(TF_bridged_DCL[,'internal.TF']) & TF_bridged_DCL[,'internal.TF']==TF_bridged_DCL[,'Gene.2'],'Gene.1'], regulation=TRUE, DCL=TRUE)
	} else {relation.DCL.oneReg.2=NULL}
	
	if (length(TF_bridged_DCL[grepl(";", TF_bridged_DCL[,'internal.TF']),'Gene.1'])!=0){
		relation.DCL.twoReg.half <- data.frame(node1=TF_bridged_DCL[grepl(";", TF_bridged_DCL[,'internal.TF']),'Gene.1'], node2=TF_bridged_DCL[grepl(";", TF_bridged_DCL[,'internal.TF']),'Gene.2'], regulation=TRUE, DCL=TRUE)
		relation.DCL.twoReg <- rbind(relation.DCL.twoReg.half,data.frame(node1=relation.DCL.twoReg.half$node2,node2=relation.DCL.twoReg.half$node1,regulation=TRUE,DCL=TRUE))
		colnames(relation.DCL.twoReg) <- c('node1','node2','regulation','DCL')
		relation.from.internal <- rbind(relation.DCL, relation.DCL.oneReg.1, relation.DCL.oneReg.2, relation.DCL.twoReg)
	}	
		
	relation.from.internal <- rbind(relation.DCL, relation.DCL.oneReg.1, relation.DCL.oneReg.2)
	relation.from.common<-rbind(data.frame(node1=TF_bridged_DCL[,'common.TF'],node2=TF_bridged_DCL[,'Gene.1'],regulation=TRUE,DCL=FALSE),data.frame(node1=TF_bridged_DCL[,'common.TF'],node2=TF_bridged_DCL[,'Gene.2'],regulation=TRUE,DCL=FALSE))
	
	tf.tar.idx<-paste(TF2target_DCL[,'TF'],TF2target_DCL[,'Gene'],sep='; ')
	G1.G2.idx<-paste(relation.from.common[,'node1'],relation.from.common[,'node2'],sep='; ')
	identify.idx<-G1.G2.idx %in% tf.tar.idx
	relation.from.common[identify.idx,'DCL']<-TRUE
	
	relation<-unique(rbind(relation.from.internal,relation.from.common))
	
	#detach(tfbridgedDCL)
	if (nrow(relation)>1000) {
		warning ("the edge number of TF_bridged_DCL(>1000) is too large to display clearly and the program maybe need long time to run.\n")
	}
	node<-unique(c(as.character(relation$node1),as.character(relation$node2)))
	DCG<-unique(c(as.character(unlist(strsplit(TF_bridged_DCL$GeneisDCG,'; '))),as.character(TF_bridged_DCL[TF_bridged_DCL[,'common.TFisDCG'],'comon.TF'])))
	TF<-unique(c(as.character(TF_bridged_DCL[!is.na(TF_bridged_DCL[,'internal.TF']),'internal.TF']),as.character(TF_bridged_DCL[,'common.TF'])))
	notinexprs<-unique(as.character(TF_bridged_DCL[!TF_bridged_DCL[,'common.TFinexprs'],'common.TF']))
	node.classes<-data.frame(TF= node %in% TF, DCG = node %in% DCG, inexprs = node %in% notinexprs)
	rownames(node.classes) <- as.character(node)
		
	f <- graph.data.frame(relation)
	if (length(V(f)) != nrow (node.classes)) {
		stop('oops - why rows of node classes not equal to nodes of graph?\n')
		} else {
		node.classes = node.classes[V(f)$name,]
		V(f)$shape <- 'circle'; V(f)$shape[node.classes$TF] <- 'square'
		V(f)$color <- 'skyblue'; V(f)$color[node.classes$DCG] <- 'pink'; 
		V(f)$color[node.classes$inexprs] <- 'gray';
		V(f)$label.color <- "black";			V(f)$label.cex <- lcex;
		V(f)$size <- vsize
		V(f)$label <- V(f)$name
		E(f)$lty <- 2; E(f)$lty[relation$DCL] <- 1
		E(f)$arrow.mode <- '-'; E(f)$arrow.mode[relation$regulation] <- '>'
		E(f)$arrow.size <- asize
		E(f)$color <- "black"
		E(f)$width <- ewidth
		pdf(figname[2])
		plot(f,layout=layout.fruchterman.reingold)
		dev.off()
		cat("The graph of", figname[2], "has been completed and saved in your working directory.\n")
	}
}


  ###############################################DRplot#########################################################
  ## visualization of DRLs
  # type: a character string to determine plot two kinds of DRLs (default 'both') or only one ('both','TF2target_DCL','TF_bridged_DCL').
  # intgenelist: a list of gene symbols, which contains only one column to display your interesting genes symbol; default is NULL.
  # figname: two character strings of figure names.
  ##############################################################################################################
"DRplot"<-
function(DRsort.res,type=c('both','TF2target_DCL','TF_bridged_DCL')[1],intgenelist=NULL,vsize=5,asize=0.25,lcex=0.3,ewidth=1,figname=c('TF2target_DCL.pdf','TF_bridged_DCL.pdf')){
	TF2target_DCL<-DRsort.res$DCLs[!is.na(DRsort.res$DCLs[,'internal.TF']),c('internal.TF','Gene.1','Gene.2','DCG')]
	TF2target_DCL.1<-TF2target_DCL[TF2target_DCL[,'internal.TF']==TF2target_DCL[,'Gene.1'],c('internal.TF','Gene.2','DCG')]
	colnames(TF2target_DCL.1)<-c('TF','Gene','DCG')
	TF2target_DCL.2<-TF2target_DCL[TF2target_DCL[,'internal.TF']==TF2target_DCL[,'Gene.2'],c('internal.TF','Gene.1','DCG')]
	colnames(TF2target_DCL.2)<-c('TF','Gene','DCG')	
	TF2target_DCL.3<-TF2target_DCL[grep("; ",TF2target_DCL[,'internal.TF']),c('Gene.1','Gene.2','DCG')]
	colnames(TF2target_DCL.3)<-c('TF','Gene','DCG')
	TF2target_DCL.4<-TF2target_DCL[grep("; ",TF2target_DCL[,'internal.TF']),c('Gene.2','Gene.1','DCG')]
	colnames(TF2target_DCL.4)<-c('TF','Gene','DCG')
	TF2target_DCL<-unique(rbind(TF2target_DCL.1,TF2target_DCL.2,TF2target_DCL.3,TF2target_DCL.4))
	TF2target_DCL.int<-TF2target_DCL
	TF_bridged_DCL<-DRsort.res$TF_bridged_DCL[,c('common.TF','Gene.1','Gene.2','internal.TF','common.TFisDCG','common.TFinexprs','GeneisDCG')]
	
	if( !is.null(intgenelist) ){
		target.gene.idx <- (TF_bridged_DCL$Gene.1 %in% as.character(intgenelist$GeneSymbol)) | (TF_bridged_DCL$Gene.2 %in% as.character(intgenelist$GeneSymbol)) | (TF_bridged_DCL$common.TF %in% as.character(intgenelist$GeneSymbol))
		target.gene.idx.tf2target <- (TF2target_DCL$TF %in% as.character(intgenelist$GeneSymbol)) | (TF2target_DCL$Gene %in% as.character(intgenelist$GeneSymbol))
		if ( length(target.gene.idx.tf2target[target.gene.idx.tf2target]) >0 ){
			TF2target_DCL.1<-TF2target_DCL[target.gene.idx.tf2target,]
		} else {
			TF2target_DCL.1 <- NULL
		}
		if ( length(target.gene.idx[target.gene.idx]) >0 ) {
			TF_bridged_DCL.1 <- TF_bridged_DCL[target.gene.idx,]
		} else { TF_bridged_DCL.1 <- NULL
		}
	}
	if ( !is.null(intgenelist) ) {
		TF2target_DCL.int <- unique(TF2target_DCL.1)
		TF_bridged_DCL<-unique(TF_bridged_DCL.1)
		#if ( nrow(DCListf2target.int)==0 & nrow(tfbridgedDCL.int)==0 ) stop ('There are no interesting genes in TF2target_DCL nor in TF_bridged_DCL. \n')
	}
	
	if(type=='both'){
		visTF2target_DCL(TF2target_DCL.int,vsize=vsize,asize=asize,lcex=lcex,ewidth=ewidth,figname)
		visTF_bridged_DCL(TF2target_DCL,TF_bridged_DCL,vsize=vsize,asize=asize,lcex=lcex,ewidth=ewidth,figname)
	}
	if(type=='TF2target_DCL'){
		visTF2target_DCL(TF2target_DCL.int,vsize=vsize,asize=asize,lcex=lcex,ewidth=ewidth,figname)
	}
	if(type=='TF_bridged_DCL'){
		visTF_bridged_DCL(TF2target_DCL,TF_bridged_DCL,vsize=vsize,asize=asize,lcex=lcex,ewidth=ewidth,figname)
	}
	
	list(TF2target_DCL=TF2target_DCL,TF_bridged_DCL=TF_bridged_DCL)
	
}

  ##########################################DRrank##################################################################
  ## Ranking Regulators 
  # tf: A data frame with 215 TFs.
  # tf2target: a data frame or matrix for regulator-to-target interaction pairs.
  # exprs_design: a data frame or matrix for displaying microarray experimant design.
  ##################################################################################################################
"DRrank"<-
function(exprs,exprs.1,exprs.2,tf,tf2target,exprs_design,p.value=0.05,DRsort.res){
############################################RIF#####################################
	reg.exprs.1<-unique(merge(exprs.1,tf,by.x="row.names",by.y="TF"))##取tf的表达谱##
	rownames(reg.exprs.1)<-reg.exprs.1[,1]
	RIF_reg_rank<-data.frame(TF=reg.exprs.1[,'Row.names'],RIFscore=FALSE)
	reg.exprs.1<-reg.exprs.1[,2:ncol(reg.exprs.1)]

	reg.exprs.2<-unique(merge(exprs.2,tf,by.x="row.names",by.y="TF"))
	rownames(reg.exprs.2)<-reg.exprs.2[,1]
	reg.exprs.2<-reg.exprs.2[,2:ncol(reg.exprs.2)]

	library(limma)
	fit<-lmFit(exprs,exprs_design)
	fit<-eBayes(fit)
	DERes<-topTable(fit,coef=colnames(exprs_design)[2],number=length(fit))
	DERes_adjp0.05<-DERes[DERes[,'adj.P.Val']<p.value,c('ID','adj.P.Val')]##取fdr即adjust p value小于0.05的DE##

	DERes_adjp0.05.exprs.1<-unique(merge(exprs.1,DERes_adjp0.05,by.x="row.names",by.y="ID",all.y=T))##取DE的表达谱##
	DERes_adjp0.05.exprs.2<-unique(merge(exprs.2,DERes_adjp0.05,by.x="row.names",by.y="ID",all.y=T))
	rownames(DERes_adjp0.05.exprs.1)<-DERes_adjp0.05.exprs.1[,'Row.names']
	ncol_DERes_adjp0.05.exprs.1<-ncol(DERes_adjp0.05.exprs.1)-1
	DERes_adjp0.05.exprs.1<-DERes_adjp0.05.exprs.1[,2:ncol_DERes_adjp0.05.exprs.1]
	mean.DERes_adjp0.05.exprs.1<-apply(DERes_adjp0.05.exprs.1,1,mean)##取DE表达平均值##
	tmp.1<-names(mean.DERes_adjp0.05.exprs.1)
	mean.DERes_adjp0.05.exprs.1.dataframe<-data.frame(tmp.1,mean.DERes_adjp0.05.exprs.1)

	DERes_adjp0.05.exprs.2<-unique(merge(exprs.2,DERes_adjp0.05,by.x="row.names",by.y="ID"))
	rownames(DERes_adjp0.05.exprs.2)<-DERes_adjp0.05.exprs.2[,'Row.names']
	ncol_DERes_adjp0.05.exprs.2<-ncol(DERes_adjp0.05.exprs.2)-1
	DERes_adjp0.05.exprs.2<-DERes_adjp0.05.exprs.2[,2:ncol_DERes_adjp0.05.exprs.2]
	mean.DERes_adjp0.05.exprs.2<-apply(DERes_adjp0.05.exprs.2,1,mean)
	tmp.2<-names(mean.DERes_adjp0.05.exprs.2)
	mean.DERes_adjp0.05.exprs.2.dataframe<-data.frame(tmp.2,mean.DERes_adjp0.05.exprs.2)

	length_DERes_adjp0.05.exprs<-nrow(mean.DERes_adjp0.05.exprs.1.dataframe)

	nreg.exprs<-nrow(reg.exprs.1)
	e1<-mean.DERes_adjp0.05.exprs.1.dataframe##[,'mean.DERes_adjp0.05.exprs.1']
	e1<-e1[,'mean.DERes_adjp0.05.exprs.1']

	e2<-mean.DERes_adjp0.05.exprs.2.dataframe##[,'mean.DERes_adjp0.05.exprs.2']
	e2<-e2[,'mean.DERes_adjp0.05.exprs.2']
	#e2<-e2[,'mean.DERes_adjp0.05.exprs.2']
	rownames(RIF_reg_rank)<-RIF_reg_rank[,1]

	for( i in 1:nreg.exprs){
		#reg.exprs.1[rownames(RIF_reg_rank[i,]),]<-as.numeric(reg.exprs.1[rownames(RIF_reg_rank[i,]),])
		cor1<-cor(t(DERes_adjp0.05.exprs.1),t(reg.exprs.1[rownames(RIF_reg_rank)[i],]),method="pearson",use="pairwise.complete.obs")
		#reg.exprs.2[rownames(RIF_reg_rank[i,]),]<-as.numeric(reg.exprs.2[rownames(RIF_reg_rank[i,]),])
		cor2<-cor(t(DERes_adjp0.05.exprs.2),t(reg.exprs.2[rownames(RIF_reg_rank)[i],]),method="pearson",use="pairwise.complete.obs")
		preRIF<-(e1*cor1)^2-(e2*cor2)^2
		preRIF<-apply(preRIF,2,sum)
		preRIFdividlength<-preRIF/length_DERes_adjp0.05.exprs
		RIF_reg_rank[i,2]<-preRIFdividlength[[names(preRIFdividlength)]]
	}

	#colnames(RIF_reg_rank)<-c("Gene","RIF")
	RIF_reg_rank<-RIF_reg_rank[order(-abs(as.numeric(RIF_reg_rank[,'RIFscore']))),]
	RIF_reg_rank<-data.frame(TF=RIF_reg_rank$TF, RIF_score=RIF_reg_rank$RIFscore, RIF_rank=matrix(1:nrow(RIF_reg_rank),nrow(RIF_reg_rank),1) )
	
#################################################TED###############################################
	gem_tf2target<-merge(tf2target,exprs,by.x='gene',by.y='row.names')
	gem_tf2target<-unique(gem_tf2target[,c('TF','gene')])
	gem_tf2target_rep<-table(gem_tf2target[,'TF']); gem_tf2target_rep<-as.data.frame(gem_tf2target_rep); colnames(gem_tf2target_rep)<-c('TF','gem_num')
	tf2target_DCG_rep<-table(DRsort.res$DCG2TF[,'TF']); tf2target_DCG_rep<-as.data.frame(tf2target_DCG_rep); colnames(tf2target_DCG_rep)<-c('TF','DCG_num')
	gem_tf2target_DCG_num<-merge(gem_tf2target_rep,tf2target_DCG_rep,by.x='TF',by.y='TF',all.x=T)
	N<-length(unique(gem_tf2target[,'gene']))
	K<-length(unique(DRsort.res$DCG2TF[,'DCG']))
	TED_reg_rank<-data.frame(gem_tf2target_DCG_num,TEDscore=0)
	for( i in 1:nrow(TED_reg_rank)){
		TED_reg_rank[i,'TEDscore']<- -log(pbinom(TED_reg_rank[i,'DCG_num'],TED_reg_rank[i,'gem_num'],K/N,lower.tail=F))
	}
	TED_reg_rank<-TED_reg_rank[order(-as.numeric(TED_reg_rank[,'TEDscore'])),]
	TED_reg_rank<-TED_reg_rank[,c('TF','TEDscore')]
#	rownames(TED_reg_rank)<-TED_reg_rank[,'TF']
	TED_reg_rank<-data.frame(TF=TED_reg_rank$TF, TED_score=TED_reg_rank$TEDscore, TED_rank=matrix(1:nrow(TED_reg_rank),nrow(TED_reg_rank),1) )

#######################################################TDD##############################################
	TDD_reg_rank<-data.frame(unique(gem_tf2target[,'TF']),TDDscore=0)
	colnames(TDD_reg_rank)<-c('TF','TDDscore')
#	DRsort.res$TF_bridged_DCL<-as.matrix(DRsort.res$TF_bridged_DCL)
#	TDD_reg_rank<-as.matrix(TDD_reg_rank)
#	gem_tf2target<-as.matrix(gem_tf2target)
	
	for( i in 1:nrow(TDD_reg_rank)) {
		k<-dim(DRsort.res$TF_bridged_DCL[DRsort.res$TF_bridged_DCL[,'common.TF']==TDD_reg_rank[i,'TF'],c(1,2,3)])[1]
		n1<-unique(DRsort.res$TF_bridged_DCL[DRsort.res$TF_bridged_DCL[,'common.TF']==TDD_reg_rank[i,'TF'],'Gene.1'])
		n2<-unique(DRsort.res$TF_bridged_DCL[DRsort.res$TF_bridged_DCL[,'common.TF']==TDD_reg_rank[i,'TF'],'Gene.2'])
		n3<-unique(gem_tf2target[gem_tf2target[,'TF']==TDD_reg_rank[i,'TF'],'gene'])
		n<-length(unique(c(n1,n2,n3)))
		if (k==0){
			TDD_reg_rank[i,'TDDscore']<-0
		} else {
			TDD_reg_rank[i,'TDDscore']<-k/(n*(n-1)/2)
		}
	}
	TDD_reg_rank[,'TDDscore']<-as.numeric(TDD_reg_rank[,'TDDscore'])/as.numeric(max(TDD_reg_rank[,'TDDscore']))
	TDD_reg_rank<-TDD_reg_rank[order(-as.numeric(TDD_reg_rank[,'TDDscore'])),]
#	rownames(TDD_reg_rank)<-TDD_reg_rank[,'TF']
#	TDD_reg_rank<-as.data.frame(TDD_reg_rank)
	TDD_reg_rank<-data.frame(TF=TDD_reg_rank$TF, TDD_score=TDD_reg_rank$TDDscore, TDD_rank=matrix(1:nrow(TDD_reg_rank),nrow(TDD_reg_rank),1) )
	if ( nrow(TED_reg_rank) > nrow(TDD_reg_rank) ){
		rank.1<-merge(TED_reg_rank,TDD_reg_rank,by.x='TF',by.y='TF',all.x=T)
	} else {
		rank.1<-merge(TED_reg_rank,TDD_reg_rank,by.x='TF',by.y='TF',all.y=T)
	}
	if ( nrow(rank.1) > nrow(RIF_reg_rank) ){
		rank<-merge(rank.1,RIF_reg_rank,by.x='TF',by.y='TF',all.x=T)
	} else {
		rank<-merge(rank.1,RIF_reg_rank,by.x='TF',by.y='TF',all.y=T)
	}
	rank<-rank[order(as.numeric(rank[,'TED_rank'])),]
	return(rank)
	
	#list(TED_reg_rank=TED_reg_rank,TDD_reg_rank=TDD_reg_rank,RIF_reg_rank=RIF_reg_rank)
}





