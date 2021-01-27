#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(dplyr)
#output: matrix for heatmap; heatmap.png; matrix.component
inputdir=args[1]
bin_bed=args[2]
#############step 1: generate n*n matrix for every bin############
inputname<-list.files(inputdir,pattern="^chr")
inputname_chr<-strsplit(inputname, split='[.]') %>% sapply(.,function(x){x[1]})
bin_bed<-read.table(bin_bed,sep='\t',fill=TRUE,stringsAsFactors=FALSE)
for(i in 1:length(inputname)){
	data<-read.table(paste(inputdir,"/",inputname[i],sep=''),sep='\t',fill=TRUE)
	rowname<-bin_bed[which(bin_bed[,1]==inputname_chr[i]),]
	rowname<-paste(rowname[,1], rowname[,4],sep='_')
	m=matrix(0,length(rowname),length(rowname))
	colnames(m)=rowname
	rownames(m)=rowname
	for(ind in 1:nrow(data)){
		m[as.character(data[ind,1]),as.character(data[ind,2])]=as.numeric(as.character(data[ind,3]))
       	}
	diag(m)=0
	write.table(m,file=paste(inputdir,"/",inputname[i],".heatmap",sep=''),sep='\t',quote=FALSE)
	###############step 2: remove 0 bin################################
	xmean<-apply(m,1,mean)
	zero_bin=which(xmean==0)
	mc=m[-zero_bin,-zero_bin]
	###############step 3: correct submatrix####################
	mm=mc
	temp_list=c()
	for(g in 2:(ncol(mm)-1)){
        	xrow=1
	        xcol=g
        	for(h in g:ncol(mm)){
                	temp_list=c(temp_list,mm[xrow,xcol])
	                xrow=xrow+1
        	        xcol=h+1
	        }
        	mean_value=mean(temp_list)
	        xrow=1
        	xcol=g
	        for(h in g:ncol(mm)){
        	        mm[xrow,xcol]=mm[xrow,xcol]/mean_value
                	mm[xcol,xrow]=mm[xrow,xcol]
	                xrow=xrow+1
        	        xcol=h+1
	        }
	        temp_list=c()
	}
        mm[nrow(mm),1]=1
        mm[1,ncol(mm)]=1
	mc=mm
	#############submatrix correlation###################
	mm=mc
	for( ii in 1:(nrow(mc)-1)){
		for(jj in (ii+1):ncol(mc)){
			row=mc[ii,]
			coll=mc[,jj]
			row=row[-c(ii,jj)]
			coll=coll[-c(ii,jj)]
			mm[ii,jj]=cor(row,coll,method="pearson")
			mm[jj,ii]=mm[ii,jj]
		}
	}	
	diag(mm)<-1
	mm[which(is.na(mm))]=rep(0,length(which(is.na(mm))))
	mc=mm
	###############PCA and outplut#####################################
	sub_binname=colnames(mc)
	m_new=matrix(0,nrow(m),ncol(m))
	colnames(m_new)=colnames(m)
	rownames(m_new)=rownames(m)
	for(row in 1:(nrow(mc)-1)){
		for(col in (row+1):ncol(mc)){
			a=sub_binname[row]
			b=sub_binname[col]
			m_new[a,b]=mc[row,col]
			m_new[b,a]=m_new[a,b]
		}
	}
	for(rr in 1:nrow(mc)){
		a=sub_binname[rr]
		m_new[a,a]=mc[rr,rr]	
	}
	q=which(rowSums(mc)==0)
	if(length(q)!=0){
        	mc=mc[-q,]
		mc=mc[,-q]
	}
	if(nrow(mc)!=0){
		p=princomp(mc,cor=TRUE)
		x=loadings(p)
		comp=x[,1:3]# get first 3 component
		ID_name=rownames(x)
		result=data.frame(ID_name,comp)	
		write.table(m_new,file=paste(inputdir,"/",inputname[i],".heatmap.correlation",sep=''),sep='\t',quote=FALSE)
		write.table(result,file=paste(inputdir,"/",inputname[i],".component",sep=''),sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
	}

}
###########END############################################
