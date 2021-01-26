#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# run in the matrix folder : Rscript --vanilla generate.chr.matrix.component.r cell_type bin 
#parameter 1: cell type
#parameter 2: bin
#output: matrix for heatmap; heatmap.png; matrix.component
cell_type=args[1]
bin=args[2]

if(cell_type=="mESC" || cell_type=="public_mESC"){
	bin_path_name=paste("/mnt/NFS/homeGene/JinLab/ssz20/zshanshan/heatmap_for_cis/mESC_name/mm10",bin,sep='/')
}
bin_path_name=paste("/mnt/rstor/genetics/JinLab/ssz20/zshanshan/heatmap_for_cis/hg19_name/",bin,sep='')
############generate matrix###########################
chr_lis=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
#chr_lis=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
#############step 1: generate n*n matrix for every bin############
inputname=paste(chr_lis,"matrix",sep='.')
outputname_one=paste(inputname,"heatmap",sep='.')
outputname_two=paste(inputname,"component",sep='.')
outputname_three=paste(inputname,"heatmap.correlation",sep='.')
outputname_four=paste(inputname,"distance.correlation",sep='.')
#row_path_name=paste(bin_path_name, paste(bin,"bin.name",sep='.'),sep='/')
row_path_name=bin_path_name
for(i in 1:length(inputname)){
	if(file.exists(inputname[i])){
		data=read.table(inputname[i],sep='\t',fill=TRUE)
		row_file_name=paste(chr_lis[i],"bin",sep='.')
		rowname=read.table(paste(row_path_name,row_file_name,sep='/'),sep='\t')
		rowname=as.character(unique(rowname[,1]))
		rowname=paste(rep(chr_lis[i],length(rowname)),rowname,sep='_')
		m=matrix(0,length(rowname),length(rowname))
		colnames(m)=rowname
		rownames(m)=rowname
		data[,1]=as.character(data[,1])
		data[,2]=as.character(data[,2])
		data[,3]=as.numeric(as.character(data[,3]))
		for(ind in 1:nrow(data)){
        		m[data[ind,1],data[ind,2]]=data[ind,3]
        		#print(ind)
        	}
		diag(m)=0
		write.table(m,file=outputname_one[i],sep='\t',quote=FALSE)
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
#write.table(mm,file=outputname_four[i],sep='\t',quote=FALSE)
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
		write.table(m_new,file=outputname_three[i],sep='\t',quote=FALSE)
		write.table(result,file=outputname_two[i],sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
	}
############plot heatmap###############################
#	n=20
#	col <- rgb(1,(n-2):0/(n-1),(n-2):0/(n-1))
#	m=m_heatmap
#	m=m[nrow(m):1,]
#	col_matrix=matrix(0,5,3)
#        colnames(col_matrix)=c("250k","500k","1M")
#        rownames(col_matrix)=c("H1","IMR90","GM12878","mESC","H9")
#        col_matrix[1,]=c(150,600,2400)
#        col_matrix[3,]=c(500,2000,12000)
#        col_matrix[2,]=c(400,1500,6000)
#       col_matrix[4,]=c(250,1000,4000)
#	col_matrix[5,]=c(400,1600,6400)
#        scale=col_matrix[cell_type,bin]
#        png(paste(outputname_one[i],"png",sep='.'))
#	par(mar=c(0,0,0,0))
#        colbreaks <- c(seq(1,scale, length=length(col)),max(m))
#        image(m,col=col,zlim=c(1, 50),breaks=colbreaks,xaxt="n",yaxt="n")
#        dev.off()
	
}
}
###########END############################################
