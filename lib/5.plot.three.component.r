#!/usr/bin/env Rscript
########This script will use score file and correlation matrix to generate three component heatmap
args = commandArgs(trailingOnly=TRUE)
library(gplots)
library(dplyr)

inputdir<-args[1]
bin_bed<-args[2]
bin<-args[3]

bin_bed<-read.table(bin_bed,sep='\t',fill=TRUE,stringsAsFactors=FALSE)
inputfile<-list.files(inputdir,pattern="chip.freq.score$")
inputfile_chr<-strsplit(inputfile, split='[.]') %>% sapply(.,function(x){x[1]})
input_heatmap<-paste(inputfile_chr,rep("matrix.heatmap.correlation"),sep='.')

for(i in 1:length(inputfile)){
		binname<-bin_bed[which(bin_bed[,1]==inputfile_chr[i]),]
        	binname<-paste(binname[,1], binname[,4],sep='_')
		data=read.table(paste(inputdir,"/",inputfile[i],sep=''),sep='\t')
		id_total=binname
		id_compo=as.character(data[,1])
		removed_id=setdiff(id_total,id_compo)
		temp_m=matrix(0,length(removed_id),ncol(data))
		temp_m[,1]=removed_id
		data=as.matrix(data)
		data=rbind(data,temp_m)
#####################################GENERATE three component heatmap############################################### 
		png(paste(input_heatmap[i],"png",sep='.'))
		nf <- layout(matrix(c(1,2,3,4),4,1,byrow=TRUE), widths=c(7,7,7,7), heights=c(1,1,1,7), TRUE)
		for(compo in c(3:5)){
			label=rep(0,nrow(data))
			label[which(is.na(data[,compo]))]<-'NA'
			a=which(as.numeric(as.character(data[,compo]))>0)
			b=which(as.numeric(as.character(data[,compo]))<0)
			c=which(as.numeric(as.character(data[,compo]))==0)
			pp=0
			qq=0
			if(length(b)!=0){
				pp=mean(as.numeric(na.omit(data[b,2])))
			}
			if(length(a)!=0){
				qq=mean(as.numeric(na.omit(data[a,2])))
			}
			
			if(pp > qq){
				label[b]=rep("A",length(b))
				label[a]=rep("B",length(a))
			}
			if(pp< qq){
				label[b]=rep("B",length(b))
				label[a]=rep("A",length(a))
			}
			label[c]=rep("C",length(c))
#######################################sort by bin number###########################################################
			new_data=data.frame(rep(1,nrow(data)),data,label)### height ID counts compo1 compo2 compo3 label 
			X=c()
			for(o in 1:nrow(new_data)){
				xs=unlist(strsplit(as.character(new_data[o,2]),"n"))[2]
				X=c(X,xs)
			}
			X=as.numeric(X)
			new_data=data.frame(new_data,X) #add bin_number to new_data
			new_data=new_data[order(as.numeric(new_data[,8])),]
			colnames(new_data)=c("height","bin_ID","chip_counts","compo1","compo2","compo3","Label","bin_number")
			compart_list=as.character(as.vector(new_data$Label))
			compart_list[which(compart_list=="A")]=rep(-1,length(which(compart_list=="A")))
			compart_list[which(compart_list=="B")]=rep(1,length(which(compart_list=="B")))
			compart_list[which(compart_list=="C")]=rep(0,length(which(compart_list=="C")))
			compart_list[which(compart_list=="0")]=rep(0,length(which(compart_list=="0")))
			par(mar=c(1,1,1,1))
			image(matrix(as.numeric(compart_list),length(compart_list),1),xaxt="n",yaxt="n",col=c("seagreen","lightgrey","red"),breaks=c(-2,-0.9,0.9,2))

		}
		
		m=read.table(input_heatmap[i],sep='\t',header=T)
		m=as.matrix(m)
		n=200
        	par(mar=c(1,1,1,1))
		col_one=colorpanel(99,"blue","white")
		first_part=col_one
		middle="lightgrey"
		col_two=colorpanel(99,"white","red")
		new_col=c(col_one,middle,col_two)
		temp_seq=unique(c(seq(-0.25,-0.0001,length=100),seq(0.0001,1,length=100),max(m)))
		colbreaks=temp_seq
		image(m,col=new_col,zlim=c(1, 50),breaks=colbreaks,xaxt="n",yaxt="n")
		dev.off()
}

