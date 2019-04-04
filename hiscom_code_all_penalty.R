perform_hiscom=function(mRNA_dat,y,gene_annot,lambda_B,lambda_W,alpha=0,err_rate=0.0001,iter_num){
	num_samp=nrow(mRNA_dat)
	num_miRNA=length(unique(gene_annot$miRNA))
	num_mRNA=nrow(gene_annot)
	miRNA_list=unique(gene_annot$miRNA)
	i=1
	genelist=as.numeric(gene_annot$miRNA==miRNA_list[i])
	if(length(miRNA_list)>1){
	for(i in 2:length(miRNA_list)){
		genelist[which(gene_annot$miRNA==miRNA_list[i])]=i
		}
	}

	library(glmnet)

	W_init=matrix(0,nrow=num_mRNA,ncol=num_miRNA)
	i=1
	a=rnorm(sum(genelist==i))
	if(a[1]<0){a=-a}

	W_init[which(genelist==i),i]=a
	if(length(miRNA_list)>1){
		for(i in 2:length(miRNA_list)){
			a=rnorm(sum(genelist==i))
			if(a[1]<0){a=-a}
			W_init[which(genelist==i),i]=a
		}
	}

	if(sum(y[,1]==0)>0){ 
	exclude=which(y[,1]==0)
	X2=as.matrix(mRNA_dat[-exclude,])
	F_init2=X2%*%W_init
		for(inner_iter in 1:ncol(F_init2)){
			if(length(table(F_init2[,inner_iter]))!=1){
			F_init2[,inner_iter]=scale(F_init2[,inner_iter])
			}
			}
	y2=y[-exclude,]
	}else{
	X2=as.matrix(mRNA_dat)
	F_init2=X2%*%W_init
		for(inner_iter in 1:ncol(F_init2)){
			if(length(table(F_init2[,inner_iter]))!=1){
			F_init2[,inner_iter]=scale(F_init2[,inner_iter])
			}
			}
	F_init2=as.matrix(F_init2)
	y2=y
	}
	glmfit=glmnet(F_init2,y2,family="cox",lambda=lambda_B, alpha=alpha)
	beta_init=rep(0,ncol(W_init))
	beta=as.numeric(glmfit$beta[,1])

	i=1
	XB_init=X2[,which(genelist==i)]*beta[i]
	if(length(miRNA_list)>1){
		for(i in 2:length(miRNA_list)){
			XB_init=cbind(XB_init,X2[,which(genelist==i)]*beta[i])
		}
	}


	glmfit=glmnet(XB_init,y2,family="cox",lambda=lambda_W ,alpha=alpha)



	WW=matrix(0,nrow=num_mRNA,ncol=num_miRNA)

	i=1
	if(glmfit$beta[which(genelist==i)][1]<0){
	WW[which(genelist==i),i]=-glmfit$beta[which(genelist==i)]
	}else{
	WW[which(genelist==i),i]=glmfit$beta[which(genelist==i)]
	}

	if(num_miRNA>1){
	for(i in 2:num_miRNA){
		if(glmfit$beta[which(genelist==i)][1]<0){
			WW[which(genelist==i),i]=-glmfit$beta[which(genelist==i)]
		}else{
			WW[which(genelist==i),i]=glmfit$beta[which(genelist==i)]
		}
	}
	}
 
err1=mean((W_init-WW)^2)
err2=mean((beta_init-beta)^2)
W_init=WW
beta_init=beta

	for(iter in 1:iter_num){

		F_init2=X2%*%W_init
		for(inner_iter in 1:ncol(F_init2)){
			if(length(table(F_init2[,inner_iter]))!=1){
			F_init2[,inner_iter]=scale(F_init2[,inner_iter])
			}
			}
		F_init2=as.matrix(F_init2)
		glmfit=glmnet(F_init2,y2,family="cox",lambda=lambda_B,alpha=alpha)
		beta=as.numeric(glmfit$beta[,1])

		i=1
		XB_init=X2[,which(genelist==i)]*beta[i]
		if(num_miRNA>1){
			for(i in 2:num_miRNA){
			XB_init=cbind(XB_init,X2[,which(genelist==i)]*beta[i])
				}
			}

		glmfit=glmnet(XB_init,y2,family="cox",lambda=lambda_W,alpha=alpha)

		WW=matrix(0,nrow=num_mRNA,ncol=num_miRNA)
	i=1
		if(glmfit$beta[which(genelist==i)][1]<0){
			WW[which(genelist==i),i]=-glmfit$beta[which(genelist==i)]
			}else{
			WW[which(genelist==i),i]=glmfit$beta[which(genelist==i)]
		}

		if(num_miRNA>1){
			for(i in 2:num_miRNA){
				if(glmfit$beta[which(genelist==i)][1]<0){
				WW[which(genelist==i),i]=-glmfit$beta[which(genelist==i)]
				}else{
				WW[which(genelist==i),i]=glmfit$beta[which(genelist==i)]
				}
			}
		}
  
	err1=max(abs(W_init-WW))
	err2=max(abs(beta_init-beta))
	if(err1<=err_rate & err2<=err_rate){
	break
	}else{
	W_init=WW
	beta_init=beta
	cat(paste(iter,"iteration :",err2,"\n"))

		}
	}
result=list(W=WW,beta=beta,err1=err1,err2=err2)
	return(result)
}

