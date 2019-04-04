
library(glmnet)

perform_cv=function(mRNA_dat,y,gene_annot,lambda_list=seq(0.05,0.5,by=0.05),err_rate=0.0001,iter_num=1000,cvnum=5,alpha=0){
	num_samp_all=nrow(mRNA_dat)

	samp=sample(1:num_samp_all,num_samp_all,replace=FALSE)
	a=c(1,as.integer(seq(0,num_samp_all,length.out=(cvnum+1)))[2:(cvnum+1)])
	loglik_list=c()

	for(lambda_iter in 1:length(lambda_list)){
	lambda=lambda_list[lambda_iter]

	cv_loglik_list=c()
	for(cv in 1:cvnum){
	train_dat=samp[-c(a[cv]:a[cv+1])]
	test_dat=samp[c(a[cv]:a[cv+1])]

	ytrain=y[train_dat,]
	xtrain=mRNA_dat[train_dat,]
	ytest=y[test_dat,]
	xtest=mRNA_dat[test_dat,]

	result_inter=perform_hiscom(xtrain,ytrain,gene_annot,lambda,lambda,iter_num=500,alpha=alpha)
	cv_loglik_list=c(cv_loglik_list,calc_likelihood(result_inter$W,	result_inter$beta,xtest,ytest))
	}
	loglik_list=c(loglik_list,sum(cv_loglik_list))
	}
	result=data.frame(lambda_list=lambda_list,loglik_list=loglik_list)
	return(result)
}