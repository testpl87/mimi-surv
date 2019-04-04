

calc_likelihood=function(W_final,beta_final,Xnew,ynew){
Fnew=Xnew%*%W_final

likelihood_mid=data.frame(time=ynew[,1],status=ynew[,2],exp=exp(Fnew%*%beta_final))

event_time=likelihood_mid$time[likelihood_mid$status==1]
event_time=unique(event_time)
event_time=event_time[order(event_time)]

loglik=0
for(iter in 1:length(event_time)){

r_t=event_time[iter]
loglik=loglik+log(sum(likelihood_mid$exp[likelihood_mid$time==r_t&likelihood_mid$status==1])/sum(likelihood_mid$exp[likelihood_mid$time>=r_t]))
}

return(loglik)
}
