plot.df <-
		dfgps_result0 <- object$dfgps_result
		coefficient <- coef.step.dfgps(dfgps_result0, candidate_index, intercept=FALSE, stand.coef= stand.coef)
		tuning2 <- apply(abs(coefficient$coefficient),2,sum)
		maxtuning2 <- max(tuning2)
		tuning2 <- tuning2 / maxtuning2
		#result_tuning <- object$dfgps_result$tuning[candidate_index]
		result_df_bind <- cbind(tuning2,result_df)
		xrange <- c(0,1)
		yrange <- c(0,object$dfgps_result$p)
		plot(result_df_bind,xlab="|beta|/max|beta|",ylab="df",type="l",xlim= xrange,ylim= yrange)
		if(object$alpha==0 &&  object$penalty=="enet"){
			dfZou <- apply(abs(coefficient$coefficient) > object$dfgps_result$delta_t*1.5 ,2,sum)
			par(new=T)
			plot(tuning2,dfZou,ann=F,xlab=",ylab=",type="l",xlim= xrange,ylim= yrange,lty=2)
		}
		
		
		