#SAC model
#mm is model.matrix to speed up computations
#
#using family="binomial" with a probit link then we will be fitting a
#spatial probit model 
#
sac.inla <- function(formula, d, W.rho, W.lambda, rho, lambda, mmatrix = NULL, 
  improve = TRUE, impacts = FALSE, fhyper = NULL, probit = FALSE, ...)
{
	#require(INLA)

	IrhoW <- Matrix::Diagonal(nrow(W.rho)) - rho * W.rho
	IlambdaW <- Matrix::Diagonal(nrow(W.lambda)) - lambda * W.lambda
	#IrhoW2<-t(IrhoW)%*%IrhoW
	IrhoW2 <- Matrix::crossprod(IlambdaW %*% IrhoW)

        #environment(formula)<-environment()
        #This is a fix to be able to use improve=TRUE later
        #environment(formula)<-environment()
        assign("IrhoW2", IrhoW2, environment(formula) )

	if(is.null(mmatrix))
		mmatrix<-model.matrix(formula, d)

	mm <- as.data.frame(as.matrix(solve(IrhoW) %*% mmatrix))
	names(mm) <- paste("x", 1:ncol(mm), sep = "")
	xnam <- names(mm)

	#mm$z<-1:nrow(mm)

	d2 <- cbind(d, mm)

	fmla <- paste(as.character(formula)[2], "~ -1+", paste(xnam, collapse= "+"))
	if(is.null(fhyper))
	fmla<-paste(fmla, "+f(idx, model=\"generic0\", Cmatrix=IrhoW2)", sep="")
	else
	fmla<-paste(fmla, "+f(idx, model=\"generic0\", Cmatrix=IrhoW2, hyper=fhyper)", sep="")
	fmla<-as.formula(fmla)


	res <- INLA::inla(fmla, data=d2, ...)


	if(improve)
		res <- INLA::inla.rerun(res)#inla.hyperpar(res, diff.logdens=20)

	#Compute log-determinat to correct the marginal-loglikelihood
	res$logdet<-as.numeric(Matrix::determinant(IrhoW2)$modulus)
	res$mlik<-res$mlik+res$logdet/2

	res$impacts<-FALSE
	if(impacts)
	{
		res$impacts<-TRUE
		#Compute weights for impacts
		if(!probit)
		{
		wtotal <- 1/(1 - rho)
		wdirect <- trIrhoWinv(W.rho, rho)/nrow(W.rho)
		}
		else
		{
		Df <- dnorm(res$summary.linear.predictor[,1])
		wtotal <- mean(Df) * 1/(1 - rho)
		wdirect <- trIrhoWinv(W.rho, rho, Df = Matrix::Diagonal(x = Df))/nrow(W.rho)
		}
		windirect <- wtotal - wdirect

		#
		#TOTAL IMPACTS
		#
		#Add summary.
		#Total impacts
		res$summary.total.impacts <- wtotal * res$summary.fixed[rownames(res$summary.fixed) != "(Intercept)",1:2]
		#Direct impacts
		res$summary.direct.impacts <- wdirect*res$summary.fixed[rownames(res$summary.fixed) != "(Intercept)",1:2]
		#INDirect impacts
		res$summary.indirect.impacts <- windirect*res$summary.fixed[rownames(res$summary.fixed) != "(Intercept)",1:2]


		#Add marginals
		res$marginals.total.impacts<-res$marginals.fixed[ names(res$marginals.fixed) !="(Intercept)"]
		res$marginals.direct.impacts<-res$marginals.fixed[ names(res$marginals.fixed) !="(Intercept)"]
		res$marginals.indirect.impacts<-res$marginals.fixed[ names(res$marginals.fixed) !="(Intercept)"]

		for(i in 1:length(res$marginals.total.impacts))
		{
	xx<-res$marginals.total.impacts[[i]]
	res$marginals.total.impacts[[i]]<-rescalemarg(xx, wtotal)
	res$marginals.direct.impacts[[i]]<-rescalemarg(xx, wdirect)
	res$marginals.indirect.impacts[[i]]<-rescalemarg(xx, windirect)
		}
	}

	return(res)
}
