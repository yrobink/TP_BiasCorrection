
library(RColorBrewer)

basic_bias = function( X , Y , lZ = NULL )##{{{
{
	method = list( mean = mean , sd = sd )
	for( i in 1:2 )
	{
		sX = apply( X , 2 , method[[i]] )
		sY = apply( Y , 2 , method[[i]] )
		
		biasd  = sY - sX
		biasr  = 100 * (sY - sX) / sY
		namesd = "Dif. X"
		namesr = "Rel. X"
		
		if( !is.null(lZ) )
		{
			for( j in 1:length(lZ) )
			{
				sZ    = apply( lZ[[j]] , 2 , method[[i]] )
				
				biasd = rbind( biasd , sY - sZ )
				biasr = rbind( biasr , 100 * (sY - sZ) / sY )
				
				namesd = base::c( namesd , paste( "Dif." , names(lZ)[[j]] ) )
				namesr = base::c( namesr , paste( "Rel." , names(lZ)[[j]] ) )
			}
		}
		
		bias = rbind( biasd , round( biasr , 2 ) )
		rownames(bias) = base::c( namesd , namesr )
		
		cat( paste( names(method)[i] , "\n" ) )
		print(bias)
		cat("\n")
	}
	
	for( p in base::c(0.01,0.1,0.25,0.5,0.75,0.9,0.99) )
	{
		sX = apply( X , 2 , quantile , probs = p )
		sY = apply( Y , 2 , quantile , probs = p )
		
		biasd  = sY - sX
		biasr  = 100 * (sY - sX) / sY
		namesd = "Dif. X"
		namesr = "Rel. X"
		
		if( !is.null(lZ) )
		{
			for( j in 1:length(lZ) )
			{
				sZ    = apply( lZ[[j]] , 2 , quantile , probs = p )
				
				biasd = rbind( biasd , sY - sZ )
				biasr = rbind( biasr , 100 * (sY - sZ) / sY )
				
				namesd = base::c( namesd , paste( "Dif." , names(lZ)[[j]] ) )
				namesr = base::c( namesr , paste( "Rel." , names(lZ)[[j]] ) )
			}
		}
		
		bias = rbind( biasd , round( biasr , 2 ) )
		rownames(bias) = base::c( namesd , namesr )
		
		
		cat( paste( "Quantile" , p , "\n" ) )
		print(bias)
		cat("\n")
	}
}
##}}}

plot_density = function( X , Y , Z = NULL )##{{{
{
	if( is.null(Z) )
	{
		bw = SBCK::bin_width_estimator(list(X,Y))
	}
	else
	{
		bw = SBCK::bin_width_estimator(list(X,Y,Z))
	}
	names(bw) = colnames(Y0)
	
	par( mfrow = base::c(2,3) )
	for( v in colnames(X) )
	{
		plot(  density( Y[,v] , bw = bw[v] ) , col = "blue" , main = v , xlab = "" )
		abline( v = min(Y[,v]) , col = "blue" , lty = 2 )
		abline( v = max(Y[,v]) , col = "blue" , lty = 2 )
		lines( density( X[,v] , bw = bw[v] ) , col = "red"  )
		abline( v = min(X[,v]) , col = "red" , lty = 2 )
		abline( v = max(X[,v]) , col = "red" , lty = 2 )
		if( !is.null(Z) )
		{
			lines( density( Z[,v] , bw = bw[v] ) , col = "green"  )
			abline( v = min(Z[,v]) , col = "green" , lty = 2 )
			abline( v = max(Z[,v]) , col = "green" , lty = 2 )
		}
	}
	leg = c("Y","X")
	col = c("red","blue")
	lty = c(1,1)
	if( !is.null(Z) )
	{
		leg = c(leg,"Z")
		col = c(col,"green")
		lty = c(lty,1)
	}
	legend( "topright" , legend = leg , col = col , lty = lty , cex = 0.8 , box.lty = 1 , box.col = "black" )
}
##}}}

delta_projection = function( X0 , X1 , lZ0 , lZ1 )##{{{
{
	method = list( mean = mean , sd = sd )
	for( i in 1:2 )
	{
		sX0 = apply( X0 , 2 , method[[i]] )
		sX1 = apply( X1 , 2 , method[[i]] )
		
		nZ     = length(lZ0)
		namesd = "Dif. X"
		namesr = "Rel. X"
		biasd  = sX1 - sX0
		biasr  = 100 * (sX1 - sX0) / sX0
		for( j in 1:nZ )
		{
			sZ0 = apply( lZ0[[j]] , 2 , method[[i]] )
			sZ1 = apply( lZ1[[j]] , 2 , method[[i]] )
			
			biasd = rbind( biasd , sZ1 - sZ0 )
			biasr = rbind( biasr , 100 * (sZ1 - sZ0) / sZ0 )
			
			namesd = base::c( namesd , paste( "Dif." , names(lZ0)[[j]] ) )
			namesr = base::c( namesr , paste( "Rel." , names(lZ0)[[j]] ) )
		}
		
		
		bias = rbind( biasd , biasr )
		rownames(bias) = base::c( namesd , namesr )
		
		cat( paste( names(method)[i] , "\n" ) )
		print(bias)
		cat("\n")
	}
	
	for( p in base::c(0.01,0.1,0.25,0.5,0.75,0.9,0.99) )
	{
		sX0 = apply( X0 , 2 , quantile , probs = p )
		sX1 = apply( X1 , 2 , quantile , probs = p )
		
		nZ     = length(lZ0)
		namesd = "Dif. X"
		namesr = "Rel. X"
		biasd  = sX1 - sX0
		biasr  = 100 * (sX1 - sX0) / sX0
		for( j in 1:nZ )
		{
			sZ0 = apply( lZ0[[j]] , 2 , quantile , probs = p )
			sZ1 = apply( lZ1[[j]] , 2 , quantile , probs = p )
			
			biasd = rbind( biasd , sZ1 - sZ0 )
			biasr = rbind( biasr , 100 * (sZ1 - sZ0) / sZ0 )
			
			namesd = base::c( namesd , paste( "Dif." , names(lZ0)[[j]] ) )
			namesr = base::c( namesr , paste( "Rel."    , names(lZ0)[[j]] ) )
		}
		
		
		bias = rbind( biasd , biasr )
		rownames(bias) = base::c( namesd , namesr )
		
		cat( paste( "Quantile" , p , "\n" ) )
		print(bias)
		cat("\n")
	}
}
##}}}

plot_projection_density = function( X0 , Y0 , X1 , lZ1 )##{{{
{
	## Find bw
	bw = SBCK::bin_width_estimator(c(list(X0),list(Y0),list(X1),lZ1))
	names(bw) = colnames(Y0)
	
	par( mfrow = base::c(2,3) )
	
	if( length(lZ1) < 5 )
	{
		colorsZ = base::c("orange","black","purple","cyan")
	}
	else
	{
		colorsZ = rev( brewer.pal( n = max(length(lZ1) ,5) , name = "YlOrRd" ) )
	}
	for( v in colnames(X0) )
	{
		plot( density( Y0[,v] , bw = bw[v] )  , col = "blue" , main = v , xlab = "" )
		abline( v = min(Y0[,v]) , col = "blue" , lty = 2 )
		abline( v = max(Y0[,v]) , col = "blue" , lty = 2 )
		lines( density(X0[,v]) , col = "red" )
		abline( v = min(X0[,v]) , col = "red" , lty = 2 )
		abline( v = max(X0[,v]) , col = "red" , lty = 2 )
		lines( density(X1[,v]) , col = "#990000" )
		abline( v = min(X1[,v]) , col = "#990000" , lty = 2 )
		abline( v = max(X1[,v]) , col = "#990000" , lty = 2 )
		for( i in 1:length(lZ1) )
		{
			lines( density(lZ1[[i]][,v]) , col = colorsZ[i] )
			abline( v = min(lZ1[[i]][,v]) , col = colorsZ[i] , lty = 2 )
			abline( v = max(lZ1[[i]][,v]) , col = colorsZ[i] , lty = 2 )
		}
	}
	leg = c("Y0","X0","X1",names(lZ1))
	col = c("blue","red","#990000",colorsZ)
	lty = c(c(1,1,1),rep( 1 , length(lZ1) ))
	legend( "topright" , legend = leg , col = col , lty = lty , cex = 0.8 , box.lty = 1 , box.col = "black" )
}
##}}}

plot_density_ETP = function( etpX0 , etpY0 , letpZ0 )##{{{
{
	par( mfrow = base::c(1,1) )
	bw = SBCK::bin_width_estimator(c(list(etpX0),list(etpY0),letpZ0))
	
	if( length(letpZ0) < 5 )
	{
		colorsZ = base::c("orange","black","purple","cyan")
	}
	else
	{
		colorsZ = rev( brewer.pal( n = max(length(lZ1) ,5) , name = "YlOrRd" ) )
	}
	
	plot(  density( etpY0 , bw = bw ) , col = "blue" , ylim = c(0,0.8) , xlab = "" , main = "evspsblpot" )
	abline( v = min(etpY0) , col = "blue" , lty = 2 )
	abline( v = max(etpY0) , col = "blue" , lty = 2 )
	lines( density( etpX0 , bw = bw ) , col = "red" )
	abline( v = min(etpX0) , col = "red" , lty = 2 )
	abline( v = max(etpX0) , col = "red" , lty = 2 )
	for( i in 1:length(letpZ0) )
	{
		etpZ0 = letpZ0[[i]]
		col   = colorsZ[i]
		lines( density( etpZ0 , bw = bw ) , col = col )
		abline( v = min(etpZ0) , col = col , lty = 2 )
		abline( v = max(etpZ0) , col = col , lty = 2 )
	}
	
	leg = c("Y0","X0",names(letpZ0))
	col = c("blue","red",colorsZ)
	lty = rep( 1 , 2 + length(lZ1) )
	legend( "topright" , legend = leg , col = col , lty = lty , cex = 0.8 , box.lty = 1 , box.col = "black" )

}
##}}}

corr_bias = function( X0 , Y0 , lZ0 = NULL , relatif = FALSE , method = "pearson" )##{{{
{
	corr = function(X) { return(stats::cor( X , method = method )) }
	
	cat( "X0\n" )
	if( relatif )
	{
		print(100 * (corr(Y0) - corr(X0)) / corr(Y0))
	}
	else
	{
		print( corr(Y0) - corr(X0) )
	}
	cat("\n")
	
	if( is.null(lZ0) )
		return()
	
	for( i in 1:length(lZ0) )
	{
		cat( names(lZ0)[i] , end = "\n" )
		if( relatif )
		{
			print(100 * (corr(Y0) - corr(lZ0[[i]])) / corr(Y0))
		}
		else
		{
			print( corr(Y0) - corr(lZ0[[i]]) )
		}
		cat("\n")
	}
}
##}}}

