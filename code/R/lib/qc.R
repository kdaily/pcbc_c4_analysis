corr_two_nominal_vars = function(v, dat){
    v=as.character(v)
    x = as.factor(as.character(dat[,v[1]]))
    y = as.factor(as.character(dat[,v[2]]))
    dd = data.frame(x=x,y=y)
    dd = dd[rowSums(is.na(dd))==0,]
    sm = chisq.test(dd$y, dd$x) # get a p-value
    r= assocstats(xtabs(~x+y,data=dd))$cramer # get a R^2 value for categorical vars
    pval = sm$p.value
    list(r,pval)
}

most_variable = function(mat, ntop=100){
    require(statmod)
    # now let's recalculate the most variable genes with the winsorized matrix (wed)
    means = rowMeans(mat); vars = apply(mat,1,var); cv2 <- vars/means^2
    useForFit <- means >= unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
    fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
    afit <- fit$coef["a1tilde"]/means+fit$coef["a0"]
    xg <- exp(seq( min(mat), max(mat), length.out=1000 ))
    vfit <- fit$coef["a1tilde"]/xg+fit$coef["a0"]

    varFitRatio <- vars/(afit*means^2)
    order(varFitRatio,decreasing=T)[1:ntop]
}
