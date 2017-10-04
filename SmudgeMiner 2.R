

bias.display <- function ( this.chip, log.standard, mask = 1:length(log.standard), chipname="This Chip"  )
###########################################################
# display several types of QC information for a single chip:
# 1. synthetic chip with range of 16-fold difference - to show obvious flaws & ridges
# 2. ratio of chip to standard with range of 2-fold difference
# 3. local background
# 4. local scale factor
# constructs ratios of probes on this chip to standard chip, then
# smooths ratio data locally, and displays lower and upper end of ratios
#
# Arguments:
#  this.chip	vector of values from one CEL file
#  log.standard	vector of trimmed mean of log2 values from all CEL files in the experiment
#  mask		vector of flags for true coding probes or non-coding control
#
# Author: Mark Reimers
# Version 1.6 -- Development
# Changes since version 1.0
#	May 28 -- changed input to log2 standard
#	Jun 15 -- inverted colors (red means high now), and added colors for off-the-scale values
#	Jun 21 -- changed jpg format to png via bitmap() function
#	Jul 19 -- minor changes to version 1.4
#	Jul 22 -- added summary metric
#	Aug 10 -- added correlation
#	Oct 14 -- added RNA Deg plot to batch -- then took it out
#	May 28 -- changed blocks to allow different chip sizes
#	June 14 -- changed distribution display
#	Aug 3	-- checked errors in CEL files
#
# Last change:Aug 5, 2005
#############################################
# How to use these programs:
# to examine a small number of chips individually
# > library(affy)
# NB. this command will try to load the cdf file for your chip type. If you cannot do this, ask the sysadmin for assistance.
# > cel.data = ReadAffy()
# you may want to trim the colnames( intensity( cel.data ) ) to be convenient
# > log.standard = construct.log.standard( intensity( cel.data ) )
# > mask = construct.mask(  cel.data [,1] )
# pick a chip and identify its index, i
# > bias.display( intensity( cel.data )[,i] , log.standard, mask )
#
#  OR to examine a batch of chips at once, with the results written to PNG files
# > batch.qc ( )
#  if you want high-resolution (large) plots, you may call
# > batch.qc( res="high")
#
# Interpretation of Plots
# The upper left plot shows how the raw scanned image would look, if we could distort the scale
#   to a logarithmic scale. This brings out the detail in the low range (intensity values between 50 & 150);
#   these values would be compressed into the same color on a linear scale, but contain more than
#   half the probes on a typical chip. This plot often shows striations where
#   The white squares are non-coding control probes with no sequence recorded in the Chip Definition file
# The upper right plot represents how this chip differs from the rest of the chips in the experiment
#   Each pixel of the image represents one probe and the color is graded from red to white depending on
#   the ratio of that probe on this chip to the average value for that same probe across all chips
# The lower left plot represents the effective background -- the lowest levels achieved by probes
#   this chip relative to the lowest levels achieved by probes in the same region in other chips
#   in the same experiment. The numbers represented here are the local means of probes at the low end
#   relative to the means of the same probes across all the chips.
# The lower right plot represents the effective sensitivity or local scale factor -- the ratio of probes
#   (that give high signal through most of the experiment) to their values across the other chips
# Since probes for genes are widely distributed, there is no biological reason for spatial patterns
# of departures from a constant log-ratio. Frequently aberrations in the ratio plot (upper-right)
# that affect a large number (but not all) probes scattered among a specific region
# are resolved in the lower plots into either a background aberration (affecting only low signal probes)
# or a change in sensitivity (affecting mostly probes of moderate to high signal).
###########################################################

{
    N = length( log.standard ) # number of probes
    chip.size = sqrt( N )
    w = 10 # width of 'window' for resolving background & scale factor
    m = floor( chip.size / w ) # number of
    offset = chip.size - m*w  # offset
    
    intercepts = slopes = matrix(NA,m,m)
    
    # set up range of indices in top left corner, which will be translated to generate indices for all squares
    range0 =  numeric(w^2)
    for ( i in 1:w) {
        for ( j in 1:w) {
            range0 [i+(j-1)*w] = i+(j-1)*chip.size
        }
    }
    
    this.chip = log2( this.chip)
    bin.points = quantile ( log.standard, seq(10,90,10)/100 )
    
    # Compute intercept and scale for all subsquares (except in offset region)
    for ( i in 1:m )  # i is block row
    for ( j in 1:m )  # j is block column
    {
        range  = range0 + offset + (i-1)*w +(offset +(j-1)*w)*chip.size
        these.probes = this.chip[range]
        std.probes = log.standard[range]
        
        which.probes = which( std.probes < bin.points[2] )  # lowest 20 %
        intercepts[i,j] = mean(these.probes[which.probes] - std.probes[which.probes] , trim = .2)
        which.probes = which( std.probes > bin.points[7] )  # upper 30 %
        slopes[i,j] = mean(these.probes[which.probes] - std.probes[which.probes] , trim = .2)
    }
    loess.control( surface = "interpolate",  statistics = "approximate" ,
    trace.hat = "approximate" , cell = 0.2, iterations = 4)
    rows = matrix( rep(1:m,m), nr=m, nc=m)
    columns = matrix( rep(1:m,m), nr=m, nc=m, byrow=T)
    inputs = data.frame( rows = c(rows), cols = c(columns) )
    loess.obj = loess( c(intercepts) ~ rows + cols, data=inputs , normalize=F , enp.target = length(log.standard) / 1000 ) # fit about one parameter for every 1000  probes
    intercepts.fitted = predict(loess.obj)
    loess.obj = loess( c(slopes) ~ rows + cols, data=inputs , normalize=F , enp.target = length(log.standard)/1000 ) # fit about one parameter for every 1000  probes
    slopes.fitted = predict(loess.obj)
    intercept.mad = mad( intercepts, na.rm=T )
    slopes.mad = mad( slopes, na.rm = T)
    #  cat ("\nIntercept Bias:", intercept.mad, "\tSlopes Bias:", slopes.mad)
    
    med.this.chip = median( this.chip )
    median.log.std = median( log.standard )
    my.heat.colors = c( heat.colors(9), rainbow(5,start=.4,end=.7))[c(14,12:3,1)] # invert usual order so that red means high, and make extremes stand out with big gaps
    
    opar=par( no.readonly = T)
    layout( matrix( c(1:8),nc=2,byrow=T),heights=c(4,1,4,1)) # arrange 2x2 display with legends at bottom
    par(mar=c(1,1,2,1))
    
    #Synthetic image
    synth.image = trim.image( this.chip, 2 )
    synth.image[mask==0] = NA
    synth.image = matrix( synth.image, nr = chip.size, nc=chip.size)[ , chip.size:1 ] # construct matrix out of vector, and invert row order to show image right-side up; by default the image plots column 1 of the intensity matrix (which is actually row 1 of the chip) along the bottom row
    
    plot( x=1:chip.size, y=1:chip.size, type="n",bty="n",	main = "Synthetic Image (log2 scale)",  xlab="",ylab="",xaxt="n",yaxt="n", bty="n")
    rect( .5, .5, chip.size+.5, chip.size+.5, col= "lightblue" )
    image ( x=1:chip.size, y=1:chip.size,
    z=synth.image, add = T, bty="n",
    zlim = c(med.this.chip-2.401, med.this.chip+2.401), # zlim just outside trimmed range
    col=my.heat.colors ) # the scale limits allow the extreme colors to correspond to the off-scale values
    
    
    # Ratio to standard chip
    ratio.image = trim.image( this.chip - log.standard, .5 )
    ratio.image[mask==0] = NA
    ratio.image = matrix( ratio.image, nr = chip.size, nc=chip.size)[ , chip.size:1 ] # invert image
    plot( x=1:chip.size, y=1:chip.size, type="n",bty="n",	main = "Ratio to Standard",  xlab="",ylab="",xaxt="n",yaxt="n", bty="n")
    image ( x=1:chip.size, y=1:chip.size,
    z=ratio.image ,  add = T, bty="n",
    zlim = c(med.this.chip-median.log.std-.601, med.this.chip-median.log.std+.601),
    col=my.heat.colors ) # the scale limits allow the extreme colors to correspond to the off-scale values
    
    mm = 2:(N-1)
    xx = this.chip - log.standard
    nbr.cor = cor( xx[mm] , ( xx[mm-1]+xx[mm+1] )/2 )
    mm = (2*chip.size+1):(N-2*chip.size)
    nbr.cor = cor( xx[mm], ( xx[mm+2*chip.size] + xx[mm-2*chip.size])/2) # cor vertically
    #  cat("\nCorrelation of vert nbr. ratios:", nbr.cor)
    
    
    # now add the legends for these top 2 images
    par( mar = c( 3,1,1,1))
    dummy.x <- seq( med.this.chip-2,med.this.chip+2, length=10)
    dummy.z <- matrix(dummy.x,ncol=1)
    image(x=dummy.x, y=1, z=dummy.z, yaxt="n", col=my.heat.colors[2:11])
    # this really should be set up to show the off-scale values
    dummy.x <- 2^seq( med.this.chip-median.log.std -.5, med.this.chip-median.log.std +.5, length=10)
    dummy.z <- matrix(dummy.x,ncol=1)
    image(x=dummy.x, y=1, z=dummy.z, yaxt="n", col=my.heat.colors[2:11])
    
    # record images of intercepts and slopes across array
    par(mar=c(1,1,2,1))
    med.int = median(intercepts, na.rm=T)
    trimmed.intercepts = trim.image( intercepts, .5 )[ , m:1 ] # invert to display right-side up
    trimmed.intercepts[ is.na(trimmed.intercepts) ] = med.int
    # N.B. a value is  missing if no probes in the log standard chip in that square are in the bottom 20%. This happens sometimes but for QC purposes we don't care about these, and filling in a neutral color isn't distracting
    image( x = w*1:m, y = w*1:m, z=trimmed.intercepts, zlim=c(med.int-.601, med.int+.601), main="Local background (log2 scale)", col=my.heat.colors, xaxt="n", yaxt="n")
    
    trimmed.slopes = trim.image( slopes, .5 )[ , m:1 ] # invert image
    med.slopes = median ( slopes, na.rm=T )
    trimmed.slopes[ is.na(trimmed.slopes) ] = med.slopes # see note above for trimmed.intercepts
    image( x = w*1:m, y = w*1:m, z=trimmed.slopes, zlim=c(med.slopes-.601, med.slopes+.601), main="Local scale factor (log2 scale)", col=my.heat.colors , xaxt="n", yaxt="n")
    
    # fit loess on 1 df per 100 probes
    rows = matrix( 1:m, nr=m,nc=m)
    cols = matrix( 1:m, nr=m,nc=m, byrow=T)
    int.loess = loess( c(intercepts) ~ c(rows) + c(cols), deg = 1 , enp.target = N / 100)
    slopes.loess = loess( c(slopes) ~ c(rows) + c(cols), deg = 1, enp.target = N / 100)
    SS.int = 10000* sum( int.loess$residuals^2 )
    SS.slope = 10000* sum( slopes.loess$residuals^2 )
    TSS = N*var( this.chip - log.standard )
    #  R.sq.int = 1 - SS.int / TSS
    #  R.sq.slope = 1 - SS.slope / TSS
    
    # now add the legends
    par( mar = c( 3,1,1,1))
    dummy.x <- seq( med.int-.5, med.int+.5, length=10)
    dummy.z <- matrix(dummy.x, ncol=1)
    image(x=dummy.x, y=1, z=dummy.z, yaxt="n", col=my.heat.colors[2:11])
    
    dummy.x <- seq(med.slopes -.5, med.slopes + .5, length=10 )
    dummy.z <- matrix(dummy.x,ncol=1)
    image(x=dummy.x, y=1, z=dummy.z, yaxt="n", col=my.heat.colors[2:11])
    
    #  par(opar)
    
    # calculate row means (relevant for chips with similar probes in rows )
    these.row.means = apply ( matrix( this.chip, nr=chip.size), 2, mean, trim=.2 )
    std.row.means   = apply ( matrix( log.standard,  nr=chip.size), 2, mean, trim=.2 )
    if ( std.row.means[100] > std.row.means[101] )
    { pm.rows =  seq(2, chip.size, 2) ; mm.rows = pm.rows -1 }  # awkward way to find PM rows
    else
    { pm.rows = seq(1, chip.size-1, 2); mm.rows = pm.rows +1  }
    row.cors = cor( these.row.means[pm.rows], std.row.means[pm.rows], use="pair")
    cat("\t", nbr.cor,"\t",row.cors,"\t", 2^( mean( these.row.means[pm.rows] - these.row.means[mm.rows] , trim=.2 ) ) )
}
trim.image = function( data, cutoff )
# this function prepares a 'trimmed' data set suitable for image plots, where values further than cutoff are set to a value just beyond cutoff, to facilitate visualizing extreme values
{
    m = median( data, na.rm=T )
    data[ data < m - cutoff ] = m - 1.2*cutoff
    data[ data > m + cutoff ] = m + 1.2*cutoff
    data
}

construct.log.standard <- function ( intensity.matrix, trim=.2 )
{
    log.std = apply( log2( intensity.matrix), 1, mean, trim = trim )
    log.std
}

construct.mask <- function ( cel.data.one.chip )
{
    mask = cel.data.one.chip
    intensity( mask )[,1] = 0
    pm( mask ) =  1
    mm( mask ) = -1
    intensity(mask)
}

batch.qc <- function ( res = c("low", "high"), use.mask=F, bias=T, ratios=T, deg=F )
# Reads in a set of CEL files using the ReadAffy() function
# Then computes a standard chip image, and a mask
# Then goes through the whole set, to produce a bias plot for each chip
# Then produces a log-ratios plot for each chip
# Then produces RNA degradation plot
# Argument res refers to the resolution of the png files made
{
    if( !require(affy)) warning("Did not complete loading package affy")
    cel.data=ReadAffy()
    colnames(intensity(cel.data)) = sub(".*/","",colnames(intensity(cel.data))) # remove leading path
    colnames(intensity(cel.data)) = sub(".CEL","",colnames(intensity(cel.data))) # remove "CEL" or "cel"
    colnames(intensity(cel.data)) = sub(".cel","",colnames(intensity(cel.data))) #
    colnames(intensity(cel.data)) = gsub(" ",".",colnames(intensity(cel.data))) # collapse any spaces in filenames
    
    # Construct mask file flagging non-coding probes
    if ( use.mask ) mask = construct.mask( cel.data[,1])
    else mask= rep( 1, dim(intensity(cel.data))[1])
    
    # check for bad values
    bads <- which( apply( intensity( cel.data) <= 0, 2, sum ) > 0 )
    if ( length( bads ) > 0 )
    stop( paste( "Cel files", colnames( intensity( cel.data ))[bads]," have negative entries"))
    
    # construct a synthetic standard chip (unless one exists already)
    if ( !exists( "log.standard" ) )
    log.standard = apply( log2( intensity(cel.data)),1, mean, trim=.2, na.rm=T)
    
    cat("\nID\tNeighbor Corr.\tRow Avg. Corr.\tPM to MM ratio")
    if ( bias )
    for ( i in 1:length(cel.data))
    {
        if ( is.null( try( png( paste( colnames(intensity(cel.data))[i],".png",sep=""), res=ifelse( res=="high", 288, 144)),silent=T))){}
        else
        bitmap( paste( colnames(intensity(cel.data))[i],".png",sep=""), res=ifelse( res=="high", 288, 144)) # make hi resolution images on demand
        cat( paste("\n", colnames(intensity(cel.data))[i] ) )
        bias.display( intensity(cel.data)[,i], log.standard, mask, colnames(intensity(cel.data))[i] )
        dev.off()
    }
    if ( ratios )
    show.ratio.densities ( cel.data, log.standard )
    if ( deg )
    plot.RNA.deg ( cel.data )
}

#### show.ratio.distributions
show.ratio.densities <- function( cel.data, log.standard )
{
    N=dim(intensity(cel.data))[2]
    for ( ii in 1:N )
    {
        this.chip = log2( intensity( cel.data)[,ii])
        ratios <- this.chip - log.standard - median(this.chip - log.standard)
        if ( ii %% 8  == 1 ) {  # start a new graph
            if ( is.null( try( png(paste("chips",ii,ii+7,"ratio.distn.png",sep=".")),silent=T))){}
            else
            bitmap(paste("chips", ii, ii+7, "rna.ratio.distn.png", sep=".") )
            plot( density( ratios), xlim = c(-1,1),main="Ratios to Standard",ylim=c(0,2.5))
            legend(-1, par("usr")[4], colnames(intensity(cel.data))[ii:min(ii+7,N)], lty=1, col=ii:min(ii+7,N),bty="n",cex=.7)
        }
        lines( density( ratios ), col=ii )
        if ( ii %% 8 == 0 || ii == N )
        {
            dev.off() # close off a graph
        }
    }
}

plot.RNA.deg <- function( cel.data) {
    library(affy)
    # will try to load chip CDF environment -- if this fails, need to catch it
    deg <- try(AffyRNAdeg(cel.data, log.it=T),silent=F)
    if ( class( deg ) == "try-error" ) warning("Problem obtaining RNA degradation")
    means = deg$means.by.number  # mean signal at each position
    nprobes = dim(means)[2]
    means <- sweep( means, 1, apply( means, 1, median), "-")
    M = dim(intensity(cel.data))[2]
    for ( ii in 1:M)  
    { 
        if ( ii %% 8  == 1 ) 
        {  # start a new graph
            upper = min( ii+7, M)
            if ( is.null( try( png(paste("chips", ii, upper, "rna.position.profile.png", sep=".") ),silent=T))){}
            else
            bitmap(paste("chips", ii, upper, "rna.position.profile.png", sep=".") )
            
            matplot( t(means[ii:upper,]), type="l", xlim=c(1,nprobes), xlab="Probe Position", ylab="Mean difference in log2 Signal", main="RNA Degradation Plot",
            lty=1, col=1:6, lwd=2)
            legend( nprobes/2, 0, colnames(intensity(cel.data))[ii:upper], lty=1, col=1:6, bty="n", cex=.7 ,ncol=2)
            dev.off()
        }
    }
}



