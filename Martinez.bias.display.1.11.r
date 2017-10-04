bias.display <- function ( this.chip, log.standard, chipname="This Chip" )
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
#  this.chip	vector of log2 values from one CEL file
#  log.standard	vector of trimmed mean of log2 values from all CEL files in the experiment
#
# Author: Mark Reimers
# Version 1.11 -- Release
# Last change: December 2010
#############################################
# How to use these programs:
# to examine a small number of chips individually
# > library(affy)
# NB. this command will try to load the cdf file for your chip type. If you cannot do this, ask the sysadmin for assistance.
# > cel.data = ReadAffy()
# you may want to trim the colnames( intensity( cel.data ) ) to be convenient
# > log.standard = construct.log.standard( intensity( cel.data ) )
# pick a chip and identify its index, i
# > bias.display( intensity( cel.data )[,i] , log.standard )
#
#  OR to examine a batch of chips at once, with the results written to PNG files
# > batch.qc ( )
#  if you want high-resolution (large) plots, you may call 
# > batch.qc( res="high")
#
# Interpretation of Plots
# The upper left plot shows the ratios of this chip over the standard intensity versus the standard intensity. Ideally it should be a flat line. Although this sort of deviation can be 'corrected' by quantile or lowess normalization the chips remain aberrant.
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
  # may want to allow for non-square chips
  w = 10 # width of 'window' for resolving background & scale factor
  m = floor( chip.size / w ) # number of windows
  offset = chip.size - m*w  # offset of windows

  intercepts = slopes = matrix(NA,m,m)

  # set up range of indices in top left corner, 
  # these will be translated to generate indices for all windows
  range0 =  numeric(w^2)
  for ( i in 1:w) {
    for ( j in 1:w) {
      range0 [i+(j-1)*w] = i+(j-1)*chip.size
    }
  }

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
  intercept.mad = mean( abs( intercepts - mean(intercepts, trim=.1, na.rm=T)), na.rm=T )
  slopes.mad = mean( abs( slopes - mean( slopes,  trim=.1, na.rm=T)), na.rm = T)

  med.this.chip = median( this.chip ) 
  median.log.std = median( log.standard ) 
  my.colors = c( heat.colors(9), rainbow(5,start=.4,end=.7))[c(14,12:3,1)] 
   # red means high, and make extremes stand out with big gaps    

  opar=par( no.readonly = T)
  layout( matrix( c(1:6),nc=3,byrow=T), heights=c(1,1), widths=c(4,4,1.2)) # arrange 2x2 display with legends at bottom

  # M-A plot to show extent of intensity-dependent bias for this chip; 
  ratios <- this.chip - log.standard
  mm <- 10*(1:floor(N/10)) # select a subsample of probes
  par( mar=c(2,2,2,1))
  plot( log.standard[mm], ratios[mm], col="mediumblue", pch=".", 
   xlab="log2 Standard", ylab="log2 Ratios", ylim=c(-1.2,1.2) )
  title(main="Ratios vs. intensity")
  RI.fit <- lowess( log.standard[mm], ratios[mm])
  quenching.deviation <- mean( abs(RI.fit$y-median(RI.fit$y) ))
  quenching.reg <- lsfit( RI.fit$x, RI.fit$y)$coef[2]
  mm <- 50*(1:floor(N/50)) # select a smaller subsample of probes
  points( log.standard[mm], ratios[mm],col="darkblue",pch=".") # give an impression of density
  abline(h=0,col=2)
  lines( RI.fit, col=7, lwd=2)
  rm(ratios,mm,RI.fit); gc()

  
  par(mar=c(.5,.2,2,.2)) # generally thin margins for the rest of the images

  # Image of ratio to standard across the chip
  ratio.image = trim.image( this.chip - log.standard, .5 )
  ratio.image = matrix( ratio.image, nr = chip.size, nc=chip.size)[ , chip.size:1 ] # invert image 
  plot( x=1:chip.size, y=1:chip.size, type="n",bty="n",	main = "Ratio to Standard",  xlab="",ylab="",xaxt="n",yaxt="n", bty="n")
  image ( x=1:chip.size, y=1:chip.size, 
	z=ratio.image ,  add = T, bty="n",
	zlim = c(med.this.chip-median.log.std-.601, med.this.chip-median.log.std+.601), 
	col=my.colors ) # the scale limits allow the extreme colors to correspond to the off-scale values

  mm = 2:(N-1)
  xx = this.chip - log.standard
  nbr.cor = cor( xx[mm] , ( xx[mm-1]+xx[mm+1] )/2 )
  mm = (2*chip.size+1):(N-2*chip.size)
  nbr.cor = cor( xx[mm], ( xx[mm+2*chip.size] + xx[mm-2*chip.size])/2) # cor vertically
#  cat("\nCorrelation of vert nbr. ratios:", nbr.cor)


  # now add the legend for the images
  par( mar = c(.2,1.4,.2,.2 ))
   # this really should be set up to show the off-scale values
  dummy.y <- 2^seq( med.this.chip-median.log.std -.5, med.this.chip-median.log.std +.5, length=10)
  dummy.z <- matrix(dummy.y,nrow=1)
  image(x=1, y=dummy.y, z=dummy.z, xaxt="n", col=my.colors[2:11])

  # record images of intercepts and slopes across array
  par(mar=c(.5,.2,2,.2)) # generally thin margins for the rest of the images
  med.int = median(intercepts, na.rm=T)
  trimmed.intercepts = trim.image( intercepts, .5 )[ , m:1 ] # invert to display right-side up
  trimmed.intercepts[ is.na(trimmed.intercepts) ] = med.int 
	# N.B. a value is  missing if no probes in the log standard chip in that square are in the bottom 20%. This happens sometimes but for QC purposes we don't care about these, and filling in a neutral color isn't distracting
  image( x = w*1:m, y = w*1:m, z=trimmed.intercepts, zlim=c(med.int-.601, med.int+.601), main="Local background (log2 scale)", col=my.colors, xaxt="n", yaxt="n")

  trimmed.slopes = trim.image( slopes, .5 )[ , m:1 ] # invert image
  med.slopes = median ( slopes, na.rm=T )
  trimmed.slopes[ is.na(trimmed.slopes) ] = med.slopes # see note above for trimmed.intercepts
  image( x = w*1:m, y = w*1:m, z=trimmed.slopes, zlim=c(med.slopes-.601, med.slopes+.601), main="Local scale factor (log2 scale)", col=my.colors , xaxt="n", yaxt="n")

 # now add the legends

# calculate row means (relevant for chips with similar probes in rows )
  these.row.means = apply ( matrix( this.chip, nr=chip.size), 2, mean, trim=.2 )
  std.row.means   = apply ( matrix( log.standard,  nr=chip.size), 2, mean, trim=.2 )
  if ( std.row.means[100] > std.row.means[101] )
	{ pm.rows =  seq(2, chip.size, 2) ; mm.rows = pm.rows -1 }  # awkward way to find PM rows
  else
	{ pm.rows = seq(1, chip.size-1, 2); mm.rows = pm.rows +1  }
  row.cors = cor( these.row.means[pm.rows], std.row.means[pm.rows], use="pair")
  pm.mm.diff <- mean( these.row.means[pm.rows] - these.row.means[mm.rows])
  cat("\t", nbr.cor,"\t",row.cors,"\t", pm.mm.diff ) 
  cat(file="qc.out.txt","\t", nbr.cor,"\t",row.cors,"\t", pm.mm.diff , append=T )
}
trim.image = function( data, cutoff )
# this function prepares a 'trimmed' data set suitable for image plots, where values further than cutoff are set to a value just beyond cutoff, to facilitate visualizing extreme values
{
  m = median( data, na.rm=T )
  data[ data < m - cutoff ] = m - 1.2*cutoff
  data[ data > m + cutoff ] = m + 1.2*cutoff
  data
}

batch.qc <- function ( res = c("low", "high"), bias=T, ratios=F, deg=F ) 
# Reads in a set of CEL files using the ReadAffy() function
# Then computes a standard chip image
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
 colnames(intensity(cel.data)) = sub(".gz","",colnames(intensity(cel.data))) # collapse any spaces in filenames 

 # construct a synthetic standard chip (unless one exists already)
 log.intensity <- log2(intensity(cel.data))
 sample.medians <- apply( log.intensity, 2, median)
 log.intensity <- t( t( log.intensity) - sample.medians + median(sample.medians)) # subtract column medians
 log.standard = apply( log.intensity, 1, mean, trim=.2, na.rm=T)
 
cat("\n\tID\tNbr Corr\tRow Corr\tlog(PM/MM) avg")
cat(file="qc.out.txt","\tID\tNbr Corr\tRow Corr\tlog(PM/MM) avg")
if ( bias ) 
 for ( i in 1:length(cel.data)) 
 {
  if ( is.null( try( png( paste( colnames(intensity(cel.data))[i],".png",sep=""), res=ifelse( res=="high", 288, 144)),silent=T))){}
  else
    bitmap( paste( colnames(intensity(cel.data))[i],".png",sep=""), res=ifelse( res=="high", 288, 144)) # make hi resolution images on demand
  cat( paste("\n", colnames(intensity(cel.data))[i] ) )
  cat( paste("\n", colnames(intensity(cel.data))[i] ),file="qc.out.txt",append=T )  
  bias.display( log.intensity[,i], log.standard, colnames(intensity(cel.data))[i] )
  dev.off()
 }
 if ( ratios )
   show.ratio.densities ( log.intensity, log.standard )
 if ( deg )
   plot.RNA.deg ( cel.data )
}

#### show.ratio.distributions
show.ratio.densities <- function( log.intensity, log.standard )
{ 
  N=dim(log.intensity)[2]
  for ( ii in 1:N ) 
  { 
    this.chip = log.intensity[,ii]
    ratios <- this.chip - log.standard - median(this.chip - log.standard) 
    if ( ii %% 8  == 1 ) {  # start a new graph
	if ( is.null( try( png(paste("chips",ii,ii+7,"ratio.distn.png",sep=".")),silent=T))){}
	else
	  bitmap(paste("chips", ii, ii+7, "rna.ratio.distn.png", sep=".") )	
	plot( density( ratios), xlim = c(-1,1),main="Ratios to Standard",ylim=c(0,2.5))
        legend(-1, par("usr")[4], colnames(log.intensity)[ii:min(ii+7,N)], lty=1, col=ii:min(ii+7,N),bty="n",cex=.7)
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
 if ( class( deg ) == "try-error" ) { warning("Problem obtaining RNA degradation") ; return() }
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


