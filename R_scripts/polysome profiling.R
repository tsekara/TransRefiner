#!/usr/bin/Rscript

## change things here to alter script output
input.dir <- "/home/tsekara/Downloads/Polysomefiles/";
filenames=dir(input.dir,pattern="*.csv")
for( i in 1:length(filenames) )
{
# name of input file to process
#input.file <- sprintf("%s%s",input.dir,"100315_2_100U.csv");
input.file <- sprintf("%s%s",input.dir,filenames[i]);
sprintf("%s%s",input.dir,filenames[i]);
show.raw <- TRUE; # show raw data in final plot
show.smoothed <- TRUE; #  show smoothed data in final plot
smooth.window.size <- 11; # usually 11, set to 1 if already smooth
rescale.data <- FALSE; # scale raw data to range 0..100%
output.directory <- "/home/tsekara/Downloads/polysome_output/";
## local minima window size in minutes
minima.window.size <- (10/60); # 10 seconds
time.range <- c(0,20);
detector.range <- c(0,100);

## command-line processing
argPos <- 1;
if(length(commandArgs(TRUE)) > 0){
  ## this makes sure that there is no default input file if command
  ## line arguments are used
  input.file <- "";
  output.directory <- "./";
}
while(argPos <= length(commandArgs(TRUE))){
  argument <- commandArgs(TRUE)[argPos];
  if(argument == "-input"){
    input.file <- commandArgs(TRUE)[argPos+1];
    argPos <- argPos + 1;
  }
  else if(argument == "-rescale"){
    rescale.data <- TRUE;
  }
  else if(argument == "-noraw"){
    show.raw <- FALSE;
  }
  else {
    cat("Command not understood:", argument);
  }
  argPos <- argPos + 1;
}

# remove directory from input file name
input.base <- sub("^.*/","",input.file);
# ensure output.directory has one '/' at the end
output.directory <- sub("/$","",output.directory);
output.directory <- sub("$","/",output.directory);

## accessory functions
## getArea -- returns the area between two time points
getArea <- function(startTime, stopTime, input.df){
  if(startTime > stopTime){
    return(data.frame(start.time = startTime,
                      end.time = stopTime,
                      base.area = NA, rect.area = NA, trap.area = NA));
  }
  sample.readings <- subset(input.df, Time >= startTime & Time <= stopTime);
  num.readings <- dim(sample.readings)[1];
  ## the area is calculated by summing up tiny trapezoids between data points
  base.area <- sample.readings$widths * sample.readings$averages;
  rect.area <- sample.readings$widths * (sample.readings$averages - min(sample.readings$averages));
  ## subtract triangle from rectangular base
  trap.area <- sum(rect.area) - sum(sample.readings$widths) *
    abs(sample.readings$averages[1] - sample.readings$averages[num.readings]) / 2;
  return(data.frame(start.time = sample.readings$Time[1],
                    end.time = sample.readings$Time[num.readings],
                    base.area = sum(base.area), rect.area = sum(rect.area), trap.area = trap.area));
}

## runmin -- Determines local minima from a series of points with
## a given window size
## note: window size is in minutes
runmin <- function(points, times, window.size){
  cat("Hunting for local minima...");
  min.locs <- NULL;
  next.loc <- 1;
  unique.points <- 0;
  for (pos in 1:length(points)){
    current.time <- times[pos];
    current.window <-
      which(
        (times >= (current.time - window.size)) &
          (times <= (current.time + window.size)));
    min.points <- which(points[current.window] == min(points[current.window])) + min(current.window) - 1;
    if(length(min.points) < length(current.window)){
      if(points[pos] == min(points[current.window]) &&
        (points[current.window[1]] != points[pos]) &&
        (points[current.window[length(current.window)]] != points[pos])){
        min.locs[next.loc] <- trunc(mean(min.points));
        if(length(unique(min.locs)) > unique.points){
          unique.points <- length(unique(min.locs));
          cat("", unique.points);
        }
        next.loc <- next.loc + 1;
      }
    }
    1  }
  cat(" done\n");
  return(unique(min.locs));
}

plotData <- function(raw = TRUE, smoothed = TRUE, adjust.time = FALSE){
  x.label <- "Time";
  time.values <- data.df$Time;
  if(adjust.time){
    time.values <- data.df$AdjustedTime;
    x.label <- "Time (adjusted)";
  }
  ## print a dummy plot so that grid can be drawn
  plot(x = time.values, y = data.df$Reading, type = "n",
       xlab = x.label, ylab = "Absorbance", xlim = time.range,
       ylim = detector.range);
  legend.names <- NULL;
  if(raw){
    ## actual point data -- raw
    points(x = time.values, y = data.df$Reading, type = "l",
           col = "red");
    legend.names <- c(legend.names, "raw data");
  }
  if(smoothed){
    legend.names <- c(legend.names, "smoothed data", "local minima");
    ## actual point data -- smoothed data
    points(x = time.values, y = data.df$smoothed, type = "l", col = "blue");
    ## local minima
    points(x = time.values[local.minima], y = data.df$smoothed[local.minima],
           type = "o", col = "magenta");
  }
  if(length(legend.names)>1){
    ## legend for plot
    legend("topright",legend=legend.names,
           fill=c("red","blue","magenta")[1:length(legend.names)],
           inset = 0.05, bg = "white");
  }
}

data.df <- read.csv(input.file, row.names = NULL);
data.df <- subset(data.df, Time <= max(time.range));
## order by time to make points more predictable
data.df <- data.df[order(data.df$Time),];
## if using raw voltage readings, multiply by 100
if(max(data.df$Reading) < 2){
  data.df$Reading <- data.df$Reading * 100;
}
## if readings are in seconds, convert to minutes
if(max(data.df$Time) > 120){
  data.df$Time <- data.df$Time / 60;
}
## smooth out by running median of 11 points
data.df$smoothed <- runmed(data.df$Reading, smooth.window.size);
## find local minima with a window size of minima.window.size
local.minima <- unique(c(1,runmin(data.df$smoothed,
                                  data.df$Time,minima.window.size),
                         dim(data.df)[1]));
## normalise to absorbance relative to highest peak, lowest local minima
min.point <- min(data.df$Reading[local.minima]);
max.point <- max(data.df$Reading);
if(rescale.data){
  data.df$Reading <- (data.df$Reading - min.point) / (max.point - min.point) * 100;
}
## smooth data using new values
data.df$smoothed <- runmed(data.df$Reading, smooth.window.size);
smoothed.range <- range(data.df$smoothed);
time.range <- range(data.df$Time);

## work out widths/heights of trapeziums between each point
## (assume baseline = 0)
data.df$widths <- c(0,data.df$Time[-1] - data.df$Time[-dim(data.df)[1]]);
data.df$averages <- c(0,(data.df$smoothed[-1] + data.df$smoothed[-dim(data.df)[1]])/2);
## determine t0 for adjusted time
## t0 is set to the start of the largest jump over 20 readings
readings.diff <- c(data.df$Reading[-(1:20)],rep(0,20)) - data.df$Reading;
maxjump.pos <- which(readings.diff == max(readings.diff))[1];
time.t0 <- data.df$Time[maxjump.pos];
data.df$AdjustedTime <- data.df$Time - time.t0;

pdf(paste(output.directory,input.base,".pdf", sep = ""), width = 9.6,
    height = 5.4);
par(mar = c(5.1,5.1,1,1));
## draw plots -- raw data only
plotData(raw = TRUE, smoothed = FALSE, adjust.time = FALSE);
plotData(raw = TRUE, smoothed = FALSE, adjust.time = TRUE);
plotData(raw = show.raw, smoothed = show.smoothed, adjust.time = FALSE);
## create region labels, and produce area table
areas.df <- NULL;
for(x in 2:length(local.minima)){
  minima.window <- local.minima[x-1]:local.minima[x];
  maxVal <- max(data.df$smoothed[minima.window]);
  xpos = data.df$Time[mean(which(data.df$smoothed[minima.window] == maxVal)) + local.minima[x-1] - 1];
  ypos <- maxVal + (max(smoothed.range) * 0.05);
  if(ypos > (max(smoothed.range))){
    ypos <- (max(smoothed.range) * 0.95);
  }
  text(x = xpos, y = ypos, labels = (x-1), cex = 0.5);
  areas.df <-
    rbind(areas.df,cbind(region = x-1,
                         getArea(data.df$Time[local.minima[x-1]],
                                 data.df$Time[local.minima[x]], data.df)));
}
write.csv(areas.df,file = paste(output.directory,"areas_",input.base,sep=""),
          row.names = FALSE);
dummy <- dev.off(); 
}