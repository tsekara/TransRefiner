

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000087","mitotic M phase",0.004,12.7645,0.981,0.000,"mitotic M phase"),
c("GO:0006396","RNA processing",3.210,25.6253,0.799,0.000,"RNA processing"),
c("GO:0006261","DNA-dependent DNA replication",0.576,9.7496,0.811,0.281,"RNA processing"),
c("GO:0051246","regulation of protein metabolic process",1.551,1.6345,0.800,0.416,"RNA processing"),
c("GO:0022616","DNA strand elongation",0.041,2.6737,0.857,0.420,"RNA processing"),
c("GO:0006366","transcription from RNA polymerase II promoter",1.430,3.5157,0.804,0.401,"RNA processing"),
c("GO:0006368","transcription elongation from RNA polymerase II promoter",0.082,3.9393,0.843,0.293,"RNA processing"),
c("GO:0016070","RNA metabolic process",15.951,22.9431,0.790,0.464,"RNA processing"),
c("GO:0016071","mRNA metabolic process",0.798,5.4283,0.834,0.373,"RNA processing"),
c("GO:0006383","transcription from RNA polymerase III promoter",0.079,2.3778,0.846,0.292,"RNA processing"),
c("GO:0044249","cellular biosynthetic process",30.048,5.0910,0.878,0.188,"RNA processing"),
c("GO:0006354","DNA-templated transcription, elongation",0.202,4.2581,0.835,0.320,"RNA processing"),
c("GO:0044267","cellular protein metabolic process",14.293,5.6655,0.826,0.450,"RNA processing"),
c("GO:0010608","posttranscriptional regulation of gene expression",0.719,4.6635,0.805,0.400,"RNA processing"),
c("GO:0016126","sterol biosynthetic process",0.075,1.5200,0.881,0.156,"RNA processing"),
c("GO:0001510","RNA methylation",0.739,7.8697,0.790,0.469,"RNA processing"),
c("GO:0032446","protein modification by small protein conjugation",0.607,1.7825,0.843,0.459,"RNA processing"),
c("GO:0034470","ncRNA processing",2.222,22.7747,0.769,0.426,"RNA processing"),
c("GO:0010501","RNA secondary structure unwinding",0.025,3.9031,0.873,0.263,"RNA processing"),
c("GO:0031047","gene silencing by RNA",0.089,2.9957,0.743,0.331,"RNA processing"),
c("GO:0000018","regulation of DNA recombination",0.045,3.3233,0.789,0.423,"RNA processing"),
c("GO:0019538","protein metabolic process",18.489,4.8827,0.890,0.339,"RNA processing"),
c("GO:0006517","protein deglycosylation",0.008,2.6946,0.885,0.109,"RNA processing"),
c("GO:0010467","gene expression",19.671,18.8327,0.891,0.222,"RNA processing"),
c("GO:0043412","macromolecule modification",9.785,3.4647,0.901,0.285,"RNA processing"),
c("GO:0051052","regulation of DNA metabolic process",0.245,2.2807,0.764,0.487,"RNA processing"),
c("GO:0006417","regulation of translation",0.692,4.6383,0.741,0.311,"RNA processing"),
c("GO:0040029","regulation of gene expression, epigenetic",0.130,2.1838,0.831,0.342,"RNA processing"),
c("GO:0006401","RNA catabolic process",0.395,4.4067,0.819,0.344,"RNA processing"),
c("GO:0048662","negative regulation of smooth muscle cell proliferation",0.007,1.5157,0.854,0.462,"RNA processing"),
c("GO:0009451","RNA modification",1.778,13.8861,0.794,0.413,"RNA processing"),
c("GO:0034660","ncRNA metabolic process",3.407,21.2967,0.809,0.453,"RNA processing"),
c("GO:0009296","(obsolete) flagellum assembly",0.143,1.7825,0.994,0.000,"(obsolete) flagellum assembly"),
c("GO:0009987","cellular process",63.780,8.1798,0.998,0.000,"cellular process"),
c("GO:0050896","response to stimulus",12.210,2.1574,0.995,0.000,"response to stimulus"),
c("GO:0051179","localization",18.495,1.6655,0.995,0.000,"localization"),
c("GO:0051649","establishment of localization in cell",1.679,10.4976,0.837,0.000,"establishment of localization in cell"),
c("GO:0015031","protein transport",2.251,5.5560,0.848,0.361,"establishment of localization in cell"),
c("GO:0006818","hydrogen transport",1.149,1.6345,0.859,0.332,"establishment of localization in cell"),
c("GO:0006839","mitochondrial transport",0.182,1.6345,0.924,0.274,"establishment of localization in cell"),
c("GO:0016192","vesicle-mediated transport",1.085,2.3152,0.912,0.330,"establishment of localization in cell"),
c("GO:0033036","macromolecule localization",3.030,5.0110,0.905,0.377,"establishment of localization in cell"),
c("GO:0015672","monovalent inorganic cation transport",1.824,1.9245,0.905,0.367,"establishment of localization in cell"),
c("GO:0051641","cellular localization",2.041,10.3615,0.909,0.345,"establishment of localization in cell"),
c("GO:0032259","methylation",3.103,6.1024,0.982,0.020,"methylation"),
c("GO:0009056","catabolic process",4.820,2.1152,0.981,0.021,"catabolism"),
c("GO:0009058","biosynthetic process",31.611,4.8928,0.978,0.033,"biosynthesis"),
c("GO:0006457","protein folding",0.903,1.7773,0.953,0.040,"protein folding"),
c("GO:0007049","cell cycle",1.885,18.2604,0.866,0.044,"cell cycle"),
c("GO:0055086","nucleobase-containing small molecule metabolic process",4.917,1.9172,0.768,0.489,"cell cycle"),
c("GO:0048518","positive regulation of biological process",1.744,1.5784,0.887,0.301,"cell cycle"),
c("GO:0043697","cell dedifferentiation",0.001,1.7825,0.879,0.382,"cell cycle"),
c("GO:0043696","dedifferentiation",0.001,1.7825,0.922,0.345,"cell cycle"),
c("GO:0051301","cell division",1.230,8.4078,0.871,0.223,"cell cycle"),
c("GO:0009116","nucleoside metabolic process",2.917,4.6073,0.762,0.346,"cell cycle"),
c("GO:0042180","cellular ketone metabolic process",0.423,1.5638,0.857,0.356,"cell cycle"),
c("GO:0006730","one-carbon metabolic process",0.328,4.6861,0.863,0.193,"cell cycle"),
c("GO:0010525","regulation of transposition, RNA-mediated",0.004,3.1046,0.865,0.132,"cell cycle"),
c("GO:0045454","cell redox homeostasis",0.861,4.8239,0.802,0.214,"cell cycle"),
c("GO:0007114","cell budding",0.016,1.8182,0.869,0.147,"cell cycle"),
c("GO:0006013","mannose metabolic process",0.031,3.9172,0.919,0.232,"cell cycle"),
c("GO:0007059","chromosome segregation",0.476,5.3809,0.880,0.201,"cell cycle"),
c("GO:0007017","microtubule-based process",0.658,3.7773,0.877,0.208,"cell cycle"),
c("GO:0007520","myoblast fusion",0.011,3.6904,0.854,0.143,"cell cycle"),
c("GO:0022402","cell cycle process",1.053,11.9172,0.808,0.219,"cell cycle"),
c("GO:0044281","small molecule metabolic process",15.138,8.0150,0.910,0.206,"cell cycle"),
c("GO:0030029","actin filament-based process",0.398,5.1726,0.882,0.197,"cell cycle"),
c("GO:0006487","protein N-linked glycosylation",0.076,3.3696,0.788,0.483,"cell cycle"),
c("GO:0001654","eye development",0.105,2.7496,0.883,0.450,"cell cycle"),
c("GO:0033554","cellular response to stress",2.967,11.0297,0.872,0.047,"cellular response to stress"),
c("GO:0006301","postreplication repair",0.023,3.3125,0.822,0.494,"cellular response to stress"),
c("GO:0032008","positive regulation of TOR signaling",0.016,4.0768,0.824,0.347,"cellular response to stress"),
c("GO:0043902","positive regulation of multi-organism process",0.036,2.3152,0.875,0.478,"cellular response to stress"),
c("GO:0008286","insulin receptor signaling pathway",0.024,3.7258,0.817,0.358,"cellular response to stress"),
c("GO:0019985","translesion synthesis",0.012,2.3615,0.819,0.470,"cellular response to stress"),
c("GO:0000196","MAPK cascade involved in cell wall organization or biogenesis",0.002,2.6946,0.776,0.300,"cellular response to stress"),
c("GO:0006996","organelle organization",3.595,19.6289,0.718,0.048,"organelle organization"),
c("GO:0043244","regulation of protein complex disassembly",0.160,3.9393,0.648,0.499,"organelle organization"),
c("GO:0035404","histone-serine phosphorylation",0.006,2.1152,0.725,0.372,"organelle organization"),
c("GO:0001578","microtubule bundle formation",0.027,4.9031,0.695,0.421,"organelle organization"),
c("GO:0030447","filamentous growth",0.026,2.1152,0.939,0.067,"filamentous growth"),
c("GO:0051186","cofactor metabolic process",3.985,1.7825,0.926,0.089,"cofactor metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
