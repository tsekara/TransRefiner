

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
revigo.data <- rbind(c("GO:0000003","reproduction",0.769,7.6188,1.000,0.000,"reproduction"),
c("GO:0007155","cell adhesion",0.544,9.3787,0.962,0.000,"cell adhesion"),
c("GO:0007610","behavior",0.170,2.2331,0.994,0.000,"behavior"),
c("GO:0007623","circadian rhythm",0.057,2.1219,0.994,0.000,"circadian rhythm"),
c("GO:0016265","(obsolete) death",0.143,6.6162,0.994,0.000,"(obsolete) death"),
c("GO:0022610","biological adhesion",0.550,9.3787,0.994,0.000,"biological adhesion"),
c("GO:0023052","signaling",6.765,19.5609,0.995,0.000,"signaling"),
c("GO:0023060","(obsolete) signal transmission",0.143,19.2361,0.994,0.000,"(obsolete) signal transmission"),
c("GO:0030005","(obsolete) cellular di-, tri-valent inorganic cation homeostasis",0.143,2.0504,0.994,0.000,"(obsolete) cellular di-, tri-valent inorganic cation homeostasis"),
c("GO:0030319","(obsolete) cellular di-, tri-valent inorganic anion homeostasis",0.143,1.8608,0.994,0.000,"(obsolete) cellular di-, tri-valent inorganic anion homeostasis"),
c("GO:0032483","regulation of Rab protein signal transduction",0.001,17.9424,0.780,0.000,"regulation of Rab protein signal transduction"),
c("GO:0048519","negative regulation of biological process",1.984,3.2191,0.815,0.332,"regulation of Rab protein signal transduction"),
c("GO:0048518","positive regulation of biological process",1.744,6.2247,0.817,0.311,"regulation of Rab protein signal transduction"),
c("GO:0034395","regulation of transcription from RNA polymerase II promoter in response to iron",0.001,1.7227,0.832,0.186,"regulation of Rab protein signal transduction"),
c("GO:0006810","transport",17.616,9.7501,0.939,0.454,"regulation of Rab protein signal transduction"),
c("GO:0043687","post-translational protein modification",0.028,2.9054,0.940,0.416,"regulation of Rab protein signal transduction"),
c("GO:0010646","regulation of cell communication",0.929,13.7541,0.772,0.306,"regulation of Rab protein signal transduction"),
c("GO:0032879","regulation of localization",0.726,4.1877,0.803,0.280,"regulation of Rab protein signal transduction"),
c("GO:0006873","cellular ion homeostasis",0.324,6.9112,0.722,0.257,"regulation of Rab protein signal transduction"),
c("GO:0043647","inositol phosphate metabolic process",0.074,1.7394,0.853,0.422,"regulation of Rab protein signal transduction"),
c("GO:0007268","chemical synaptic transmission",0.187,6.5775,0.791,0.445,"regulation of Rab protein signal transduction"),
c("GO:0006650","glycerophospholipid metabolic process",0.543,1.9454,0.835,0.420,"regulation of Rab protein signal transduction"),
c("GO:0042127","regulation of cell proliferation",0.313,7.9683,0.777,0.256,"regulation of Rab protein signal transduction"),
c("GO:0007267","cell-cell signaling",0.407,8.9560,0.801,0.479,"regulation of Rab protein signal transduction"),
c("GO:0050789","regulation of biological process",19.373,32.2834,0.774,0.453,"regulation of Rab protein signal transduction"),
c("GO:0007264","small GTPase mediated signal transduction",0.485,2.5054,0.720,0.487,"regulation of Rab protein signal transduction"),
c("GO:0009605","response to external stimulus",1.370,3.4613,0.916,0.475,"regulation of Rab protein signal transduction"),
c("GO:0016192","vesicle-mediated transport",1.085,2.9605,0.952,0.291,"regulation of Rab protein signal transduction"),
c("GO:0007205","protein kinase C-activating G-protein coupled receptor signaling pathway",0.018,5.0462,0.780,0.367,"regulation of Rab protein signal transduction"),
c("GO:0007167","enzyme linked receptor protein signaling pathway",0.279,7.2395,0.725,0.462,"regulation of Rab protein signal transduction"),
c("GO:0065009","regulation of molecular function",1.726,9.3073,0.824,0.165,"regulation of Rab protein signal transduction"),
c("GO:0042221","response to chemical",3.071,5.5699,0.910,0.251,"regulation of Rab protein signal transduction"),
c("GO:0043363","nucleate erythrocyte differentiation",0.000,2.0991,0.589,0.447,"regulation of Rab protein signal transduction"),
c("GO:0040008","regulation of growth",0.172,2.6934,0.849,0.242,"regulation of Rab protein signal transduction"),
c("GO:0046777","protein autophosphorylation",0.077,1.9270,0.919,0.456,"regulation of Rab protein signal transduction"),
c("GO:0006468","protein phosphorylation",4.137,7.1637,0.893,0.415,"regulation of Rab protein signal transduction"),
c("GO:0060345","spleen trabecula formation",0.000,1.7227,0.683,0.500,"regulation of Rab protein signal transduction"),
c("GO:0051174","regulation of phosphorus metabolic process",0.580,8.9627,0.796,0.274,"regulation of Rab protein signal transduction"),
c("GO:0014070","response to organic cyclic compound",0.227,5.6496,0.910,0.386,"regulation of Rab protein signal transduction"),
c("GO:0032501","multicellular organismal process",2.373,27.3144,0.994,0.000,"multicellular organismal process"),
c("GO:0032502","developmental process",2.812,13.0958,0.994,0.000,"developmental process"),
c("GO:0040011","locomotion",0.997,5.1520,0.994,0.000,"locomotion"),
c("GO:0048610","(obsolete) cellular process involved in reproduction",0.143,8.4953,0.994,0.000,"(obsolete) cellular process involved in reproduction"),
c("GO:0048856","anatomical structure development",2.540,11.7188,0.623,0.000,"anatomical structure development"),
c("GO:0010004","gastrulation involving germ band extension",0.002,2.5576,0.666,0.478,"anatomical structure development"),
c("GO:0008594","photoreceptor cell morphogenesis",0.001,2.9054,0.612,0.462,"anatomical structure development"),
c("GO:0022038","corpus callosum development",0.003,2.3117,0.639,0.493,"anatomical structure development"),
c("GO:0002076","osteoblast development",0.003,2.0991,0.633,0.498,"anatomical structure development"),
c("GO:0050877","neurological system process",0.502,10.7086,0.702,0.000,"neurological system process"),
c("GO:0051179","localization",18.495,11.3322,0.995,0.000,"localization"),
c("GO:0055061","(obsolete) di-, tri-valent inorganic anion homeostasis",0.143,1.8608,0.994,0.000,"(obsolete) di-, tri-valent inorganic anion homeostasis"),
c("GO:0055066","(obsolete) di-, tri-valent inorganic cation homeostasis",0.143,2.2576,0.994,0.000,"(obsolete) di-, tri-valent inorganic cation homeostasis"),
c("GO:0065007","biological regulation",20.498,28.1021,0.995,0.000,"biological regulation"),
c("GO:0000289","nuclear-transcribed mRNA poly(A) tail shortening",0.018,2.9490,0.956,0.016,"nuclear-transcribed mRNA poly(A) tail shortening"),
c("GO:0031330","negative regulation of cellular catabolic process",0.021,2.5035,0.802,0.387,"nuclear-transcribed mRNA poly(A) tail shortening"),
c("GO:0007154","cell communication",7.219,10.7086,0.961,0.029,"cell communication"),
c("GO:0009566","fertilization",0.037,15.7689,0.878,0.043,"fertilization"),
c("GO:0008037","cell recognition",0.067,13.1336,0.881,0.054,"cell recognition"),
c("GO:0048103","somatic stem cell division",0.007,3.1281,0.880,0.111,"cell recognition"),
c("GO:0034329","cell junction assembly",0.044,7.2493,0.825,0.320,"cell recognition"),
c("GO:0034330","cell junction organization",0.056,8.4585,0.842,0.126,"cell recognition"),
c("GO:0050808","synapse organization",0.070,2.1899,0.840,0.394,"cell recognition"),
c("GO:0043062","extracellular structure organization",0.061,2.0297,0.841,0.390,"cell recognition"),
c("GO:0030952","establishment or maintenance of cytoskeleton polarity",0.015,2.5054,0.836,0.299,"cell recognition"),
c("GO:0015908","fatty acid transport",0.043,2.5054,0.888,0.203,"cell recognition"),
c("GO:0017145","stem cell division",0.013,1.9315,0.884,0.470,"cell recognition"),
c("GO:0016197","endosomal transport",0.131,2.4477,0.959,0.221,"cell recognition"),
c("GO:0009988","cell-cell recognition",0.010,14.1205,0.886,0.113,"cell recognition"),
c("GO:0071577","zinc II ion transmembrane transport",0.028,1.9325,0.955,0.348,"cell recognition"),
c("GO:0070997","neuron death",0.053,7.3317,0.863,0.126,"cell recognition"),
c("GO:0010876","lipid localization",0.296,1.9325,0.950,0.500,"cell recognition"),
c("GO:0008219","cell death",0.458,6.6162,0.864,0.187,"cell recognition"),
c("GO:0051674","localization of cell",0.633,4.9084,0.955,0.253,"cell recognition"),
c("GO:0030032","lamellipodium assembly",0.013,2.0687,0.834,0.432,"cell recognition"),
c("GO:0030030","cell projection organization",0.608,8.0105,0.815,0.387,"cell recognition"),
c("GO:0006928","movement of cell or subcellular component",0.973,9.2216,0.856,0.156,"cell recognition"),
c("GO:0048193","Golgi vesicle transport",0.297,1.8362,0.955,0.253,"cell recognition"),
c("GO:0045175","basal protein localization",0.000,2.5054,0.962,0.335,"cell recognition"),
c("GO:0016477","cell migration",0.293,7.4593,0.810,0.179,"cell recognition"),
c("GO:0019695","choline metabolic process",0.016,1.7227,0.990,0.067,"choline metabolism"),
c("GO:0006793","phosphorus metabolic process",13.507,4.2904,0.951,0.070,"phosphorus metabolism"),
c("GO:0042541","hemoglobin biosynthetic process",0.003,2.0991,0.963,0.074,"hemoglobin biosynthesis"),
c("GO:0000271","polysaccharide biosynthetic process",0.532,2.0612,0.898,0.252,"hemoglobin biosynthesis"),
c("GO:0006024","glycosaminoglycan biosynthetic process",0.559,2.0687,0.970,0.167,"hemoglobin biosynthesis"),
c("GO:0020027","hemoglobin metabolic process",0.004,2.5054,0.964,0.075,"hemoglobin metabolism"),
c("GO:0006629","lipid metabolic process",3.522,2.3012,0.894,0.077,"lipid metabolism"),
c("GO:0043112","receptor metabolic process",0.037,2.2048,0.959,0.085,"receptor metabolism"));

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
