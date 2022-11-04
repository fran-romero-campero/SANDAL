# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
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
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0019752","carboxylic acid metabolic process",4.4462817509531,2.81530558540903,0.392822333478267,0,"carboxylic acid metabolic process"),
c("GO:0001505","regulation of neurotransmitter levels",0.00459326627164577,1.60799551052543,0.964241947889512,0.43940066,"carboxylic acid metabolic process"),
c("GO:0005975","carbohydrate metabolic process",5.18579762068807,1.6703340141206,0.743387428676697,0.16110445,"carboxylic acid metabolic process"),
c("GO:0006418","tRNA aminoacylation for protein translation",0.275595976298746,2.42457563370858,0.315609373818932,0.63160935,"carboxylic acid metabolic process"),
c("GO:0009059","macromolecule biosynthetic process",6.71076202287447,1.35157943607677,0.456831553272343,0.60411289,"carboxylic acid metabolic process"),
c("GO:0009069","serine family amino acid metabolic process",0.307748840200266,1.40894865094696,0.461409737901018,0.66163687,"carboxylic acid metabolic process"),
c("GO:0016051","carbohydrate biosynthetic process",1.55252399981627,2.41550321297496,0.433097466252197,0.13440917,"carboxylic acid metabolic process"),
c("GO:0018130","heterocycle biosynthetic process",3.97317532497359,1.55751055248051,0.476151218744264,0.55323376,"carboxylic acid metabolic process"),
c("GO:0019438","aromatic compound biosynthetic process",4.3360433604336,1.77379656289051,0.471664024160979,0.4569607,"carboxylic acid metabolic process"),
c("GO:0034637","cellular carbohydrate biosynthetic process",1.04267144366359,2.50986423065067,0.410302721526566,0.64288953,"carboxylic acid metabolic process"),
c("GO:0034645","cellular macromolecule biosynthetic process",4.69891139589362,1.40861797880273,0.433015540944422,0.56861054,"carboxylic acid metabolic process"),
c("GO:0035556","intracellular signal transduction",2.69165403518442,1.40218669780612,0.921590800596915,0.2303724,"carboxylic acid metabolic process"),
c("GO:0042133","neurotransmitter metabolic process",3.12829656337943,1.60799551052543,0.723209740062977,0.16260321,"carboxylic acid metabolic process"),
c("GO:0043038","amino acid activation",0.289375775113683,2.29275603613353,0.466851455772504,0.6554175,"carboxylic acid metabolic process"),
c("GO:0044271","cellular nitrogen compound biosynthetic process",6.7842542832208,1.32814619820322,0.437243962045042,0.61895329,"carboxylic acid metabolic process"),
c("GO:0044272","sulfur compound biosynthetic process",0.71654953837674,1.42905551849866,0.57042468616189,0.3665967,"carboxylic acid metabolic process"),
c("GO:0044281","small molecule metabolic process",7.84989205824262,2.77836508763452,0.778102908310183,0.14321761,"carboxylic acid metabolic process"),
c("GO:1901362","organic cyclic compound biosynthetic process",4.94694777456249,1.70472486938813,0.476160358446283,0.54314467,"carboxylic acid metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = T,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  #bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

