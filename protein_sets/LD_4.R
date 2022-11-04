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
revigo.data <- rbind(c("GO:0006790","sulfur compound metabolic process",2.39446939931405,2.62167643682623,0.809333103597943,0,"sulfur compound metabolic process"),
c("GO:0009081","branched-chain amino acid metabolic process",0.462096478442625,1.59911431823676,0.547220002528114,0.07610544,"branched-chain amino acid metabolic process"),
c("GO:0006082","organic acid metabolic process",8.96817105384944,1.53793651051609,0.536279027031364,0.65748781,"branched-chain amino acid metabolic process"),
c("GO:0019752","carboxylic acid metabolic process",8.63921761436941,1.55121699008965,0.469152838804948,0.56696945,"branched-chain amino acid metabolic process"),
c("GO:0043043","peptide biosynthetic process",5.18431313154083,1.96092391527736,0.377102350592567,0.48596431,"branched-chain amino acid metabolic process"),
c("GO:0043603","cellular amide metabolic process",6.72366287092487,2.01267529302303,0.603686072286542,0.11260293,"branched-chain amino acid metabolic process"),
c("GO:0044271","cellular nitrogen compound biosynthetic process",13.6693277614853,1.73422699812461,0.491567666533709,0.57708587,"branched-chain amino acid metabolic process"),
c("GO:1901566","organonitrogen compound biosynthetic process",14.6944990362876,2.11576007748892,0.5113710661501,0.22739788,"branched-chain amino acid metabolic process"),
c("GO:0044281","small molecule metabolic process",15.1962015088436,1.8539727914393,0.837218040594346,0.07044348,"small molecule metabolic process"));

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
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

