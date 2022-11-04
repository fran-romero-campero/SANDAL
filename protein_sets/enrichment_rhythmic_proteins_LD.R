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
revigo.data <- rbind(c("GO:0006108","malate metabolic process",0.0826787928896238,2.66859528060254,0.738841597262421,0,"malate metabolic process"),
c("GO:0009117","nucleotide metabolic process",1.40094621285196,1.4881199296299,0.602091344292803,0.63538194,"malate metabolic process"),
c("GO:0019752","carboxylic acid metabolic process",4.4462817509531,1.91283507084552,0.625175355834492,0.55713715,"malate metabolic process"),
c("GO:0043648","dicarboxylic acid metabolic process",0.436360295806348,1.88845786342838,0.699835371968233,0.66556582,"malate metabolic process"),
c("GO:0055086","nucleobase-containing small molecule metabolic process",1.88323917137476,1.31065579504451,0.631997197433989,0.66191227,"malate metabolic process"),
c("GO:0006260","DNA replication",0.675210141931928,2.56666264813171,0.792365418960641,0.08842434,"DNA replication"),
c("GO:0001505","regulation of neurotransmitter levels",0.00459326627164577,2.38995320520681,0.97092114841269,0.43940066,"DNA replication"),
c("GO:0005975","carbohydrate metabolic process",5.18579762068807,1.87191521487905,0.866287777556541,0.14722108,"DNA replication"),
c("GO:0006790","sulfur compound metabolic process",1.84189977492995,1.44963473706621,0.860572170185384,0.143392,"DNA replication"),
c("GO:0019637","organophosphate metabolic process",2.93050388131,1.40258441052798,0.801983608038446,0.1529118,"DNA replication"),
c("GO:0042133","neurotransmitter metabolic process",3.12829656337943,2.38995320520681,0.824737696342239,0.12638798,"DNA replication"),
c("GO:0044281","small molecule metabolic process",7.84989205824262,2.26019183433259,0.88222310898327,0.1348386,"DNA replication"),
c("GO:0072527","pyrimidine-containing compound metabolic process",0.298562307656975,1.88845786342838,0.812869875298591,0.32250485,"DNA replication"),
c("GO:0072528","pyrimidine-containing compound biosynthetic process",0.215883514767351,1.60785554978551,0.818021608654045,0.31311476,"DNA replication"),
c("GO:0016192","vesicle-mediated transport",2.42065132515732,3.06200405558916,0.618423216906132,0,"vesicle-mediated transport"),
c("GO:0015031","protein transport",3.74810527766295,1.97311072005111,0.495214438552703,0.52578184,"vesicle-mediated transport"),
c("GO:0015833","peptide transport",0.234256579853934,1.91926207730473,0.668566427494455,0.59459352,"vesicle-mediated transport"),
c("GO:0033036","macromolecule localization",5.54407238987644,1.82261280306103,0.58930544491423,0.58253109,"vesicle-mediated transport"),
c("GO:0042886","amide transport",0.330715171558495,1.91926207730473,0.658223408320978,0.61738691,"vesicle-mediated transport"),
c("GO:0051641","cellular localization",4.34982315924854,1.65792912426216,0.57108363001552,0.54700372,"vesicle-mediated transport"),
c("GO:0071705","nitrogen compound transport",5.32359560883744,1.76926085830689,0.581788162601663,0.60273194,"vesicle-mediated transport"));

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

