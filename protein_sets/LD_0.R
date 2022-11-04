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
revigo.data <- rbind(c("GO:0006091","generation of precursor metabolites and energy",3.17560569937111,1.31289659955724,0.941876241878982,0.08654616,"generation of precursor metabolites and energy"),
c("GO:0016052","carbohydrate catabolic process",1.39338790428508,1.37253333218944,0.951687096770014,0.0646863,"carbohydrate catabolic process"),
c("GO:0005975","carbohydrate metabolic process",5.82619132132875,1.38007645744064,0.944512842134143,0.10295284,"carbohydrate catabolic process"),
c("GO:0043648","dicarboxylic acid metabolic process",1.04559755206008,2.5068514180486,0.853823691765462,0,"dicarboxylic acid metabolic process"),
c("GO:0006108","malate metabolic process",0.0542426908363648,2.2869077692087,0.880968225713351,0.37041205,"dicarboxylic acid metabolic process"),
c("GO:0009150","purine ribonucleotide metabolic process",2.36005048154109,1.32371549437121,0.785720597265359,0.38468667,"dicarboxylic acid metabolic process"),
c("GO:0071705","nitrogen compound transport",4.58568019970729,3.33357069617169,0.579366676865238,0,"nitrogen compound transport"),
c("GO:0015031","protein transport",2.72570819953152,3.67170844833107,0.40215182249732,0.4191796,"nitrogen compound transport"),
c("GO:0015833","peptide transport",0.346155973030943,3.63223301849061,0.612582087305758,0.63460052,"nitrogen compound transport"),
c("GO:0016192","vesicle-mediated transport",1.57215851663744,3.01049040788697,0.625233831394174,0.38729748,"nitrogen compound transport"),
c("GO:0033036","macromolecule localization",3.61581119208329,3.194815274056,0.597942980125443,0.42261954,"nitrogen compound transport"),
c("GO:0042886","amide transport",0.401364748179044,3.63223301849061,0.607714767944658,0.64489699,"nitrogen compound transport"),
c("GO:0051641","cellular localization",2.42254470970673,2.69033390106819,0.60896163378462,0.39773809,"nitrogen compound transport"),
c("GO:0071702","organic substance transport",5.84243123323571,2.81549155862077,0.567686413188827,0.47314391,"nitrogen compound transport"),
c("GO:0072525","pyridine-containing compound biosynthetic process",0.422521648339285,1.32058026018695,0.791347776941325,0.06913715,"pyridine-containing compound biosynthetic process"));

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

