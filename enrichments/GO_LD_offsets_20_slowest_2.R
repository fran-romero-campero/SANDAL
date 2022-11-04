revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006108","malate metabolic process",0.0542426908363648,0.791560038485553,0.00867955,"malate metabolic process"),
c("GO:0006511","ubiquitin-dependent protein catabolic process",0.745391180522826,0.755988311547229,0.05582463,"tRNA aminoacylation for protein translation"),
c("GO:0006396","RNA processing",3.95590749649223,0.740937799478007,0.18106872,"tRNA aminoacylation for protein translation"),
c("GO:0006418","tRNA aminoacylation for protein translation",0.983965528172663,0.555628162659437,0.5890698,"tRNA aminoacylation for protein translation"),
c("GO:0006508","proteolysis",5.13675685219537,0.733906425183149,0.35122827,"tRNA aminoacylation for protein translation"),
c("GO:0006520","cellular amino acid metabolic process",5.49669770567442,0.657517978628289,0.59483906,"tRNA aminoacylation for protein translation"),
c("GO:0006979","response to oxidative stress",0.552298974215833,0.90695898912996,0,"response to oxidative stress"),
c("GO:0007165","signal transduction",7.83941060962102,0.886561995347656,0.4289045,"response to oxidative stress"),
c("GO:0007018","microtubule-based movement",0.38076880019981,0.984284864120268,0,"microtubule-based movement"),
c("GO:0015031","protein transport",2.72570819953152,1,0,"protein transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

treemap(
  stuff,
  index = c("representative","description"),
  vSize = "uniqueness",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = T,
  lowerbound.cex.labels = 0,
  position.legend = "none"
)

