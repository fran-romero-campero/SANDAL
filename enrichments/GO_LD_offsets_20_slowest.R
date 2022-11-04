
revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006418","tRNA aminoacylation for protein translation",0.983965528172663,2,0.603843430405122,0,"tRNA aminoacylation for protein translation"),
c("GO:0006108","malate metabolic process",0.0542426908363648,2,0.840637516084543,0.36855849,"tRNA aminoacylation for protein translation"),
c("GO:0006396","RNA processing",3.95590749649223,2,0.719475428625882,0.5890698,"tRNA aminoacylation for protein translation"),
c("GO:0006508","proteolysis",5.13675685219537,2,0.793282613783876,0.36407788,"tRNA aminoacylation for protein translation"),
c("GO:0007165","signal transduction",7.83941060962102,2,0.984457381767146,0.01585203,"signal transduction"),
c("GO:0008152","metabolic process",60.5596218240454,2,1,0,"metabolic process"),
c("GO:0015031","protein transport",2.72570819953152,2,1,0,"protein transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "",
  inflate.labels = T,
  lowerbound.cex.labels = 0,
  position.legend = "none"
)

