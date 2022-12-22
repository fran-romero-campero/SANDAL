
revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006260","DNA replication",1.47064002254612,0.759068420232301,0.06698262,"DNA replication"),
c("GO:0005975","carbohydrate metabolic process",5.82619132132875,0.874038714002005,0.10373946,"DNA replication"),
c("GO:0006414","translational elongation",0.456556209988225,0.742206423635257,0.30635107,"DNA replication"),
c("GO:0006470","protein dephosphorylation",0.442878672241423,0.742631602766077,0.23615188,"DNA replication"),
c("GO:0006810","transport",18.081452126523,1,0,"transport"),
c("GO:0015979","photosynthesis",0.217587118211011,0.909629373401142,0,"photosynthesis"));

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
  inflate.labels = T,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,
  position.legend = "none"
)

