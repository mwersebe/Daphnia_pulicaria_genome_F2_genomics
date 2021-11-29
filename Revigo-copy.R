# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0007017","\n \n Microtubule-based \n process",0.697,4.426,7.203,5.195,0.993,0.000),
c("GO:0008150","\n \n \n Biological process",100.000,-3.795,7.509,7.352,1.000,0.000),
c("GO:0008152","Metabolic \n process",63.536,0.184,2.762,7.155,1.000,0.000),
c("GO:0009987","Cellular \n process",82.432,-6.752,-1.573,7.268,1.000,0.000),
c("GO:0016042","\n Lipid catabolic \n process",0.532,0.253,5.692,5.078,0.985,0.000),
c("GO:0051179","Localization",18.989,-2.107,-0.670,6.631,1.000,0.000),
c("GO:0051606","Detection \n of stimulus",0.399,-5.871,3.826,4.953,0.886,0.000),
c("GO:0098876","Vesicle-mediated \n transport to \n membrane",0.039,5.878,1.827,3.947,0.859,0.005),
c("GO:0000280","Nuclear division",0.230,3.666,-4.461,4.714,0.804,0.006),
c("GO:0007018","Microtubule-based movement",0.340,0.363,8.7,4.883,0.926,0.006),
c("GO:1903047","Mitotic cell \n cycle process",0.349,-5.648,-3.290,4.895,0.811,0.006),
c("GO:0006928","Movement cell or \n subcellular",1.001,2.489,7.869,5.353,0.993,0.007),
c("GO:0007049","Cell cycle",1.774,6.791,-2.041,5.601,0.992,0.008),
c("GO:0006730","\n \n One-carbon \n metabolic",0.397,-1.814,8.034,4.950,0.980,0.031),
c("GO:0044281","\n \n Small molecule metabolic process",16.281,0.093,-5.779,6.564,0.973,0.047),
c("GO:0071897","\n DNA \n biosynthetic \n process",0.767,1.409,-0.707,5.237,0.977,0.066),
c("GO:0044238","Primary metabolic \n process",50.800,-2.328,-5.401,7.058,0.960,0.133),
c("GO:0046907","Intracellular transport",1.785,6.295,2.912,5.604,0.788,0.236),
c("GO:0071704","Organic substance metabolic process",56.056,-2.673,-4.628,7.101,0.958,0.263),
c("GO:0009581","Detection of external stimulus",0.064,-6.032,2.169,4.160,0.832,0.278),
c("GO:0044237","Cellular metabolic process",56.282,-1.710,-4.744,7.102,0.948,0.287),
c("GO:0009416","Response to light stimulus",0.181,-5.085,3.657,4.609,0.791,0.302),
c("GO:0009628","Response to abiotic stimulus",0.584,-6.194,3.089,5.118,0.883,0.336),
c("GO:0009605","Response to external stimulus",1.716,-5.513,2.948,5.587,0.873,0.389),
c("GO:0006810","Transport",18.119,5.637,3.176,6.610,0.832,0.485),
c("GO:0048285","Organelle fission",0.275,4.054,-4.041,4.791,0.859,0.507),
c("GO:0070925","Organelle assembly",0.627,4.074,-4.420,5.149,0.844,0.553),
c("GO:1990126","Retrograde transport, \n endosome to plasma membrane",0.006,6.448,2.056,3.103,0.832,0.561),
c("GO:0006811","Ion transport",10.217,5.842,3.480,6.361,0.843,0.697));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = uniqueness, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$uniqueness), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 6.5 );
p1 <- p1 + labs (y = "Semantic Space X", x = "Semantic Space Y", size = "Number of GO Terms")+ theme(axis.text = element_text(size = 10))+ theme(axis.title = element_text(size = 10));
p1 <- p1 + theme(legend.key = element_blank()) + theme(legend.text = element_text(size = 10));
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1 + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20)) +theme(legend.title = element_text(size = 20)) +theme(legend.text = element_text(size = 20))


# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

