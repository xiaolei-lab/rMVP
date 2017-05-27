MVP.Version <-
function(){
##############################################################################################
# MVP: A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For GWAS
# Designed by Xiaolei Liu
# Writen by Lilin Yin, Zhiwu Zhang, and Xiaolei Liu
# Build date: Aug 30, 2017
# Last update: May 25, 2017
##############################################################################################
version <- "0.1"

cat(paste("#", paste(rep("-", 30), collapse=""), "Welcome to MVP", paste(rep("-", 30), collapse=""), sep=""), "\n")

cat("#", "A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For GWAS","\n")
cat(paste("#   ", "Version: ", version, sep=""), "\n")
cat(paste("#  ", "Authors: Lilin Yin, Zhiwu Zhang, and Xiaolei Liu", "\n"))
cat(paste("#  ", "Contact: xiaoleiliu@mail.hzau.edu.cn", paste(rep(" ", 9), collapse=""), "\n"))
cat(paste(rep("*", 80), collapse=""), "\n")
}#end of MVP.version function
