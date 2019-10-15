args<-commandArgs(trailingOnly=TRUE)
sink("sink.txt")

print(as.numeric(args))
