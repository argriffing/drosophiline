.PHONY: all
all : read-coverage rabc-to-likelihood

read-coverage : read-coverage.o xgstats.o phylo.o dgrpstats.o
	g++ xgstats.o dgrpstats.o phylo.o read-coverage.o -o read-coverage
read-coverage.o : read-coverage.C
	g++ -c read-coverage.C

rabc-to-likelihood : rabc-to-likelihood.o xgstats.o phylo.o dgrpstats.o
	g++ xgstats.o dgrpstats.o phylo.o rabc-to-likelihood.o -o rabc-to-likelihood
rabc-to-likelihood.o : rabc-to-likelihood.C
	g++ -c rabc-to-likelihood.C


xgstats.o : xgstats.C xgstats.H
	g++ -c xgstats.C

dgrpstats.o : dgrpstats.C dgrpstats.H
	g++ -c dgrpstats.C

phylo.o : phylo.C phylo.H
	g++ -c phylo.C


clean :
	rm *.o
