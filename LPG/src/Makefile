OBJECTS = LPG.o ModelInfo.o Point.o ConstraintConsensus.o FeasibilityVector.o TimeKeeper.o MatlabScriptWriter.o ClusterBuilder.o Cluster.o 
CC = g++
CFLAGS = -Wall -g -c -O2 
LFLAGS = ../../solvers/funcadd0.o ../../solvers/amplsolver.a  -lm  


main:$(OBJECTS) 
	$(CC) $(OBJECTS) $(LFLAGS) -o ../bin/LPG

 
	
clean:
	\rm LPG.o ModelInfo.o Point.o ConstraintConsensus.o FeasibilityVector.o TimeKeeper.o MatlabScriptWriter.o ClusterBuilder.o Cluster.o 
	\rm ../tmp/*.m ../tmp/*.data