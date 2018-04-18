#include "precice/SolverInterface.hpp"
#include "../TempLeft/configFile.cpp"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "../TempLeft/heatSolver.cpp"

using std::cout;
using std::endl;

using namespace precice;
using namespace precice::constants;

int main(int argc, char** argv){

  std::string solverName = "TEMPRIGHT";

  SolverInterface interface(solverName, 0, 1);
  interface.configure(configFileName);

  cout << "Configure preCICE..." << endl;

  
  int row, col;

  std::vector<double> rightNodes_n((N+2)*N, temp_bound), rightNodes((N+2)*N,temp_bound), 
                      midNodes(N+2,temp_bound), rightBound(N+2,temp_bound), rightBound_n(N+2,temp_bound);

  int dimensions        = interface.getDimensions();
  int meshID            = interface.getMeshID("rightNodes");
  int leftDataID 		    = interface.getDataID("tempMid", meshID);
  int rightDataID       = interface.getDataID("rightBound", meshID);
  

  std::vector<int> vertexIDs(N+2);
  double* grid = new double[(N+2)*dimensions];

  for(row=1; row < N; row++){
  	for(col=1; col < N+1; col++){
  		rightNodes[row*(N+2)+col] = temp_init;
  	}
  	midNodes[row] = temp_init;
  	rightBound[row] = temp_init;
  }

  for (int i = 0; i < N+2; i++) {
    for (int dim = 0; dim < dimensions; dim++)
      grid[i * dimensions + dim] = i * (1 - dim);
        }

  interface.setMeshVertices(meshID, N + 2, grid, vertexIDs.data());

  cout << "Initialize preCICE..." << endl;
  interface.initialize();

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(rightDataID, N + 2, vertexIDs.data(), rightBound.data());
    interface.fulfilledAction(actionWriteInitialData());
  }
  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(leftDataID, N + 2, vertexIDs.data(), midNodes.data());
  }
  
  while (interface.isCouplingOngoing()) {
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    heatSolver(rightNodes, rightNodes_n, rightBound, midNodes, N, dt, 0);

    interface.writeBlockScalarData(rightDataID, N + 2, vertexIDs.data(), rightBound.data());
    
    interface.advance(dt);
    
    interface.readBlockScalarData(leftDataID, N + 2, vertexIDs.data(), midNodes.data());
    
    if (interface.isActionRequired(actionReadIterationCheckpoint())) { 
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else{
      rightNodes = rightNodes_n;
    }
  }

  interface.finalize();

  delete [] grid;

  return 0;
}