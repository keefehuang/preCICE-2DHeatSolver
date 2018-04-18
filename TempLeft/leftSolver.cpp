#include "precice/SolverInterface.hpp"
#include "configFile.cpp"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "heatSolver.cpp"


using std::cout;
using std::endl;

using namespace precice;
using namespace precice::constants;

int main(int argc, char** argv){

  std::string solverName = "TEMPLEFT";
  
  SolverInterface interface(solverName, 0, 1);
  interface.configure(configFileName);

  int row, col;

  std::vector<double> leftNodes_n((N+2)*N, temp_bound), leftNodes((N+2)*(N+1),temp_bound), 
                      midNodes(N+2,temp_bound), rightBound(N+2,temp_bound), midNodes_n(N+2,temp_bound);

  int dimensions        = interface.getDimensions();
  int meshID            = interface.getMeshID("leftNodes");
  int leftDataID        = interface.getDataID("tempMid", meshID);
  int rightDataID       = interface.getDataID("rightBound", meshID);
  

  std::vector<int> vertexIDs(N+2, 0);
  double* grid = new double[(N+2)*dimensions];  

  for(row=1; row <= N; row++){
    for(col=2; col <= N+1; col++){
      leftNodes[row*(N+2)+col] = temp_init;
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
    interface.writeBlockScalarData(leftDataID, N + 2, vertexIDs.data(), midNodes.data());
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(rightDataID, N + 2, vertexIDs.data(), rightBound.data());
  }
  
  while (interface.isCouplingOngoing()) {
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }
    
    heatSolver(leftNodes, leftNodes_n, midNodes, rightBound, N, dt, 1);
    
    interface.writeBlockScalarData(leftDataID, N + 2, vertexIDs.data(), midNodes.data());
    
    interface.advance(dt);
    
    interface.readBlockScalarData(rightDataID, N + 2, vertexIDs.data(), rightBound.data());

    if (interface.isActionRequired(actionReadIterationCheckpoint())) {
      
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else{
      leftNodes = leftNodes_n;
    }
  }
  interface.finalize();

  delete [] grid;

  return 0;
}