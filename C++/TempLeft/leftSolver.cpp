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
  // Solver name must match participant definition in configuration file.
  std::string solverName = "TEMPLEFT";
  
  // Only one instantiation of SolverInterface required for parallel implementation.
  SolverInterface interface(solverName, 0, 1);
  interface.configure(configFileName);

  int row, col;

  // Vectors used for computation defined and initialized.
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

  // Mesh geometries are defined.
  interface.setMeshVertices(meshID, N + 2, grid, vertexIDs.data());

  // preCICE is initialized: data communication and data structures are setup
  cout << "Initialize preCICE..." << endl;
  interface.initialize();

  // Data is initialized before the start of the simulation
  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(leftDataID, N + 2, vertexIDs.data(), midNodes.data());
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(rightDataID, N + 2, vertexIDs.data(), rightBound.data());
  }
  
  // While loop is used for simulation process. Auxilliary method isCouplingOngoing is used as a break condition
  while (interface.isCouplingOngoing()) {
    // Basic structure here to allow easy implementation of checkpointing here. Checkpointing is inherent in the update
    // step here and does not use the preCICE API
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }
    
    // Arbitrary function implemented here to perform update step on data
    heatSolver(leftNodes, leftNodes_n, midNodes, rightBound, N, dt, 1);
    // Updated data is written to preCICE. Unique VertexIDs are referenced to write data to correct location.
    interface.writeBlockScalarData(leftDataID, N + 2, vertexIDs.data(), midNodes.data());
    
    // Time step is advanced.
    interface.advance(dt);
    interface.readBlockScalarData(rightDataID, N + 2, vertexIDs.data(), rightBound.data());

    // Basic structure here to allow easy implementation of checkpointing here. Checkpointing is inherent in the update
    // step here and does not use the preCICE API
    if (interface.isActionRequired(actionReadIterationCheckpoint())) {
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else{
      // Rudimentary form of checkpointing performed here. Updates to local data are only saved if convergence occurs
      leftNodes = leftNodes_n;
    }
  }

  // Finalize() called to perform cleanup and close communication channels
  interface.finalize();

  delete [] grid;

  return 0;
}