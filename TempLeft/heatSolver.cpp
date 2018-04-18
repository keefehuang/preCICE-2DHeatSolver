#ifndef HSOLVER_H_
#define HSOLVER_H_


#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include "configFile.cpp"

#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;

void printArray(std::vector<double> temp, int N, int solverNum){

  int row, col;

    for(row = 0; row < N+solverNum-1; row++){
      std::cout << "\n" << std::endl;
      for(col = 0; col < N+2; col++){
        std::cout << temp[row*(N+2)+col] << " ";
      }
    }

    std::cout << "\n"; 

}

void heatSolver(std::vector<double> &temp, 
                    std::vector<double> &temp_n, 
                    std::vector<double> &nearBound, 
                    std::vector<double> &farBound, 
                    int N, 
                    double dt,
                    int solverNum){
  int row, col;
  std::vector<double> nearBound_n(nearBound);
  std::vector<double> last_temp(temp), last_bound(nearBound);
  
  double offDiag = -dt*(N*N), diag = 1+4*dt*(N*N), precision = 1e-3, error = 10;

  while(error > precision){
  	std::copy(temp_n.begin(), temp_n.end(), last_temp.begin());
    std::copy(nearBound_n.begin(), nearBound_n.end(), last_bound.begin());
    for(row = 1; row < N+solverNum-1; row++){
      for(col = 1; col < N+1; col++){
        temp_n[row*(N+2)+col] = (temp[row*(N+2)+col] - offDiag*temp_n[(row+1)*(N+2)+col] 
          - offDiag*temp_n[(row-1)*(N+2)+col] - offDiag*temp_n[row*(N+2)+col+1] - offDiag*temp_n[row*(N+2)+col-1])/diag;
      }
    }

    for(col = 1; col < N+1; col++){
      nearBound_n[col] = (nearBound[col] - offDiag*temp_n[(N+solverNum-1)*(N+2)+col] - offDiag*nearBound_n[col-1] 
        - offDiag*nearBound_n[col+1] - offDiag*farBound[col])/diag;
      temp_n[(N+solverNum-1)*(N+2)+col] = (temp[(N+solverNum-1)*(N+2)+col] - offDiag*temp_n[(N+solverNum-2)*(N+2)+col] 
        - offDiag*temp_n[(N+solverNum-1)*(N+2)+col+1] - offDiag*temp_n[(N+solverNum-1)*(N+2)+col-1] - offDiag*nearBound_n[col])/diag;
      
    }
    error = 0;
    for(row = 1; row < N+solverNum-1; row++){
      for(col = 1; col < N+1; col++){
        error += pow(temp_n[row*(N+2)+col] - last_temp[row*(N+2)+col], 2);
      }
    }

    for(col = 1; col < N+1; col++){
    	error += pow(nearBound_n[col] - last_bound[col], 2);
    }
    error = sqrt(error);

  }

  nearBound = nearBound_n;

return;

}

void exportScalarData(std::ofstream& outFile, int N_slices, double* data, std::string dataname)
{  
  outFile << "SCALARS " << dataname << " float" << std::endl;
  outFile << "LOOKUP_TABLE default" << std::endl;

  for(int i = 0; i < N_slices; i++)
  { 
    outFile << data[i] << std::endl;      
  }
  
  outFile << std::endl;
}

void exportVectorData(std::ofstream& outFile, int N_slices, double* data, const char* dataname)
{
  outFile << "VECTORS " << dataname << " float" << std::endl;

  for(int i = 0; i < N_slices; i++)
  {   
    double vx = data[i]; 
    double vy = 0.0; 
    double vz = 0.0;           
    outFile << vx << "  " << vy << "  " << vz << std::endl;      
  }  
  
  outFile << std::endl;          
}

void initializeWriting(std::ofstream&filestream)
{
  filestream.setf(std::ios::showpoint);
  filestream.setf(std::ios::scientific);
  filestream << std::setprecision(16);
}

void writeHeader(std::ostream& outFile)
{
  outFile << "# vtk DataFile Version 2.0" << std::endl << std::endl
          << "ASCII" << std::endl << std::endl
          << "DATASET UNSTRUCTURED_GRID" << std::endl << std::endl;
}

void exportMesh(std::ofstream& outFile, int N_slices, double* grid)
{  
  outFile << "POINTS " << N_slices << " float "<<std::endl << std::endl;
  
  for (int i = 0; i<N_slices; i++)
  {
    double x = grid[2 * i + 0]; 
    double y = grid[2 * i + 1];
    double z = 0.0;
    outFile << x << "  " << y << "  " << z << std::endl;
  }
  outFile << std::endl;
}

void writeVtk(double t, int iteration, const char* filename_prefix, int N_slices, double* grid, double* templeft, double* tempright){
  std::stringstream filename_stream;
  filename_stream << filename_prefix <<"_"<< iteration <<".vtk";
  std::string filename = filename_stream.str();
  printf("writing timestep at t=%f to %s\n", t, filename.c_str());

  std::ofstream outstream(filename);  

  initializeWriting(outstream);
  writeHeader(outstream);
  exportMesh(outstream, N_slices, grid);
  
  outstream << "POINT_DATA " << N_slices << std::endl;
  outstream << std::endl;  
    
  exportScalarData(outstream, N_slices, templeft, "templeft");
  exportScalarData(outstream, N_slices, tempright, "tempright");
  
  outstream.close();    
}

// int main(){

// 	std::vector<double> leftNodes_n((N+2)*(N+1), temp_bound), leftNodes((N+2)*(N+1),temp_bound), 
//                       midNodes(N+2,temp_bound), rightBound(N+2,temp_bound), midNodes_n(N+2,temp_bound);

// 	for(int row=1; row <= N; row++){
// 		for(int col=1; col < N+1; col++){
// 			leftNodes[row*(N+2)+col] = temp_init;
// 			}
// 			midNodes[row] = temp_init;
// 			rightBound[row] = temp_init;
// 	}

//   // for(auto item : midNodes){
//   //   std::cout << item << std::endl;
//   // }

//   printArray(leftNodes_n, N, 1);

// 	heatSolver(leftNodes, leftNodes_n, midNodes, rightBound, N, dt, 1);

// 	// for(auto item : leftNodes_n){
// 	// 	std::cout << item << std::endl;
// 	// }
//   printArray(leftNodes_n, N, 1);

// 	heatSolver(leftNodes, leftNodes_n, midNodes, rightBound, N, dt, 1);

//   printArray(leftNodes_n, N, 1);
//   // for(auto item : leftNodes_n){
//   //   std::cout << item << std::endl;
//   // }
// }


#endif