#ifndef HEADERS_HPP
#define HEADERS_HPP

#include <iostream>
#include <fstream>
#include <Eigen>
#include <mpi.h>
#include <math.h>
#include <array>

using namespace std;

//================================================================================
// SPECIAL TYPES
//================================================================================

// Types for dense storage
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ScaVector;
typedef Eigen::Matrix<int,    Eigen::Dynamic, 1> IntVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ScaMatrix;
typedef Eigen::Matrix<int,    Eigen::Dynamic, Eigen::Dynamic> IntMatrix;

// Type for sparse storage
typedef Eigen::SparseMatrix<double> SpMatrix;

// Structure for mesh
struct Mesh
{
  int nbOfNodes;              // number of nodes
  int nbOfTri;                // number of triangles
  ScaMatrix coords;           // coordinates for each node          (Size: nbOfNodes x 3)
  IntMatrix triNodes;         // nodes for each triangle            (Size: nbOfTri x 3)
  IntVector triNum;           // gmsh number for each triangle      (Size: nbOfTri)
  IntVector triPart;          // partition number for each triangle (Size: nbOfTri)
  
  // Infos for parallel computations
  IntVector numNodesToExch;   // number of nodal values to exchanges between the current proc and each other proc  (Size: nbTasks)
  IntMatrix nodesToExch;      // list of nodal values to exchanges between the current proc and each other proc    (Size: nbTasks x max(numNodesToExch) )

};

// Structure for problem
struct Problem
{
  SpMatrix K;    // stiffness matrix
  SpMatrix M;    // mass matrix
  SpMatrix A;    // system matrix
  ScaVector b;      // right-hand side vector
};

//================================================================================
// FUNCTIONS
//================================================================================

//==== Functions in 'mesh.cpp'

// Read the mesh from a gmsh-file (.msh) and store in a mesh-structure 'mesh'
void readMsh(Mesh& mesh, string fileName);

// Write a field 'vec' in a gmsh-file (.msh)
void exportFieldMsh(ScaVector& vec, Mesh& mesh, string viewName, string fileName);

//==== Functions in 'parallel.cpp'

// Build the local numbering and list of nodes for MPI communications
void buildListsNodesMPI(Mesh& mesh);

// MPI-parallel exchange/add the interface terms
void exchangeAddInterfMPI(ScaVector& vec, Mesh& mesh);

double dotProductMPI(ScaVector& vec1, ScaVector& vec2, IntVector& NodesToRemove);

//==== Functions in 'problem.cpp'

// Compute the matrices of the linear wgsystem
void buildProblem(Problem& p, Mesh& mesh, double alpha, ScaVector& f);

//==== Functions in 'solver.cpp'

// Solution of the system Au=b with Jacobi
void jacobi(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh, double tol, int maxit);

void ConjGrad(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh, double tol, int maxit);

#endif /* HEADERS_HPP */
