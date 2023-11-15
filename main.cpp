#include "headers.hpp"

int myRank;
int nbTasks;

int main(int argc, char* argv[])
{
  
  // 1. Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);
 
  // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
  Mesh mesh;
  readMsh(mesh, "benchmark/mesh.msh");
  buildListsNodesMPI(mesh);
  
  // 3. Build problem (vectors and matrices)
  ScaVector uNum(mesh.nbOfNodes);
  ScaVector uExa(mesh.nbOfNodes);
  ScaVector f(mesh.nbOfNodes);
  for(int i=0; i<mesh.nbOfNodes; ++i){
    double x = mesh.coords(i,0);
    double y = mesh.coords(i,1);
    uNum(i) = 0.;
    uExa(i) = x+y;
    f(i) = x+y;
  }


  
  Problem pbm;
  double alpha = 1;
  buildProblem(pbm,mesh,alpha,f);
  
  // 4. Solve problem
  double tol = 1e-9; // (Currently useless)
  int maxit = 1000;
  //jacobi(pbm.A, pbm.b, uNum, mesh, tol, maxit);
  ConjGrad(pbm.A,pbm.b, uNum,mesh,tol,maxit);
  ScaVector uErr = uNum - uExa;
  
  double L2tot = uErr.transpose()*pbm.M * uErr;  

  double L2exa = uExa.transpose()*pbm.M * uExa;  

  double sumerr;
  double sumexa;

  MPI_Allreduce( &L2tot , &sumerr , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
  MPI_Allreduce( &L2exa , &sumexa , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);

  double L2err = sqrt(sumerr/sumexa);

  exportFieldMsh(uNum, mesh, "solNum", "benchmark/solNum.msh");
  exportFieldMsh(uExa, mesh, "solRef", "benchmark/solExa.msh");
  exportFieldMsh(uErr, mesh, "solErr", "benchmark/solErr.msh");
  // 6. Finilize MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}
