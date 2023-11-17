#include "headers.hpp"

int myRank;
int nbTasks;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);

  array<std::string, 5> texts = {"benchmark/mesh02.msh", "benchmark/mesh01.msh", "benchmark/mesh005.msh", "benchmark/mesh0025.msh", "benchmark/mesh00125.msh",};
    // ^ An array of 3 elements with the type std::string

  for(const auto& text : texts) {   // Range-for!
    
    cout << "USING MESH : "<< text << "\n" <<"\n";
    
    // 1. Initialize MPI
    
  
    // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
    Mesh mesh;
    readMsh(mesh, text);
    buildListsNodesMPI(mesh);
    
    // 3. Build problem (vectors and matrices)
    ScaVector uNum(mesh.nbOfNodes);
    ScaVector uExa(mesh.nbOfNodes);
    ScaVector f(mesh.nbOfNodes);
    for(int i=0; i<mesh.nbOfNodes; ++i){
      double x = mesh.coords(i,0);
      double y = mesh.coords(i,1);
      uNum(i) = 0.;
      uExa(i) = cos(4*M_PI*x)*cos(M_PI*y);
      f(i) = (1+17*M_PI*M_PI)*uExa(i);
    }
    
    Problem pbm;
    double alpha = 1;
    buildProblem(pbm,mesh,alpha,f);
    
    // 4. Solve problem
    double tol = 1e-6; // (Currently useless)
    int maxit = 100000;
    ScaVector test = 0*uNum;

      MPI_Barrier(MPI_COMM_WORLD);
      double starttime = MPI_Wtime();


    //jacobi(pbm.A, pbm.b, uNum, mesh, tol, maxit);

    ConjGrad(pbm.A, pbm.b, uNum, mesh, tol, maxit);

      MPI_Barrier(MPI_COMM_WORLD);
      double endtime = MPI_Wtime();

      if (myRank==0) {
          cout << "   -> time elapsed : " << endtime-starttime << " s "<< endl;
      }
    
    // 5. Compute error and export fields

    ScaVector uErr = uNum - uExa;  
    ScaVector MuErr = pbm.M * uErr;
    ScaVector MuExa = pbm.M * uExa;  
    double L2totLocal = MuErr.dot(uErr);
    double L2exaLocal = MuExa.dot(uExa);
    double L2tot, L2exa;
    MPI_Reduce(&L2totLocal, &L2tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&L2exaLocal, &L2exa, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double L2err = sqrt(L2tot/L2exa);

    if(myRank == 0) printf("\nFinal relative L2 error: %.3e\n",L2err);

    exportFieldMsh(uNum, mesh, "solNum", "benchmark/solNum.msh");
    exportFieldMsh(uExa, mesh, "solRef", "benchmark/solExa.msh");
    exportFieldMsh(uErr, mesh, "solErr", "benchmark/solErr.msh");

    // 6. Finalize MPI
    

  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}