#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Solution of the system Au=b with Jacobi
//================================================================================

void jacobi(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh, double tol, int maxit)
{
  if(myRank == 0)
    cout << "== jacobi" << endl;
  
  // Compute the solver matrices
  int size = A.rows();
  ScaVector Mdiag(size);
  SpMatrix N(size, size);
  for(int k=0; k<A.outerSize(); ++k){
    for(SpMatrix::InnerIterator it(A,k); it; ++it){
      if(it.row() == it.col())
        Mdiag(it.row()) = it.value();
      else
        N.coeffRef(it.row(), it.col()) = -it.value();
    }
  }
  exchangeAddInterfMPI(Mdiag, mesh);

  // Create list of nodes that don't need to be computed in the residu
  IntVector NodesToRemove(mesh.nbOfNodes);
  NodesToRemove = 0*NodesToRemove;
  for (int nTask = 0; nTask<myRank; nTask++){
    int i = 0;
    for (int j = 0; j<mesh.numNodesToExch(nTask); j++){
      i = mesh.nodesToExch(nTask,j);
      if (NodesToRemove(i) == 0) {
        NodesToRemove(i) = 1;
      }
    }
  }

  // Compute the 2-norm of b, for the residu
  double NormB = sqrt(dotProductMPI(b,b,NodesToRemove));

  // Jacobi solver
  double residuNorm = 1e2;
  int it = 0;
  while (residuNorm > tol && it < maxit){
    
    // Compute N*u
    ScaVector Nu = N*u;
    exchangeAddInterfMPI(Nu, mesh);
    
    // Update field
    for(int i=0; i<size; i++){
      u(i) = 1/Mdiag(i) * (Nu(i) + b(i));
    }
    
    // Update residual and iterator
    if((it % 500) == 0){
      
      double residuNormLocal = 0;
      ScaVector Au = A*u;
      exchangeAddInterfMPI(Au, mesh);
      ScaVector residu = Au-b;
      for (int i=0; i<size; i++){
        if (NodesToRemove(i) == 0) {
          residuNormLocal += residu(i)*residu(i);
        } 
      }
      MPI_Allreduce(&residuNormLocal, &residuNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      residuNorm = sqrt(residuNorm)/NormB;

      if(myRank == 0)
        printf("   [%d] residual: %.3e \n",it,residuNorm);
    }
    it++;
  }
  
  if(myRank == 0){
    cout << "   -> final iteration: " << it << " (prescribed max: " << maxit << ")" << endl;
    cout << "   -> final residual: " << residuNorm << " (prescribed tol: " << tol << ")" << endl;
  }
}

void ConjGrad(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh, double tol, int maxit)
{

  double r0Norm2, rNorm2,r_newNorm2, alpha, beta;
  ScaVector Ap, r_new, p_new, x_new;
  
  ScaVector x = u;
  ScaVector Ax = A*x;
  exchangeAddInterfMPI(Ax,mesh);
  ScaVector r  = b - Ax; 
  ScaVector p = r;
  ScaVector r0 = r;
 
  
  int it = 0;

  IntVector NodesToRemove(mesh.nbOfNodes);
  NodesToRemove = 0*NodesToRemove;
  for (int nTask = 0; nTask<myRank; nTask++){
    int i = 0;
    for (int j = 0; j<mesh.numNodesToExch(nTask); j++){
      i = mesh.nodesToExch(nTask,j);
      if (NodesToRemove(i) == 0) {
        NodesToRemove(i) = 1;
      }
    }
  }

  //r0Norm2 = r0.transpose()*r0;
  r0Norm2 = dotProductMPI(r0,r0,NodesToRemove);
  
  r_newNorm2 = r0Norm2;
  do
    {
       
      //rNorm2 = r.transpose()*r;

      rNorm2 = dotProductMPI(r,r,NodesToRemove);
      Ap = A*p;

      exchangeAddInterfMPI(Ap, mesh);
      
      //alpha = rNorm2 / (p.transpose()*A*p);

      alpha = rNorm2 / dotProductMPI(p,Ap,NodesToRemove);


      x_new = x + alpha*p;
      
      r_new = r - alpha*Ap;

      //r_newNorm2 = r_new.transpose()*r_new;

      r_newNorm2 = dotProductMPI(r_new,r_new,NodesToRemove);
      beta = r_newNorm2/rNorm2;

      p_new = r_new + beta*p;
      
      x = x_new; r = r_new; p = p_new;
      it ++;
      
    }
    while (sqrt(r_newNorm2) > (tol * sqrt(r0Norm2)) && it < maxit);

  u = x_new;
  
  if(myRank == 0){
    printf("\r   -> final iteration: %i ( max: %i)\n", it, maxit);
    printf("   -> final residual: %e ( tol: %e)\n", sqrt(r_newNorm2), tol);
  }

}