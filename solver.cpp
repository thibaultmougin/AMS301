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
    if((it % 100) == 0){
      if(myRank == 0)
        cout << "   [" << it << "] residual: " << residuNorm << endl;
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
  ScaVector r  = b - Ax; 
  ScaVector p = r;
  ScaVector r0 = r;
 
  r0Norm2 = r0.transpose()*r0;
  int it = 0;

  r_newNorm2 = r0Norm2;
  do
    {
       
      rNorm2 = r.transpose()*r;
     
      Ap = A*p;

      alpha = rNorm2 / (p.transpose()*A*p);

      
      x_new = x + alpha*p;
      
      r_new = r - alpha*Ap;

      r_newNorm2 = r_new.transpose()*r_new;

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