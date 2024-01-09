#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>


// Calculates the shielding factor for capm concentric cylindrical shields.
// See eq. (18) of https://arxiv.org/abs/1310.8242

int main (void)
{
  // Define the dimension n of the matrix
  // and the signum s (for LU decomposition)
  int capm = 4; // number of shields.
  int rank = 2*capm; // matrix dimension

  int n=1; // multipole order
  // 1 = uniform field

  int i,j; // indices of matrix
  int m; // index corresponding to which shield
  int s;
  double element;
  double sumrow;
  double summat;
  double sfact;


  // Define geometry
  double r[rank];
  double t[capm];
  double capr1[capm];

  // Define mu
  double mur[capm];

  // Define all the used matrices
  gsl_matrix * a = gsl_matrix_alloc (rank, rank);
  gsl_matrix * inverse = gsl_matrix_alloc (rank, rank);
  gsl_permutation * perm = gsl_permutation_alloc (rank);

  //Define parameters of the problem
  double k=0.025; // scale factor in Fig. 2
  
  capr1[0]=0.5; // R1 = 0.5 m
  for(i=0;i<capm;i++){
    capr1[i]=capr1[0]*pow((1+k),i);
  }

  for(i=0;i<capm;i++){
    t[i]=1./16.*0.0254; // 1/16th of an inch
    mur[i]=40000;
    printf("%d %f %f %f\n",i,capr1[i],t[i],mur[i]);
  }

  // Eq. (17)
  j=0;
  for(i=0;i<capm;i++){
    r[j]=capr1[i];
    printf("r%d %f\n",j,r[j]);
    j++;
    r[j]=capr1[i]+t[i];
    printf("r%d %f\n",j,r[j]);
    j++;
  }

  printf("Fill the Matrix\n");
  
  // Fill the matrix m
  for (i = 0; i < rank; i++){
    for (j = 0; j < rank; j++) {
      if(j<i)element=pow(r[j]/r[i],2*n);
      else if(j>i)element=-1;
      else if(j==i) {
	// Note:  complicated by C starting at 0 index
	if((i+1)%2==1){
	  m=((i+1)+1)/2-1;
	  element=-(mur[m]+1)/(mur[m]-1);
	} else {
	  m=(i+1)/2-1;
	  element=(mur[m]+1)/(mur[m]-1);
	}
	//	printf("Check %d %d %d\n",i,j,m);
      }
      gsl_matrix_set (a, i, j, element);
      printf("%f ",element);
    }
    printf("\n");
  }

  printf("Print the matrix\n");
  // print the matrix
  for (i = 0; i < rank; i++){
    for (j = 0; j < rank; j++) {
      element=gsl_matrix_get (a, i, j);
      printf("%f ",element);
    }
    printf("\n");
  }



  printf("Inverse\n");
  
  // Make LU decomposition of matrix m
  gsl_linalg_LU_decomp (a, perm, &s);
  
  // Invert the matrix m
  gsl_linalg_LU_invert (a, perm, inverse);

  for (i = 0; i < rank; i++){
    for (j = 0; j < rank; j++) {
      element=gsl_matrix_get (inverse, i, j);
      printf("%f ",element);
    }
    printf("\n");
  }

  // Act on vector of 1's
  // i.e., add up along a row.
  // then we add up all the rows anyway.
  // just sum the whole matrix, then.
  summat=0;
  for (i = 0; i < rank; i++){
    sumrow=0;
    for (j = 0; j < rank; j++) {
      element=gsl_matrix_get (inverse, i, j);
      sumrow=sumrow+element;
    }
    printf("%d sumrow %f\n",i,sumrow);
    summat=summat+sumrow;
  }
  printf("summat %f\n",summat);

  sfact=1./(1.+summat);

  printf("sfact %f\n",sfact);


}

