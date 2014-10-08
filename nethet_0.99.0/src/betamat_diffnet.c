#include <R.h>
#include <math.h>


void betamat_diffnet(double *betamat, int *ind1, int *ind2, int *uptrirownr, int *uptricolnr, int *lind1, int *lind2, double *sig1, double *sig2, double *sig, int *k){

  int i,j,colnri,rownri,colnrj,rownrj; int l1=*lind1; int l2=*lind2; int kk=*k;

  for (i=0;i<l1;i++){
    colnri=uptricolnr[ind1[i]]-1; /*8!indexing*/
    rownri=uptrirownr[ind1[i]]-1; 
    for (j=0;j<l2;j++){
      colnrj=uptricolnr[ind2[j]]-1;
      rownrj=uptrirownr[ind2[j]]-1;
      betamat[l1*j+i]=sig[colnri+kk*rownri]*sig[colnrj+kk*rownrj]-sig1[colnri+kk*rownri]*sig[colnrj+kk*rownrj]-sig[colnri+kk*rownri]*sig2[colnrj+kk*rownrj]+sig1[colnri+kk*rownri]*sig2[colnrj+kk*rownrj]+sig[colnri+kk*rownrj]*sig[rownri+kk*colnrj]+sig[colnri+kk*colnrj]*sig[rownri+kk*rownrj];
   
    }
  }
}
