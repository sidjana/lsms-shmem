

/*     subroutine plglmax(lmax,x,plm)
c     ****************************************************************
c     Associated Legendre function....................................
c     Calclates all the p(l,m)'s up to lmax...........................
c     based on the formulae given in "Numerical Recipes" pages 180-183
c     (Equations 6.6.7, 6.6.8 and 6.6.9...............................
c     W. H. Press, B. P. Flannery, S A Teukolsky and W. T. Vetterling.
c     Cambridge Univ Press 1986.......................................
c     ****************************************************************
c
c Inputs:   lmax    integer scalar, max l for Plm
c           x       real*8 scalar, argument of Plm(x)
c Returns:  plm     real*8 array of ((lmax+1)*(lmax+2)/2), Plm(x)
c
*/

__device__
__inline__ void plglmax_device(int lmax, double x, double *plm)
{

  double zero=0;
  double one=1;
  double two=2;
  double tol=1e-12;

  if(one-fabs(x) <= tol) {
    int lend=(lmax+1)*(lmax+2)/2;
    for(int l=threadIdx.x;l<lend;l+=blockDim.x) //Parallel
    {
      plm[l]=0;
    }
    if(x < zero) {
      for(int l=threadIdx.x;l<=lmax;l+=blockDim.x) //Parallel
      {
        int i=(l+1)*l/2; 
        plm[i]=one-two*(l%2);
      }
    } else {
      for(int l=threadIdx.x;l<=lmax;l+=blockDim.x) //Parallel
      {
        int i=(l+1)*l/2; 
        plm[i]=one;
      }
    }
    return;
  }
  if(threadIdx.x!=0) return; //do rest in serial
  //begin calculation of p(l,m)'s...................................
  if(lmax == zero) {
    //special case lmax=0..........................................
    plm[0]=one;
  } else {

    //        =============================================================
    //       minus sign added to be consistant with Numerical Recipes
    //       which has (-1)^m factor in plm :  July 97  by xgz.............
    //        =============================================================
    double somx2=-sqrt((one-x)*(one+x));
    if(lmax == one ) {
      //         ==========================================================
      //         special case lmax=1.......................................
      plm[0]=one;
      plm[1]=x;
      plm[2]=somx2;
    } else {
      for(int m=0;m<=lmax;m++) {

        //              =======================================================
        //                                       m       m
        //              calculate the first two P   and P
        //                                       m       m+1
        //              =======================================================

        if(m==zero) {
          plm[0]=one;
          plm[1]=x;
        } else {
          double pmm=somx2;
          double fact=one;
          for(int i=2;i<=m;i++) {
            fact=fact+two;
            pmm=pmm*fact*somx2;
          }
          int mm=(m+1)*(m+2)/2;
          plm[mm-1]=pmm;

          if( mm+m+1 <= (lmax+1)*(lmax+2)/2 ) {
            plm[mm+m]=x*(2*m+1)*pmm;
          }

        } //end else
        //              =======================================================
        //                                 m        m
        //              calculate the rest P     to P
        //                                  m+2      lmax
        //              =======================================================
        int ll=(m+2)*(m+1)/2;
        double fact=(two*m+one)*x;
        for(int l=m+2; l<=lmax; l++) {
          double pmm=(l+m-1)*plm[ll-1];
          fact=fact+two*x;
          ll=ll+l-1;
          plm[ll+l-1]=( fact*plm[ll-1] - pmm )/double(l-m);
        }

      } //end for m
    } //end else
  }
}
