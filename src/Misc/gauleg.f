      subroutine gauleg(x1,x2,x,w,n)
c     ================================================================
c     ****************************************************************
c     generates gaussian points and weights to be used in ............
c     Gauss-Legendre quadrature.......................................
c     See:: Numerical Recipes, page 125...............................
c     ****************************************************************
c Input:  x1,x2   real*8 scalars, end points of integration
c         n       integer scalar, number of points
c Returns: x      real*8 array of (n), node points of Gaussian quadrature
c          w      real*8 array of (n), weights associated with x.
c
      implicit   none
c
      integer    n,m
      integer    i,j
c
      real*8     x(n)
      real*8     w(n)
      real*8     x1,x2
      real*8     xm
      real*8     xl
      real*8     z
      real*8     dz
      real*8     p1,p2,p3
      real*8     pp
      real*8     zero
      real*8     fourth
      real*8     half
      real*8     one
      real*8     two
c
      parameter (zero=0.d0)
      parameter (fourth=0.25d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
      real*8     pi
      parameter (pi=3.141592653589793d0)
c
      m=(n+1)/2
      xm=half*(x2+x1)
      xl=half*(x2-x1)
      do i=1,m
        z=cos(pi*(i-fourth)/(n+half))    ! initial guess for the roots
	dz=1.d+10
1       continue                             ! return here for Newton refinement
          p1=one
          p2=zero
          do j=1,n                           ! calculate P(n,z)
            p3=p2
            p2=p1
c recursion relation
            p1=((two*j-one)*z*p2-(j-one)*p3)/j
          enddo
c calculate P'(n,z) scaled by (z^2-1)
          pp=n*(z*p1-p2)
	if(z+dz/64.d0.eq.z) then
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=two*xl*(one-z*z)/(pp*pp)
        w(n+1-i)=w(i)
	else
	  dz=(z*z-one)*p1/pp
          z=z-dz
	  goto 1
	endif
        enddo
      return
      end
