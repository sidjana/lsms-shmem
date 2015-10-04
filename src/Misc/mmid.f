      subroutine mmid(y,dydx,nvar,nr,xs,htot,nstep,yout,rn,ecomp,vr,
     *     vso,nsizm,issym,eunit,df,b,allp1)
c     ===============
c
c------ numerical recipes subroutine mmid
c
      implicit real*8 (a-h,o-z)
c
      parameter (nmax=16)
c
      complex*16 ecomp
      dimension y(nvar),dydx(nvar),yout(nvar),ym(nmax),yn(nmax)
      dimension rn(nr),vr(3,nr)
      real*8 eunit
c
c     =====================================================================
c     dimensioned to stop array bounds problems............................
c     =====================================================================
      integer issym(nsizm)
      real*8  vso(nsizm,nsizm)
c
      external df
c
      h=htot/nstep
      do 11 i=1,nvar
      ym(i)=y(i)
      yn(i)=y(i)+h*dydx(i)
   11 continue
      x=xs+h
      call df(x,yn,yout,nvar,nr,rn,ecomp,vr,vso,nsizm,issym,eunit,
     >        b,allp1)
      h2=h+h
      do 13 n=2,nstep
      do 12 i=1,nvar
      swap=ym(i)+h2*yout(i)
      ym(i)=yn(i)
      yn(i)=swap
   12 continue
      x=x+h
      call df(x,yn,yout,nvar,nr,rn,ecomp,vr,vso,nsizm,issym,eunit,
     >        b,allp1)
   13 continue
      do 14 i=1,nvar
      yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
   14 continue
      return
      end
