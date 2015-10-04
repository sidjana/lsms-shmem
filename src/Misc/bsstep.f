      subroutine bsstep(y,dydx,nv,nr,x,htry,s1oeps,yscal,hdid,hnext,
     &     idid,rn,ecomp,vr,vso,nsizm,issym,eunit,df,
     >     b,allp1)
c     =================
c
!      implicit real*8 (a-h,o-z)
      implicit none
      integer nmax,imax,nuse,idid,nsizm
      real*8 one,shrink,grow
      parameter (nmax=16,imax=11,nuse=7,one=1.0d0,shrink=0.95d0,
     & grow=1.2d0)
c
      complex*16 ecomp
      real*8 x,y,dydx,yscal,yerr,ysav,yseq,s1oeps,rn,vr,dysav
      dimension y(nv),dydx(nv),yscal(nv),yerr(nmax),ysav(nmax),
     & dysav(nmax),yseq(nmax),nseq(imax)
      dimension rn(nr),vr(3,nr)
      real*8 eunit,b,allp1,hdid,hnext
      real*8 h,xfrom,xto,den,errmax,xsav,htry,xest
      integer i,ii,nseq,nv,nr
c     ============================================================
c     dimensioned to stop array bounds problems...................
c     ============================================================
      integer   issym(nsizm)
      real*8    vso(nsizm,nsizm)
      real*8 x_extrapol(imax),d_extrapol(nmax,imax)
c
      external df
c
!     save 
      data nseq /2,4,6,8,12,16,24,32,48,64,96/
c
      xfrom=x
      xto=xfrom+htry
      h=htry
      xsav=x
      if(nv.gt.nmax) stop 'bssteps: nv>nmax'
      do i=1,nv
        ysav(i)=y(i)
        dysav(i)=dydx(i)
      enddo
    1 continue
      do i=1,imax
         call mmid(ysav,dysav,nv,nr,xsav,h,nseq(i),yseq,rn,ecomp,vr,
     *        vso,nsizm,issym,eunit,df,b,allp1)
        xest=(h/nseq(i))**2
        call rzextr(i,xest,yseq,y,yerr,nv,nuse,x_extrapol,d_extrapol)
        errmax=0.0d0
        do ii=1,nv,2
          den=abs(yscal(ii))+abs(yscal(ii+1))
          if(den.gt.1.0e-10) then 
            errmax=errmax+(abs(yerr(ii))+abs(yerr(ii+1)))/den
          endif
        enddo
        errmax=errmax*s1oeps
        if (errmax.lt.one) then
          x=x+h
          hdid=h
          idid=i
          call df (x,y,dydx,nv,nr,rn,ecomp,vr,vso,nsizm,issym,eunit,
     >             b,allp1)
          if (i.eq.nuse) then
            hnext=h*shrink
          elseif (i.eq.nuse-1) then
            hnext=h*grow
          else
            hnext=(h*nseq(nuse-1))/nseq(i)
          endif
        return
        endif
      enddo
      h=0.25d0*h/2.d0**(dble(imax-nuse)/2.d0)
      if (x+h.eq.x) then
        write(6,*) xfrom,xto,x,h
        stop 'uflow'
      endif
      goto 1
      end
