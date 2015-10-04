c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine genpot(mynod,num_atoms,n_spin_pola,n_spin_cant,mtasa,
     >                  ztotss,qtotws,mtotws,evec,afm,
     >                  rho,vrold,vrnew,vrms,qrms,
     >                  vx,enxc,vdif,i_vdif,
     >                  madmat,
     >                  omegint,omegmt,omegws,
     >                  u0_const,v0_const,alpha_mad,
     >                  r_mesh,jmt,rins,rmt,vshift,
     >                  emad,emadp,u0,madterm,dq_site,mm_site,qsub,
     &                  ntype,n_per_type,iexch,
     >                  pi4,iprint,istop)
c     ================================================================
c     rmt represents last r of spherical potential
      use MPPModule, only : GlobalSum, GlobalCollect
      use MPPModule, only : NumPEs
c
      implicit   none
c
      include    'atom_param.h'
c
      character  sname*32
      character  istop*32
c                        
      integer    mynod
      integer    num_atoms
      integer    n_spin_pola
      integer    n_spin_cant
      integer    mtasa
      integer    afm
      integer    i_vdif
      integer    iexch
      integer    i
      integer    iprint
      integer    jmt
      integer    is
      integer    isold
      integer    ntype,n_per_type
      integer    m,m1,lm,lm1
c
      real*8     alpgga
      real*8     epcorr
      real*8     epcorrave
      real*8     qtotmt
      real*8     mtotmt
      real*8     ztotss
      real*8     evec(3)
      real*8     rho(iprpts,n_spin_pola),
      real*8     vrold(iprpts,n_spin_pola)
      real*8     vrnew(iprpts,n_spin_pola)
      real*8     vrms(n_spin_pola)
      real*8     qrms(n_spin_pola)
      real*8     vx(iprpts,n_spin_pola)
      real*8     enxc(iprpts)
      real*8     vdif
      real*8     madmat(ntype)
      real*8     omegint
      real*8     omegmt,omegws
      real*8     u0_const
      real*8     v0_const
      real*8     alpha_mad
      real*8     r_mesh(iprpts),rtmp(0:iprpts)
      real*8     rins,rmt,vshift,r_getqm
      real*8     emad(n_spin_pola)
      real*8     emadp(n_spin_pola)
      real*8     u0
      real*8     madterm
c
      real*8     qtotws
      real*8     mtotws
      real*8     qint_local,qint,qint1,qint1_total
      real*8     mint
      real*8     dq_site
      real*8     mm_site
      real*8     qsub(ntype)
      real*8     vmt1
      real*8     vmt
      real*8     vmt0
      real*8     ro3
      real*8     dz
      real*8     sp
      real*8     vxout(ipspin)
      real*8     rhojmt,rhojmt1,rhojmt2,drhot
      real*8     pi4
      real*8     sfac
      real*8     excout
      real*8     zero
      real*8     half
      real*8     one
      real*8     two
      real*8     sqr2
      real*8     three
      real*8     third
      real*8     rhoint
      real*8     work(4)
      real*8     mint_vec(3)
      real*8     msgbuf(4)
c
      parameter  (zero=0.0d+00)
      parameter  (half=0.5d+00)
      parameter  (one=1.0d+00)
      parameter  (two=2.0d+00)
      parameter  (three=3.0d+00)
      parameter  (third=one/three)
      parameter  (sname='genpot')
c
      if(mtasa.eq.0)r_getqm=rins
      if(mtasa.eq.1)r_getqm=rmt
      if(mtasa.eq.2)r_getqm=rmt
      qint1_total=zero
      qint1=zero
      call zeroout(emad,n_spin_pola)
      call zeroout(emadp,n_spin_pola)
      if(n_spin_pola.eq.1) then
         sfac=one
      else
         sfac=half
      endif
      sqr2=sqrt(two)
      rtmp(0)=0.d0
      do i=1,jmt+2
	rtmp(i)=sqrt(r_mesh(i))
      enddo
c
c     ================================================================
c     calculate qtotmt and mtotmt.....................................
c     ----------------------------------------------------------------
      if(rho(jmt+1,1).eq.0.d0)rho(jmt+1,1)=rho(jmt,1)
      if(rho(jmt+2,1).eq.0.d0)rho(jmt+2,1)=rho(jmt,1)
      if(rho(jmt+1,n_spin_pola).eq.0.d0)rho(jmt+1,n_spin_pola)
     >=rho(jmt,n_spin_pola)
      if(rho(jmt+2,n_spin_pola).eq.0.d0)rho(jmt+2,n_spin_pola)
     >=rho(jmt,n_spin_pola)
      if(iprint.ge.1)write(*,*)rho(jmt,1),rho(jmt+1,1),rho(jmt+2,1)
      if(iprint.ge.1)write(*,*)
     >rho(jmt,n_spin_pola),rho(jmt+1,n_spin_pola),rho(jmt+2,n_spin_pola)
ctest call getqm_mt(n_spin_pola,jmt+2,rins,rtmp,rho,iprpts,
      call getqm_mt(n_spin_pola,jmt  ,rins,rtmp,rho,iprpts,
     >              mtasa,qtotmt,mtotmt,r_getqm,iprint)
c     ----------------------------------------------------------------
c
      if(iprint.ge.0) then
         write(6,'(/,'' GENPOT::'',
     >               '' Total charge and moment in W-S cell:'')')
         write(6,'(10x,''qtotws'',t40,''='',f18.11)') qtotws
         write(6,'(10x,''mtotws'',t40,''='',f18.11)') mtotws
      endif
c
c     ================================================================
c     calculate interstitial charge, moment & charge density..........
c     ================================================================
      qint_local = qtotws - qtotmt
      if(mtasa.eq.1) then                    ! ASA case
         msgbuf(1)=qint_local*n_per_type
         msgbuf(2)=omegmt*n_per_type
c        -------------------------------------------------------------
         call GlobalSum(msgbuf,2)
c        -------------------------------------------------------------
	 rhoint=msgbuf(1)/msgbuf(2)
	 if(abs(rhoint).gt.1.d-10) then
	   rhoint=pi4*rhoint/n_spin_pola
	   do is=1,n_spin_pola
	   do i=1,jmt
	     rho(i,is)=rho(i,is)+rhoint*r_mesh(i)*r_mesh(i)
	   enddo
	   enddo
c     ----------------------------------------------------------------
      if(rho(jmt+1,1).eq.0.d0)rho(jmt+1,1)=rho(jmt,1)
      if(rho(jmt+2,1).eq.0.d0)rho(jmt+2,1)=rho(jmt,1)
      if(rho(jmt+1,n_spin_pola).eq.0.d0)rho(jmt+1,n_spin_pola)
     >=rho(jmt,n_spin_pola)
      if(rho(jmt+2,n_spin_pola).eq.0.d0)rho(jmt+2,n_spin_pola)
     >=rho(jmt,n_spin_pola)
ctest call getqm_mt(n_spin_pola,jmt+2,rins,rtmp,rho,iprpts,
      call getqm_mt(n_spin_pola,jmt  ,rins,rtmp,rho,iprpts,
     >              mtasa,qtotmt,mtotmt,r_getqm,iprint)
c     ----------------------------------------------------------------
	 endif
         qint=zero
         mint_vec(1)=zero
         mint_vec(2)=zero
         mint_vec(3)=zero
         rhoint=zero
      else                 ! Muffin-tin or ASA-MT cases
ctest    mint = (mtotws - mtotmt)*i_vdif
         mint=0.d0
         if(i_vdif.ne.0)mint = (mtotws - mtotmt)
         mint_vec(1)=mint*evec(1)
         mint_vec(2)=mint*evec(2)
         mint_vec(3)=mint*evec(3)
         msgbuf(1)=qint_local*n_per_type
         msgbuf(2)=mint_vec(1)*n_per_type
         msgbuf(3)=mint_vec(2)*n_per_type
         msgbuf(4)=mint_vec(3)*n_per_type
c        -------------------------------------------------------------
         call GlobalSum(msgbuf,4)
c        -------------------------------------------------------------
         qint=msgbuf(1)
         if(qint.le.zero) then
            write(6,'('' GENPOT:: negative qint:'',1d20.13)')qint
            call kill(0,9)
            call fstop(sname)
         endif
         mint_vec(1)=msgbuf(2)
         mint_vec(2)=msgbuf(3)
         mint_vec(3)=msgbuf(4)
         rhoint=msgbuf(1)/omegint
         if(iprint.ge.0) then
            write(6,'(10x,''qint'',t40,''='',f18.11)') qint
            write(6,'(10x,''mint'',t40,''='',f18.11)') mint
         endif
      endif
      if(n_spin_cant.eq.2) then
         mint=sqrt(mint_vec(1)*mint_vec(1)+mint_vec(2)*mint_vec(2)+
     >             mint_vec(3)*mint_vec(3))
      else
         mint=mint_vec(3)
      endif
      if(abs(mint).lt.1.0d-5) then
         mint=zero
         mint_vec(1)=zero
         mint_vec(2)=zero
         mint_vec(3)=zero
      endif
c
      if(mtasa.ge.2) then
      if(rho(jmt+1,1).eq.0.d0)rho(jmt+1,1)=rho(jmt,1)
      if(rho(jmt+2,1).eq.0.d0)rho(jmt+2,1)=rho(jmt,1)
      if(rho(jmt+1,n_spin_pola).eq.0.d0)rho(jmt+1,n_spin_pola)
     >=rho(jmt,n_spin_pola)
      if(rho(jmt+2,n_spin_pola).eq.0.d0)rho(jmt+2,n_spin_pola)
     >=rho(jmt,n_spin_pola)
ctest call getqm_mt(n_spin_pola,jmt+2,rins,rtmp,rho,iprpts,
      call getqm_mt(n_spin_pola,jmt  ,rins,rtmp,rho,iprpts,
     >              1,qtotmt,mtotmt,r_getqm,iprint)
         qint1 = qtotws - qtotmt+rhoint*(omegmt-omegws)
c     ================================================================
c     calculate qsub on this node.....................................
c     ================================================================
c   Make sure that the charge within each WS cell remain the same
c   as qtotws when calculating using qtotmt+rhoint*vol_int
      qsub(mynod+1)=ztotss-qtotws+rhoint*omegws   !(23)
      dq_site=qtotws-ztotss
      work(1)=mtotws*evec(1)
      work(2)=mtotws*evec(2)
      work(3)=mtotws*evec(3)
      else
      qsub(mynod+1)=ztotss-qtotmt+rhoint*omegmt
c
c     ================================================================
c     calculate the no. of excess electrons and the moment on the site
c     ================================================================
      dq_site=qtotmt-ztotss+qint/dble(num_atoms)
      work(1)=mtotmt*evec(1)+mint_vec(1)/dble(num_atoms)
      work(2)=mtotmt*evec(2)+mint_vec(2)/dble(num_atoms)
      work(3)=mtotmt*evec(3)+mint_vec(3)/dble(num_atoms)
      endif
      mm_site=sqrt(work(1)*work(1)+work(2)*work(2)+work(3)*work(3))
      if(mtotws.lt.zero) then
         mm_site=-mm_site
      endif
c
c     ================================================================
c     obtain the qsub from all other nodes............................
c     ----------------------------------------------------------------
      call GlobalCollect(qsub)
c     ----------------------------------------------------------------
c
c     ================================================================
c     calculate muffin-tin zero potential and its contribution to the
c     Coulomb energy..................................................
c     note: lattice constant factor is included in madmat.............
c     ================================================================
      call getvmt(ntype,mynod+1,mtasa,qsub,madmat,rmt,rhoint,
     &     qtotmt-ztotss,qint1,vshift,vmt,vmt1,u0)
      madterm=-(vmt1-alpha_mad*rhoint)
c        =============================================================
c        calculate vmt, the electro-static contribution to the muffin-
c        tin zero potential through a global sum......................
c        vmt1 is the site-dependent constant potential................
c        =============================================================
      if(mtasa.eq.0) then                ! Muffin-tin case
         msgbuf(1)=vmt*n_per_type
         msgbuf(2)=u0*n_per_type
c        -------------------------------------------------------------
         call GlobalSum(msgbuf,2)
c        -------------------------------------------------------------
         vmt=msgbuf(1)/omegint
         u0=msgbuf(2)
c        =============================================================
c        calculate the exchange-correlation potential related parameters
c        =============================================================
         ro3=( ( pi4/three)*rhoint )**(-third)
         dz=mint/qint
         if(dz.lt.-1.d0)write(*,*)'genpot,dz,mint,qint',dz,mint,qint
         rhojmt=pi4*rmt**2*rhoint
      else if(mtasa.ge.2) then                ! ASA-MT case
         msgbuf(1)=vmt*n_per_type
         msgbuf(2)=u0*n_per_type
         msgbuf(3)=qint1*n_per_type
         msgbuf(4)=qint1*(vmt1+two*(ztotss-qtotmt)/rmt)*n_per_type
c        -------------------------------------------------------------
         call GlobalSum(msgbuf,4)
c        -------------------------------------------------------------
         vmt=msgbuf(1)/omegint
         u0=msgbuf(2)
	 qint1_total=msgbuf(3)
	 vmt1=vmt1-msgbuf(4)/qint1_total

c        =============================================================
c        calculate the exchange-correlation potential related parameters
c        =============================================================
         call interp(rtmp(1),rho(1,1),jmt,
     >   sqrt(rmt),rhojmt1,drhot,.false.)
         call interp(rtmp(1),rho(1,n_spin_pola),jmt,
     >   sqrt(rmt),rhojmt2,drhot,.false.)
         rhojmt=rhojmt1+(n_spin_pola-1)*rhojmt2
         ro3=1.d10
         if(rhojmt.gt.1.d-10)
     >   ro3=(three*rmt*rmt/rhojmt)**third
         if(i_vdif.ne.0)dz=(rhojmt1-rhojmt2)/rhojmt
C        ro3=((pi4/three)*rhoint)**(-third)
C        dz=mint/qint
      else                               ! ASA case
         msgbuf(1)=vmt*n_per_type
         msgbuf(2)=u0*n_per_type
c        -------------------------------------------------------------
         call GlobalSum(msgbuf,2)
c        -------------------------------------------------------------
         vmt=msgbuf(1)/dble(num_atoms)
         u0=msgbuf(2)
c        =============================================================
c        calculate the exchange-correlation potential related parameters
c        =============================================================
         call interp(rtmp(1),rho(1,1),jmt,
     >   sqrt(rmt),rhojmt1,drhot,.false.)
         call interp(rtmp(1),rho(1,n_spin_pola),jmt,
     >   sqrt(rmt),rhojmt2,drhot,.false.)
         rhojmt=rhojmt1+(n_spin_pola-1)*rhojmt2
         ro3=1.d10
         if(rhojmt.gt.1.d-10)
     >   ro3=(three*rmt*rmt/rhojmt)**third
         if(i_vdif.ne.0)dz =(rho(jmt,1)-rho(jmt,n_spin_pola))/rhojmt
      endif
c
      if(mtasa.ge.1) vmt0=zero

      do is=1,n_spin_pola
c        =============================================================
c        calculate vxout, the exchange-correlation potential corres-
c        ponding to the interstial constant charge density, and excout,
c        the exchange-correlation energy.............................. 
c        vmtz is the muffin-tin zero potential.
c        emad will be used in the total energy calculation.
c        emadp will be used in the pressure calculation.
c        =============================================================
	 sp=3.0d0-2.0d0*is
c        -------------------------------------------------------------
         call newexchg(n_spin_pola,sp,
     >               rho(1,1),rho(1,n_spin_pola),
     >               vx(1,is),enxc,vxout(is),excout,ro3,dz,
     >               r_mesh,jmt,iexch)
c        -------------------------------------------------------------
	 if(mtasa.lt.1) then
         emad(is) = sfac*(qint+sp*mint)*excout
         emadp(is)=-sfac*(qint+sp*mint)*three*(excout-vxout(is))
	 endif
      enddo
c
      if(mtasa.ge.1) then
         do is=1,n_spin_pola
	   if(mtasa.eq.1) then
            msgbuf(is)=vxout(is)*n_per_type
	   else
            msgbuf(is)=vxout(is)*qint1*n_per_type
	   endif
         enddo
c        -------------------------------------------------------------
c        if(abs(qint1).le.1.d-10)msgbuf(1)=0.d0
c        if(abs(qint1).le.1.d-10)msgbuf(2)=0.d0
         call GlobalSum(msgbuf,n_spin_pola)
c       write(6,'(''msgbuf'',2d16.8,2i5)')msgbuf(1),msgbuf(2),mynod
c        -------------------------------------------------------------
         do is=1,n_spin_pola
	   if(mtasa.eq.1) then
            vxout(is)=msgbuf(is)/dble(num_atoms)
	   else
            vxout(is)=msgbuf(is)/qint1_total
            emad(is)=sfac*vxout(is)*qint1_total
	   endif
         enddo
      endif
c
      call zeroout(vrnew,iprpts*n_spin_pola)
      if(iexch.ge.100)then
        epcorrave=n_per_type*
     >  epcorr(rhojmt,0.d0,r_mesh(jmt),0,alpgga,50.d0)
        call GlobalSum(epcorrave)
        epcorrave=epcorrave/dble(num_atoms)
        vxout(1)=epcorrave
        vxout(n_spin_pola)=epcorrave
      endif
      do is=1,n_spin_pola
         if(afm.eq.1) then
            isold=n_spin_pola+1-is
         else
            isold=is
         endif
c        =============================================================
c        generate new potential.......................................
c        -------------------------------------------------------------
         call newpot(n_spin_pola,ztotss,
     >               rho(1,1),rho(1,n_spin_pola),
     >               vrold(1,isold),vrnew(1,is),
     >               vrms(is),vx(1,is),
     >               vmt1,vmt,vxout(is),
     >               rtmp,jmt,rins,rmt,
     >               mtasa,iexch)
c        -------------------------------------------------------------
      enddo
c
c     ================================================================
c     vdif is the difference between the spin-down muffin-tin zero and
c     the spin-up muffin-tin zero potentials..........................
c     Note that vr(ir,1) and vr(ir,2) are all relative to their own
c     muffin-tin zero.................................................
c     ================================================================
      vdif=vxout(n_spin_pola)-vxout(1)
c     shift minority potential by an amount coresponding to -i_vdif Tesla
      if(i_vdif.lt.0.and.n_spin_cant.ne.2) vdif=vdif-i_vdif*4.256e-6
c     ****************************************************************
c     Major printout..................................................
      if(iprint.ge.0) then
         write(6,'(/,10x,''Potential reconstruction:'')')
         write(6,'(10x,''MT sphere volume'',t40,''='',
     >             f18.11)')omegmt
         write(6,'(10x,''Interstitial charge density'',t40,''='',
     >             f18.11)')rhoint
	 if(qint1.ne.zero)
     &   write(6,'(10x,''Charge from shape correction'',t40,''='',
     >             f18.11)')qint1
         write(6,'(10x,''Muffin-tin zero potential'',t40,''='',
     >             3f18.11)')vmt+vxout(1),vmt,vxout(1)
         write(6,'(10x,''Madelung constant shift'',t40,''='',
     >             f18.11)')-alpha_mad*rhoint
         write(6,'(10x,''Coulomb energy due to rho0'')')
         write(6,'(10x,''and the Mad. contribution'',t40,''='',
     >             f18.11)')u0
         write(6,'(10x,''Exc-Energy for rho0'',t40,''='',
     >             f18.11)')excout
         write(6,'(10x,''Point charge'',t40,''='',
     >             f18.11)') qsub(mynod+1)
         write(6,'(10x,''VMT1'',t40,''='',f18.11)') vmt1
         do is=1,n_spin_pola
	    write(6,'(10x,''Spin'',t40,''='',1i5)') is
	    write(6,'(14x,''emad'',t40,''='',f18.11)') emad(is)
	    write(6,'(14x,''emadp'',t40,''='',f18.11)')emadp(is)
	    write(6,'(14x,''vmtz'',t40,''='',3f18.11)') vmt+vxout(is)
     >      ,vmt,vxout(is)
            write(6,'(14x,''Vrms'',t40,''='',f18.11)') vrms(is)
            write(6,'(14x,''Qrms'',t40,''='',f18.11)') qrms(is)
         enddo
      endif
c 
c     ================================================================
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end

