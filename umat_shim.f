!     umat_shim.f
!
!     CalculiX's external-behaviour loader (external.c) loads the DLL for
!     material "@CURMOD" and calls the symbol "umat_" (gfortran name for
!     a subroutine named "umat").  Our CHILE implementation lives in
!     umat_user (exported as "umat_user_"), so this shim just forwards.
!
      subroutine umat(amat,iel,iint,kode,elconloc,emec,emec0,
     &    beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &    icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &    iorien,pgauss,orab,pnewdt,ipkon)
      implicit none
      character*80 amat
      integer iel,iint,kode,ithermal,icmd,ielas,mi,nstate_
      integer ipkon,iorien
      double precision elconloc(*),emec(6),emec0(6),beta(6)
      double precision xokl(3,3),voj,xkl(3,3),vj
      double precision t1l,dtime,time,ttime,pnewdt
      double precision xstateini(*),xstate(*)
      double precision stre(6),stiff(21)
      double precision pgauss(3),orab(*)
      external umat_user
      call umat_user(amat,iel,iint,kode,elconloc,emec,emec0,
     &    beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &    icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &    iorien,pgauss,orab,pnewdt,ipkon)
      return
      end subroutine
