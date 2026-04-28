!=======================================================================
! dflux.f -- CalculiX user volumetric heat-flux subroutine
!
! Applies the cure exotherm body source for the CASH2D model:
!     q_v = (1 - Vf) * pr(alpha) * H * dalpha/dt   [W/m^3]
!
! The source is computed in umat_user.f and cached in the cure_state
! module. This routine simply looks it up by (integration point,
! element) and returns it as flux(1).
!
! Activate in the input deck with, e.g.:
!     *DFLUX
!     ELSET_COMPOSITE, BFNU, 1.0
! The 1.0 magnitude column is required syntax; the actual value comes
! from this routine and is not scaled (iscale=0).
!
! Signature: matches the CalculiX 2.20 stub. Verify against the
! prototype in your CalculiX source tree -- argument lists drift
! slightly between minor versions.
!=======================================================================
      subroutine dflux(flux,sol,kstep,kinc,time,noel,npt,coords,
     &  jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
     &  ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale)
!
      use cure_state
      implicit none
!
      character*8  lakonl
      character*20 loadtype
!
      integer kstep,kinc,noel,npt,jltyp,iscale,nmpc,
     &        konl(20),ipompc(*),nodempc(3,*),ikmpc(*),ilmpc(*)
      real*8  flux(2),sol,time(2),coords(3),temp,press,area,
     &        vold(0:4,*),co(3,*),coefmpc(*)
!
!     Default: no source, no time-fraction scaling.
      flux(1) = 0.0d0
      iscale  = 0
!
!     Only respond to body (volumetric) loads. CalculiX uses 'BFNU'
!     for body-flux non-uniform; 'BF' is the uniform variant. The
!     leading characters are what matters.
      if (loadtype(1:2).ne.'BF') return
!
!     Out-of-range guard: cure_state is sized at compile time
!     (max_el, max_ip in the module). If the mesh exceeds these,
!     bump the parameters and recompile.
      if (noel.le.0 .or. noel.gt.max_el) return
      if (npt .le.0 .or. npt .gt.max_ip) return
!
!     Read the source assembled by umat_user. This is zero before the
!     first UMAT call for this (element, IP), which is the safe
!     start-of-analysis default.
      flux(1) = q_exo_arr(npt, noel)
!
      return
      end
