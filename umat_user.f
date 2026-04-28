!=======================================================================
! umat_user.f -- CalculiX user material for thermoset composite curing
! Ported from CURING_THERMOSET_ABAQUS.f (CASH2D model).
!
! Combines the original UMATHT (cure kinetics) and UEXPAN
! (Tg-relative CTE + cure shrinkage) into the mechanical UMAT, since
! CalculiX has no UMATHT/UEXPAN equivalent. The volumetric exotherm
! source must be applied separately in dflux.f, reading the cure
! state through the cure_state module below.
!
! Input deck requirements:
!   *USER MATERIAL,CONSTANTS=40   (no *EXPANSION; eigenstrains handled
!                                  inside this routine)
!   *DEPVAR
!     37
!   *INITIAL CONDITIONS,TYPE=SOLUTION    (set alpha, T_prev, etc.)
!
! NOTE: ncmat_ in CalculiX must be >= 40. If your build was compiled
! with the default lower limit, re-build with -DNCMAT_=40 or edit
! CalculiX.h.
!
! STATE VARIABLE LAYOUT (xstate / xstateini):
!   1   alpha            degree of cure
!   2   Tg               glass transition temperature
!   3-5 cum. thermal strain (e11,e22,e33)        diagnostic only
!   6-8 cum. chemical strain (e11,e22,e33)       diagnostic only
!   9   E11
!  10   E22
!  11   curerate         dalpha/dt evaluated at end of increment
!  12   ag
!  13   aInf
!  14   VgInf
!  15   aHg
!  16   aCg
!  17   alphaf
!  18   Ef
!  19   nuf
!  20   num
!  21   Em
!  22   Vf
!  23   T_prev           used to reconstruct dT (CalculiX UMAT does
!                        not pass DTEMP)
!  24-26 aX0,aX1,aX2
!  27-30 a1..a4
!  31-35 Tp1..Tp5
!  36   Ts = T - Tg
!  37   dT/dt
!
! ELCONLOC LAYOUT (40 constants, kinetics first, mechanical second):
!   1  A           pre-exponential factor
!   2  Ea          activation energy
!   3  n           Cole-Kamal exponent 1
!   4  m           Cole-Kamal exponent 2
!   5  D           diffusion constant
!   6  ac0         initial critical DoC
!   7  acT         critical DoC temperature slope
!   8  eta         Dibenedetto coefficient
!   9  TgInf       ultimate Tg
!  10  Tg0         initial Tg
!  11  H           total enthalpy of reaction
!  12  pr0         density of uncured resin
!  13  prInf       density of cured resin
!  14  pf          density of fibres
!  15  Em0         resin uncured E
!  16  EmI         resin fully cured E
!  17  ag          gel point / load-transfer DoC
!  18  aInf        DoC at which final stiffness is reached
!  19  num         resin Poisson ratio
!  20  aHg         resin CTE heat-up glassy
!  21  aCg         resin CTE cool-down glassy
!  22  alphaf      fibre CTE (isotropic)
!  23  Ef          fibre E
!  24  G12f        fibre G12
!  25  G23f        fibre G23
!  26  Vf          fibre volume fraction
!  27  nuf         fibre Poisson ratio
!  28  VgInf       total load-transferring shrinkage (volumetric)
!  29-31 aX0,aX1,aX2  parabolic CTE-vs-DoC at high Ts
!  32-35 a1..a4    CTE piecewise slopes between Tp1..Tp5
!  36-40 Tp1..Tp5  Ts breakpoints for piecewise CTE
!=======================================================================

      module cure_state
        implicit none
!       Per (integration point, element) cure state visible to dflux.f.
!       Sized large enough to cover any reasonable curing model; bump
!       max_el if the mesh is bigger. Allocation avoided to keep the
!       module thread-friendly under CalculiX's pragma usage.
        integer, parameter :: max_el = 200000
        integer, parameter :: max_ip = 27
        real*8 :: alpha_arr (max_ip, max_el) = 0.0d0
        real*8 :: dalpha_arr(max_ip, max_el) = 0.0d0
        real*8 :: pc_arr    (max_ip, max_el) = 0.0d0
!       Pre-assembled volumetric exotherm source [W/m^3]:
!         q_exo = (1 - Vf) * pr(alpha) * H * dalpha/dt
!       Published by umat_user so dflux.f can read it directly without
!       knowing material constants.
        real*8 :: q_exo_arr (max_ip, max_el) = 0.0d0
      end module cure_state

!=======================================================================
      subroutine umat_user(amat,iel,iint,kode,elconloc,emec,emec0,
     &  beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,icmd,ielas,
     &  mi,nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab,
     &  pnewdt,depvisc)
!
      use cure_state
      implicit none
!
      character*80 amat
      integer iel,iint,kode,ithermal,icmd,ielas,mi(*),nstate_,
     &        iorien,i,j,k,kk
      real*8 elconloc(*),emec(6),emec0(6),beta(6),
     &       xokl(3,3),xkl(3,3),voj,vj,t1l,dtime,time,ttime,
     &       xstateini(nstate_,mi(1),*),xstate(nstate_,mi(1),*),
     &       stre(6),stiff(21),pgauss(3),orab(7,*),pnewdt,depvisc
!
!     ----- material constants -----
      real*8 A,Ea,nck,mck,D,ac0,acT,eta,TgInf,Tg0,H,
     &       pr0,prInf,pf
      real*8 Em0,EmI,ag,aInf,num,aHg,aCg,alphaf,Ef,G12f,G23f,Vf,nuf,
     &       VgInf,aX0,aX1,aX2,a1,a2,a3,a4,Tp1,Tp2,Tp3,Tp4,Tp5
!
!     ----- working scalars -----
      real*8 alpha_n,alpha_np1,curerate,Tg,K_arr,Kd,Em,Gm,
     &       Kf,Km,Kt,G12,G23,E11,E22,v12,v23,v21,
     &       T_prev,dT,Ts,dTdt,
     &       aHHT,aCHT,alpham,alpha1,alpha2,alpha3,
     &       deth1,deth2,deth3,DVg,dechm,dech1,dech2,dech3,
     &       Ash,q,Ef_l,G12f_l,G23f_l,
     &       pr,pc,wf,denom,dstran(6),dstran_el(6)
      real*8 c(6,6),R
      parameter (R = 8.3145d0)
!
!=======================================================================
!     1.  Unpack material constants
!=======================================================================
      A      = elconloc(1)
      Ea     = elconloc(2)
      nck    = elconloc(3)
      mck    = elconloc(4)
      D      = elconloc(5)
      ac0    = elconloc(6)
      acT    = elconloc(7)
      eta    = elconloc(8)
      TgInf  = elconloc(9)
      Tg0    = elconloc(10)
      H      = elconloc(11)
      pr0    = elconloc(12)
      prInf  = elconloc(13)
      pf     = elconloc(14)
!
      Em0    = elconloc(15)
      EmI    = elconloc(16)
      ag     = elconloc(17)
      aInf   = elconloc(18)
      num    = elconloc(19)
      aHg    = elconloc(20)
      aCg    = elconloc(21)
      alphaf = elconloc(22)
      Ef     = elconloc(23)
      G12f   = elconloc(24)
      G23f   = elconloc(25)
      Vf     = elconloc(26)
      nuf    = elconloc(27)
      VgInf  = elconloc(28)
      aX0    = elconloc(29)
      aX1    = elconloc(30)
      aX2    = elconloc(31)
      a1     = elconloc(32)
      a2     = elconloc(33)
      a3     = elconloc(34)
      a4     = elconloc(35)
      Tp1    = elconloc(36)
      Tp2    = elconloc(37)
      Tp3    = elconloc(38)
      Tp4    = elconloc(39)
      Tp5    = elconloc(40)
!
      Ash = 0.0d0
!
!=======================================================================
!     2.  Read previous-step state and integrate cure (forward Euler
!         with the rate computed at the previous converged state, same
!         scheme as the original UMATHT).
!=======================================================================
      alpha_n  = xstateini(1 ,iint,iel)
      Tg       = xstateini(2 ,iint,iel)
      curerate = xstateini(11,iint,iel)
      T_prev   = xstateini(23,iint,iel)
!
      if (alpha_n.lt.1.0d-6) alpha_n = 1.0d-6
      if (T_prev .le.0.0d0 ) T_prev  = t1l       ! first increment
!
      alpha_np1 = alpha_n + dtime*curerate
      if (alpha_np1.gt.0.999d0) alpha_np1 = 1.0d0
!
!     Re-evaluate kinetics at (alpha_{n+1}, T_{n+1})
      Tg = Tg0 + (eta*alpha_np1*(TgInf-Tg0))
     &           /(1.0d0 - (1.0d0-eta)*alpha_np1)
      K_arr = A*dexp(-Ea/(R*(t1l+273.15d0)))
      Kd    = dexp(D*(alpha_np1 - ac0 - acT*(t1l+273.15d0)))
      curerate = K_arr*(alpha_np1**mck)*((1.0d0-alpha_np1)**nck)
     &                /(1.0d0 + Kd)
      if (alpha_np1.ge.1.0d0) curerate = 0.0d0
!
!     Composite density (rule of mixtures) -- needed by dflux for the
!     exotherm source term qv = (1-Vf)*H*pr/pc * curerate
      pr = pr0*alpha_np1 + (1.0d0 - alpha_np1)*prInf
      pc = Vf*pf + (1.0d0 - Vf)*pr
!
!     Publish to the cure_state module for dflux.f
      if (iel.le.max_el .and. iint.le.max_ip) then
         alpha_arr (iint,iel) = alpha_np1
         dalpha_arr(iint,iel) = curerate
         pc_arr    (iint,iel) = pc
         q_exo_arr (iint,iel) = (1.0d0 - Vf)*pr*H*curerate
      endif
!
      dT   = t1l - T_prev
      dTdt = 0.0d0
      if (dtime.gt.0.0d0) dTdt = dT/dtime
      Ts   = t1l - Tg
!
!=======================================================================
!     3.  Cure-dependent resin modulus (CHILE)
!=======================================================================
      if (alpha_np1.lt.ag) then
        Em = Em0
      else if (alpha_np1.ge.ag .and. t1l.lt.Tg .and.
     &         alpha_np1.lt.aInf) then
        Em = Em0*(1.0d0-(alpha_np1-ag)/(aInf-ag))
     &     + EmI* (alpha_np1-ag)/(aInf-ag)
      else if (alpha_np1.ge.aInf .and. t1l.lt.Tg) then
        Em = EmI
      else
        Em = 0.01d0*EmI                ! rubbery state above Tg
      endif
      Gm = 0.5d0*Em/(1.0d0+num)
!
!     Compression / pre-gel softening of fibre props (matches original
!     UMAT lines 231-235; uses start-of-increment strain for sign).
      q     = -10.0d0
      Ef_l   = Ef
      G12f_l = G12f
      G23f_l = G23f
      if (emec0(1).lt.0.0d0 .or. alpha_np1.lt.ag) then
        Ef_l   = (0.01d0 + 0.99d0/(1.0d0+dexp(q*(alpha_np1-ag))))*Ef
        G12f_l = (0.01d0 + 0.99d0/(1.0d0+dexp(q*(alpha_np1-ag))))*G12f
        G23f_l = (0.01d0 + 0.99d0/(1.0d0+dexp(q*(alpha_np1-ag))))*G23f
      endif
!
!=======================================================================
!     4.  Halpin-Tsai / self-consistent micromechanics for orthotropic
!         lamina (transverse isotropic; identical formulae to the
!         original UMAT).
!=======================================================================
      Kf  = Ef_l/(2.0d0*(1.0d0 - nuf - 2.0d0*nuf**2))
      Km  = Em  /(2.0d0*(1.0d0 - num - 2.0d0*num**2))
      Kt  = ((Kf+Gm)*Km + (Kf-Km)*Gm*Vf)
     &      /(Kf+Gm - (Kf-Km)*Vf)
      G12 = Gm*(G12f_l+Gm + (G12f_l-Gm)*Vf)
     &       /(G12f_l+Gm - (G12f_l-Gm)*Vf)
      G23 = Gm*(Km*(Gm+G23f_l) + 2.0d0*Gm*G23f_l
     &         + Km*(G23f_l-Gm)*Vf)
     &      /(Km*(Gm+G23f_l) + 2.0d0*Gm*G23f_l
     &        - (Km+2.0d0*Gm)*(G23f_l-Gm)*Vf)
      E11 = Ef_l*Vf + Em*(1.0d0-Vf)
     &      + (4.0d0*(num-nuf**2)*Kf*Km*Gm*(1.0d0-Vf)*Vf)
     &        /((Kf+Gm)*Km + (Kf-Km)*Gm*Vf)
      v12 = nuf*Vf + num*(1.0d0-Vf)
     &      + ((num-nuf)*(Km-Kf)*Gm*(1.0d0-Vf)*Vf)
     &        /((Kf+Gm)*Km + (Kf-Km)*Gm*Vf)
      E22 = 1.0d0/(0.25d0/Kt + 0.25d0/G23 + (v12**2)/E11)
      v23 = (2.0d0*E11*Kt - E11*E22 - 4.0d0*(v12**2)*Kt*E22)
     &      /(2.0d0*E11*Kt)
      v21 = v12*E22/E11
!
!=======================================================================
!     5.  Build full 3D orthotropic stiffness in material axes.
!         CalculiX UMAT is always 3D (plane reduction is done at the
!         element level), so the NTENS=1/3/4 branches in the original
!         are dropped.
!=======================================================================
      do j=1,6
        do i=1,6
          c(i,j) = 0.0d0
        enddo
      enddo
      denom = (1.0d0 - v23 - v12*v21)*(1.0d0 + v23)
      c(1,1) = E11*(1.0d0 - v23**2)/denom
      c(1,2) = E11*(v21 + v21*v23)/denom
      c(1,3) = c(1,2)
      c(2,1) = E22*(v12 + v12*v23)/denom
      c(2,2) = E22*(1.0d0 - v12*v21)/denom
      c(2,3) = E22*(v23 + v12*v21)/denom
      c(3,1) = c(2,1)
      c(3,2) = c(2,3)
      c(3,3) = c(2,2)
      c(4,4) = G23
      c(5,5) = G12
      c(6,6) = G12
!
!=======================================================================
!     6.  Eigenstrain increments (thermal + chemical) -- folded UEXPAN.
!         Sign convention: positive expansion strains.
!=======================================================================
      aHHT = aX2*alpha_np1**2 + aX1*alpha_np1 + aX0
      aCHT = aHHT
!
      if (alpha_np1.lt.ag) then
        alpham = 0.0d0
      else
        if (dTdt.ge.0.0001d0) then                 ! heat-up
          if      (Ts.ge.Tp5) then
            alpham = aHHT
          else if (Ts.ge.Tp4) then
            alpham = a4*(Ts-Tp4) + a3*(Tp4-Tp3) + a2*(Tp3-Tp2)
     &              + a1*(Tp2-Tp1) + aHg
          else if (Ts.ge.Tp3) then
            alpham = a3*(Ts-Tp3) + a2*(Tp3-Tp2) + a1*(Tp2-Tp1) + aHg
          else if (Ts.ge.Tp2) then
            alpham = a2*(Ts-Tp2) + a1*(Tp2-Tp1) + aHg
          else if (Ts.ge.Tp1) then
            alpham = a1*(Ts-Tp1) + aHg
          else
            alpham = aHg
          endif
        else                                       ! cool-down
          if      (Ts.ge.Tp4) then
            alpham = aCHT
          else if (Ts.ge.Tp3) then
            alpham = a3*(Ts-Tp3) + a2*(Tp3-Tp2) + a1*(Tp2-Tp1) + aCg
          else if (Ts.ge.Tp2) then
            alpham = a2*(Ts-Tp2) + a1*(Tp2-Tp1) + aCg
          else if (Ts.ge.Tp1) then
            alpham = a1*(Ts-Tp1) + aCg
          else
            alpham = aCg
          endif
        endif
      endif
!
      alpha1 = (alphaf*Ef_l*Vf + alpham*Em*(1.0d0-Vf))
     &        /(Ef_l*Vf + Em*(1.0d0-Vf))
      alpha2 = (alphaf + nuf*alphaf)*Vf
     &        + (alpham + num*alpham)*(1.0d0-Vf)
     &        - (nuf*Vf + num*(1.0d0-Vf))*alpha1
      alpha3 = alpha2
!
      deth1 = alpha1*dT
      deth2 = alpha2*dT
      deth3 = alpha3*dT
!
!     Cure shrinkage increment. The factor-of-1/2 in the original
!     UEXPAN was a workaround for ABAQUS calling UEXPAN twice per
!     iteration; do NOT include it here.
      if (alpha_np1.lt.ag) then
        DVg = 0.0d0
      else
        DVg = curerate*dtime*(Ash/(aInf-ag)
     &        + (2.0d0*(VgInf-Ash)*(alpha_np1-ag))/((aInf-ag)**2))
      endif
      dechm = ((1.0d0+DVg)**(1.0d0/3.0d0)) - 1.0d0
      if (Vf.gt.0.95d0) then
        dech1 = 0.0d0
        dech2 = 0.0d0
        dech3 = 0.0d0
      else
        dech1 = (dechm*Em*(1.0d0-Vf))/(Ef_l*Vf + Em*(1.0d0-Vf))
        dech2 = (dechm + num*dechm)*(1.0d0-Vf)
     &          - (nuf*Vf + num*(1.0d0-Vf))*dech1
        dech3 = dech2
      endif
!
!=======================================================================
!     7.  Strain increment from CalculiX's mechanical strains, then
!         subtract eigenstrains. emec is delivered without thermal
!         pre-subtraction because no *EXPANSION is declared.
!=======================================================================
      do k=1,6
        dstran(k) = emec(k) - emec0(k)
      enddo
      dstran_el(1) = dstran(1) - deth1 - dech1
      dstran_el(2) = dstran(2) - deth2 - dech2
      dstran_el(3) = dstran(3) - deth3 - dech3
      dstran_el(4) = dstran(4)
      dstran_el(5) = dstran(5)
      dstran_el(6) = dstran(6)
!
!=======================================================================
!     8.  Stress update.
!=======================================================================
      do k=1,6
        do i=1,6
          stre(k) = stre(k) + c(k,i)*dstran_el(i)
        enddo
      enddo
!
!=======================================================================
!     9.  Tangent stiffness packing.
!         CalculiX expects upper-triangular column-major storage:
!           stiff(1)=c(1,1)
!           stiff(2)=c(1,2)  stiff(3)=c(2,2)
!           stiff(4)=c(1,3)  stiff(5)=c(2,3)  stiff(6)=c(3,3)
!           stiff(7)=c(1,4)  stiff(8)=c(2,4)  stiff(9)=c(3,4)  stiff(10)=c(4,4)
!           stiff(11)=c(1,5) ...                                stiff(15)=c(5,5)
!           stiff(16)=c(1,6) ...                                stiff(21)=c(6,6)
!         When icmd==3, only the stress is required and stiff is left
!         untouched.
!=======================================================================
      if (icmd.ne.3) then
        kk = 0
        do j=1,6
          do i=1,j
            kk = kk + 1
            stiff(kk) = c(i,j)
          enddo
        enddo
      endif
!
!=======================================================================
!    10.  Persist state.
!=======================================================================
      xstate(1 ,iint,iel) = alpha_np1
      xstate(2 ,iint,iel) = Tg
      xstate(3 ,iint,iel) = xstateini(3,iint,iel) + deth1
      xstate(4 ,iint,iel) = xstateini(4,iint,iel) + deth2
      xstate(5 ,iint,iel) = xstateini(5,iint,iel) + deth3
      xstate(6 ,iint,iel) = xstateini(6,iint,iel) + dech1
      xstate(7 ,iint,iel) = xstateini(7,iint,iel) + dech2
      xstate(8 ,iint,iel) = xstateini(8,iint,iel) + dech3
      xstate(9 ,iint,iel) = E11
      xstate(10,iint,iel) = E22
      xstate(11,iint,iel) = curerate
      xstate(12,iint,iel) = ag
      xstate(13,iint,iel) = aInf
      xstate(14,iint,iel) = VgInf
      xstate(15,iint,iel) = aHg
      xstate(16,iint,iel) = aCg
      xstate(17,iint,iel) = alphaf
      xstate(18,iint,iel) = Ef
      xstate(19,iint,iel) = nuf
      xstate(20,iint,iel) = num
      xstate(21,iint,iel) = Em
      xstate(22,iint,iel) = Vf
      xstate(23,iint,iel) = t1l
      xstate(24,iint,iel) = aX0
      xstate(25,iint,iel) = aX1
      xstate(26,iint,iel) = aX2
      xstate(27,iint,iel) = a1
      xstate(28,iint,iel) = a2
      xstate(29,iint,iel) = a3
      xstate(30,iint,iel) = a4
      xstate(31,iint,iel) = Tp1
      xstate(32,iint,iel) = Tp2
      xstate(33,iint,iel) = Tp3
      xstate(34,iint,iel) = Tp4
      xstate(35,iint,iel) = Tp5
      xstate(36,iint,iel) = Ts
      xstate(37,iint,iel) = dTdt
!
      return
      end
