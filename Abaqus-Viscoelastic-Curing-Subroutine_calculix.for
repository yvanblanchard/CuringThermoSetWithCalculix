!=======================================================================
! Abaqus-Viscoelastic-Curing-Subroutine_calculix.for
!=======================================================================
! CalculiX port of Abaqus-Viscoelastic-Curing-Subroutine.for.
!
! CalculiX has NO equivalents for UEXPAN, USDFLD, DISP or HETVAL.
! All Abaqus subroutine roles are consolidated here:
!
!   Abaqus USDFLD  --> Section 2 of umat_user (cure kinetics)
!   Abaqus UEXPAN  --> Section 6 of umat_user (eigenstrain)
!   Abaqus HETVAL  --> dflux subroutine (same as dflux.f)
!   Abaqus DISP    --> Not needed: use *BOUNDARY + *AMPLITUDE in .inp
!   Abaqus UMAT    --> Section 4-5 of umat_user (Prony viscoelastic)
!
! Model: transverse-isotropic carbon/epoxy, 9-term Prony series
!        with cure-shifted WLF reduced time (DiBenedetto Tg).
!        Cure kinetics: dual-Arrhenius (Kamal-Sourour / first-order).
!        Eigenstrain: orthotropic CTE + chemical shrinkage.
!
! Compilation (replaces umat_user.f + dflux.f):
!   gfortran -shared -fPIC -O2 -ffixed-form -ffixed-line-length-none
!            -Wl,--export-all-symbols
!            -J _ccx_build -I _ccx_build
!            Abaqus-Viscoelastic-Curing-Subroutine_calculix.for
!            -o UserSubroutines.dll
!
! Input deck requirements:
!   *USER MATERIAL, CONSTANTS=18
!    <18 values -- see elconloc layout below>
!   *DEPVAR
!    64
!   *INITIAL CONDITIONS, TYPE=SOLUTION
!    COMPOSITE, 1.0e-3, 253.15, 0.,0.,0., 0.,  <T_prev>, 0.,0.,0., <T>,
!               0., ... (64 values total)
!   *CONDUCTIVITY, TYPE=ORTHO
!    5.43, 0.4135, 0.4135
!   *SPECIFIC HEAT
!    862.
!   *DENSITY
!    1578.
!   *DFLUX
!    COMPOSITE, BFNU, 1.0
!
! ELCONLOC layout (18 constants):
!    1  A1_kin   Arrhenius pre-exp 1  (alpha <= 0.3 branch)   [1/min]
!    2  A2_kin   Arrhenius pre-exp 2  (alpha <= 0.3 branch)   [1/min]
!    3  A3_kin   Arrhenius pre-exp 3  (alpha > 0.3 branch)    [1/min]
!    4  DE1      Activation energy 1                          [J/mol]
!    5  DE2      Activation energy 2                          [J/mol]
!    6  DE3      Activation energy 3                          [J/mol]
!    7  B        Autocatalytic constant                       [-]
!    8  Tg0      Initial glass-transition temperature         [K]
!    9  TgInf    Ultimate glass-transition temperature        [K]
!   10  eta_dB   DiBenedetto fitting parameter lambda        [-]
!   11  RHOR     Resin density (cured)                       [kg/m^3]
!   12  HR       Total enthalpy of reaction                  [J/kg]
!   13  VF       Fibre volume fraction                       [-]
!   14  CTE1     Thermal exp. coeff., fibre direction        [1/K]
!   15  CTE2     Thermal exp. coeff., transverse             [1/K]
!   16  CCS1     Chem. shrinkage coeff., fibre direction     [-]
!   17  CCS2     Chem. shrinkage coeff., transverse          [-]
!   18  alpha_gel Gelation degree of cure                    [-]
!
! Suggested *USER MATERIAL line for AS4/3501-6 representative values:
!   2.101D9, -2.014D9, 1.960D5, 8.07D4, 7.78D4, 5.66D4, 0.47,
!   253.15, 493.15, 0.50, 1272., 473.6D3, 0.52,
!   0.5D-6, 35.3D-6, -167.D-6, -8810.D-6, 0.30
!
! XSTATE layout (64 state variables):
!    1   alpha         degree of cure                   [-]
!    2   Tg            glass-transition temperature     [K]
!    3   cum_eigen_11  cumulative eigenstrain, dir 1    [-]
!    4   cum_eigen_22  cumulative eigenstrain, dir 2    [-]
!    5   cum_eigen_33  cumulative eigenstrain, dir 3    [-]
!    6   curerate      dalpha/dt                        [1/s]
!    7   (reserved)
!    8   (reserved)
!    9   (reserved)
!   10   T_prev        temperature at end of prev incr  [K]
!   11-64 q(K,I) Prony arm hereditary integrals,
!         K=1..6 (Voigt stress index), I=1..9 (arm index)
!         stored column-major: q(1,1),q(2,1),...,q(6,1),
!                               q(1,2),q(2,2),...,q(6,9)
!=======================================================================

!=======================================================================
! cure_state module -- shared between umat_user and dflux
!=======================================================================
      module cure_state
        implicit none
        integer, parameter :: max_el = 200000
        integer, parameter :: max_ip = 27
        real*8 :: q_exo_arr(max_ip, max_el) = 0.0d0
      end module cure_state

!=======================================================================
! umat_user -- CalculiX user material (viscoelastic curing composite)
!=======================================================================
      subroutine umat_user(amat,iel,iint,kode,elconloc,emec,emec0,
     &  beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,icmd,
     &  ielas,mi,nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,
     &  orab,pnewdt,depvisc)
!
      use cure_state
      implicit none
!
      character*80 amat
      integer iel,iint,kode,ithermal,icmd,ielas,mi(*),nstate_,
     &        iorien,I,J,K,KK,Countnum
      real*8 elconloc(*),emec(6),emec0(6),beta(6),
     &       xokl(3,3),xkl(3,3),voj,vj,t1l,dtime,time,ttime,
     &       xstateini(nstate_,mi(1),*),xstate(nstate_,mi(1),*),
     &       stre(6),stiff(21),pgauss(3),orab(7,*),pnewdt,depvisc
!
!     ---- material constants from elconloc --------------------
      real*8 A1_kin,A2_kin,A3_kin,DE1,DE2,DE3,B_cat
      real*8 Tg0,TgInf,eta_dB,RHOR,HR,VF
      real*8 CTE1,CTE2,CCS1,CCS2,alpha_gel
!
!     ---- Prony series (hardcoded, 9 terms) -------------------
!     Unrelaxed (instantaneous) stiffness [Pa]  AS4/3501-6
      real*8 C11u,C12u,C22u,C23u,C44u,C55u
      parameter (C11u=127.4D9, C12u=3.88D9, C22u=10.0D9,
     &           C23u=3.0D9,   C44u=3.5D9,  C55u=5.0D9)
!     Fully relaxed: 14 % of unrelaxed
      real*8 C11inf,C12inf,C22inf,C23inf,C44inf,C55inf
      parameter (C11inf=C11u*0.14D0, C12inf=C12u*0.14D0,
     &           C22inf=C22u*0.14D0, C23inf=C23u*0.14D0,
     &           C44inf=C44u*0.14D0, C55inf=C55u*0.14D0)
!     Delta (relaxing portion)
      real*8 dC11,dC12,dC22,dC23,dC44,dC55
      parameter (dC11=C11u-C11inf, dC12=C12u-C12inf,
     &           dC22=C22u-C22inf, dC23=C23u-C23inf,
     &           dC44=C44u-C44inf, dC55=C55u-C55inf)
!     Relaxation times [s] -- Kim & White 3501-6 resin master curve,
!     taken from example-Job-Abaqus-Viscoelastic-Curing.inp (2026-07-03).
!     old (placeholder) values: 1e-8,1e-6,1e-4,1e-2,1,1e2,1e4,1e6,1e8
      real*8 Tai(9)
      data Tai /29.2D0,  2.92D3,  1.82D5,  1.10D7,  2.83D8,
     &          7.94D9,  1.95D11, 3.32D12, 4.92D14/
!     Prony amplitudes, from the same reference deck (sum = 0.994).
!     Applied to both fiber- and transverse-direction relaxing parts;
!     the relaxed floor stays C_inf = 14% of unrelaxed (hardcoded above).
!     old A1p: 0.005,0.010,0.020,0.055,0.100,0.180,0.200,0.170,0.120
!     old A2p: 0.010,0.020,0.040,0.080,0.130,0.185,0.185,0.130,0.080
      real*8 A1p(9)
      data A1p /0.059D0,0.066D0,0.083D0,0.112D0,0.154D0,
     &          0.262D0,0.184D0,0.049D0,0.025D0/
      real*8 A2p(9)
      data A2p /0.059D0,0.066D0,0.083D0,0.112D0,0.154D0,
     &          0.262D0,0.184D0,0.049D0,0.025D0/
!
!     ---- WLF parameters (standard epoxy) --------------------
      real*8 C1WLF,C2WLF
      parameter (C1WLF=17.44D0, C2WLF=51.6D0)
      real*8 R_gas
      parameter (R_gas=8.3143D0)
!
!     ---- working variables -----------------------------------
      real*8 alpha_n,alpha_np1,curerate,Tg,Tg_new,T_prev,dT
      real*8 k1,k2,k3,dadt_min,ShiftFac,DZETA,dTg,gel
      real*8 d_eigen(6),cum_eigen(6),epsilonE(6)
      real*8 depsilon(6),depsilon_el(6)
      real*8 qold(6,9),q(6,9)
      real*8 f(6)
      real*8 sigmaE(6)
      real*8 c_eff(6,6)
      real*8 xi,expxi,gammai,gammai_arr(9)
      real*8 pr,pc,pr0_kin,prInf_kin
!     ---- orientation support (2026-07-03) ---------------------
!     CalculiX passes GLOBAL strains to user materials; with
!     *ORIENTATION (iorien>0) the material frame is built from
!     orab and all constitutive work is done in the local frame.
      real*8 skl(3,3),emecl(6),emec0l(6),vloc(6),col(6),cg(6,6)
      integer IQ
!
!=======================================================================
!     1. Unpack material constants
!=======================================================================
      A1_kin   = elconloc(1)
      A2_kin   = elconloc(2)
      A3_kin   = elconloc(3)
      DE1      = elconloc(4)
      DE2      = elconloc(5)
      DE3      = elconloc(6)
      B_cat    = elconloc(7)
      Tg0      = elconloc(8)
      TgInf    = elconloc(9)
      eta_dB   = elconloc(10)
      RHOR     = elconloc(11)
      HR       = elconloc(12)
      VF       = elconloc(13)
      CTE1     = elconloc(14)
      CTE2     = elconloc(15)
      CCS1     = elconloc(16)
      CCS2     = elconloc(17)
      alpha_gel= elconloc(18)
!
!     ---- rotate strains into the material (ply) frame ---------
      do K=1,6
          emecl(K)  = emec(K)
          emec0l(K) = emec0(K)
      end do
      if (iorien.gt.0) then
          call visc_frame(orab(1,iorien),skl)
          call visc_rotsym(emecl, skl,-1)
          call visc_rotsym(emec0l,skl,-1)
      end if
!
!=======================================================================
!     2. Read previous-step state
!=======================================================================
      alpha_n   = xstateini(1 ,iint,iel)
      Tg        = xstateini(2 ,iint,iel)
      cum_eigen(1) = xstateini(3,iint,iel)
      cum_eigen(2) = xstateini(4,iint,iel)
      cum_eigen(3) = xstateini(5,iint,iel)
      cum_eigen(4) = 0.0d0
      cum_eigen(5) = 0.0d0
      cum_eigen(6) = 0.0d0
      curerate  = xstateini(6 ,iint,iel)
      T_prev    = xstateini(10,iint,iel)
!
!     Seed first increment
      if (alpha_n  .lt. 1.0d-6) alpha_n  = 1.0d-6
      if (T_prev   .le. 0.0d0 ) T_prev   = t1l
!
!=======================================================================
!     3. Cure kinetics -- dual Arrhenius (Kamal-Sourour / first-order)
!        Kinetic constants fitted in [1/min]; convert to [1/s] at end.
!        Forward Euler integration using curerate from previous step.
!=======================================================================
      alpha_np1 = alpha_n + dtime*curerate
      alpha_np1 = min(alpha_np1, 1.0d0)
!
!     Rate constants at current temperature
      k2 = A2_kin*dexp(-DE2/R_gas/t1l)
      if (alpha_np1 .le. 0.3d0) then
          k1      = A1_kin*dexp(-DE1/R_gas/t1l)
          dadt_min= (k1 + k2*alpha_np1)*(1.0d0-alpha_np1)
     &                                  *(B_cat-alpha_np1)
      else
          k3      = A3_kin*dexp(-DE3/R_gas/t1l)
          dadt_min= k3*(1.0d0-alpha_np1)
      end if
!     Convert /min to /s
      curerate = dadt_min / 60.0d0
      if (alpha_np1 .ge. 1.0d0) curerate = 0.0d0
!
!     Cure-dependent Tg (DiBenedetto equation)
      Tg_new = Tg0 + (eta_dB*alpha_np1*(TgInf-Tg0))
     &               /(1.0d0-(1.0d0-eta_dB)*alpha_np1)
!
!     Publish exotherm to cure_state for dflux subroutine
!     q_exo = (1-Vf) * rho_r(alpha) * H * dalpha/dt
!     Use simple linear density interpolation
      pr0_kin  = RHOR * 0.97d0   ! approximate uncured density
      prInf_kin= RHOR
      pr = pr0_kin*alpha_np1 + (1.0d0-alpha_np1)*prInf_kin
      pc = VF*1790.0d0 + (1.0d0-VF)*pr   ! composite density
      if (iel.le.max_el .and. iint.le.max_ip) then
          q_exo_arr(iint,iel) = (1.0d0-VF)*pr*HR*curerate
      end if
!
!=======================================================================
!     4. WLF shift factor (cure-shifted, relative to DiBenedetto Tg)
!=======================================================================
!     SIGN FIX 2026-07-03: log aT = -C1*(T-Tg)/(C2+T-Tg).
!     Above Tg (rubbery) aT << 1 -> reduced time dt/aT large -> fast
!     relaxation; below Tg aT >> 1 -> frozen. The previous +C1 exponent
!     (inherited from the reconstruction .for) froze the RUBBERY state
!     and let the glassy state relax -- physically inverted.
      dTg = t1l - Tg_new
      if (dTg .gt. -C2WLF*0.99d0) then
          ShiftFac = 10.0d0**( max(-30.0d0, min(30.0d0,
     &                -C1WLF*dTg / (C2WLF + dTg) )) )
      else
          ShiftFac = 1.0d30          ! deep glassy: frozen
      end if
      DZETA = dtime / ShiftFac       ! reduced time increment [s_reduced]
!
!=======================================================================
!     5. Pre-gelation sigmoid factor (avoids step discontinuity at gel)
!=======================================================================
      gel = 0.01d0 + 0.99d0/(1.0d0+dexp(-50.0d0*(alpha_np1-alpha_gel)))
!
!=======================================================================
!     6. Eigenstrain increments (replaces Abaqus UEXPAN)
!        Thermal: deth = CTE * dT
!        Chemical shrinkage: decs = CCS * d(alpha)  (negative on cure)
!=======================================================================
      dT = t1l - T_prev
!
      d_eigen(1) = CTE1*dT + CCS1*(alpha_np1-alpha_n)
      d_eigen(2) = CTE2*dT + CCS2*(alpha_np1-alpha_n)
      d_eigen(3) = CTE2*dT + CCS2*(alpha_np1-alpha_n)  ! transv. isotropy
      d_eigen(4) = 0.0d0
      d_eigen(5) = 0.0d0
      d_eigen(6) = 0.0d0
!
!     Update cumulative eigenstrain
      cum_eigen(1) = cum_eigen(1) + d_eigen(1)
      cum_eigen(2) = cum_eigen(2) + d_eigen(2)
      cum_eigen(3) = cum_eigen(3) + d_eigen(3)
!
!=======================================================================
!     7. Mechanical strain increment and total elastic strain
!        emec(k)  = total strain at end of increment (no *EXPANSION)
!        emec0(k) = total strain at start
!        Elastic = total - cumulative eigenstrain
!=======================================================================
      do K=1,6
          depsilon(K)    = emecl(K) - emec0l(K)
          depsilon_el(K) = depsilon(K) - d_eigen(K)
          epsilonE(K)    = emecl(K) - cum_eigen(K)
      end do
!
!=======================================================================
!     8. Read Prony arm internal variables from xstateini
!        Stored at xstate(11..64) column-major: q(K=1..6, I=1..9)
!=======================================================================
      Countnum = 11
      do I=1,9
          do K=1,6
              qold(K,I)= xstateini(Countnum,iint,iel)
              Countnum = Countnum + 1
          end do
      end do
!
!=======================================================================
!     9. Prony arm update (Wineman-Pipkin recursive algorithm)
!        Reduced time: dzeta = dt / aT(alpha,T)
!        Transverse-isotropic coupling (same as Abaqus UMAT):
!          sigma_1 depends on eps_1 (C11) and eps_2+eps_3 (C12)
!          sigma_2 depends on eps_1 (C12), eps_2 (C22), eps_3 (C23)
!          sigma_3 same by transverse isotropy
!          sigma_4 (shear 23): C44
!          sigma_5 (shear 13): C55
!          sigma_6 (shear 12): C55
!=======================================================================
      do I=1,9
          if (Tai(I) .gt. 1.0d-30) then
              xi = DZETA / Tai(I)
              if (xi .gt. 1.0d-8) then
                  expxi  = dexp(-xi)
                  gammai = (1.0d0 - expxi) / xi
              else
                  expxi  = 1.0d0 - xi + 0.5d0*xi**2
                  gammai = 1.0d0 - 0.5d0*xi + xi**2/6.0d0
              end if
          else
              expxi  = 0.0d0
              gammai = 0.0d0
          end if
          gammai_arr(I) = gammai
!
!         Scale Prony contribution by gel factor (pre-gel = near zero)
          q(1,I)= expxi*qold(1,I) + gammai*gel*(
     &            A1p(I)*dC11*depsilon_el(1)
     &          + A1p(I)*dC12*(depsilon_el(2)+depsilon_el(3)))
!
          q(2,I)= expxi*qold(2,I) + gammai*gel*(
     &            A2p(I)*dC12*depsilon_el(1)
     &          + A2p(I)*dC22*depsilon_el(2)
     &          + A2p(I)*dC23*depsilon_el(3))
!
          q(3,I)= expxi*qold(3,I) + gammai*gel*(
     &            A2p(I)*dC12*depsilon_el(1)
     &          + A2p(I)*dC23*depsilon_el(2)
     &          + A2p(I)*dC22*depsilon_el(3))
!
          q(4,I)= expxi*qold(4,I) +
     &            gammai*gel*A2p(I)*dC44*depsilon_el(4)
!
          q(5,I)= expxi*qold(5,I) +
     &            gammai*gel*A1p(I)*dC55*depsilon_el(5)
!
          q(6,I)= expxi*qold(6,I) +
     &            gammai*gel*A1p(I)*dC55*depsilon_el(6)
      end do
!
!=======================================================================
!    10. Accumulate hereditary integrals (sum over all 9 arms)
!=======================================================================
      do K=1,6
          f(K) = 0.0d0
      end do
      do I=1,9
          do K=1,6
              f(K) = f(K) + q(K,I)
          end do
      end do
!
!=======================================================================
!    11. Effective stress update (long-time stiffness + hereditary)
!        sigma = C_inf * epsilon_elastic_total + sum(q_I)
!        Scaled by gel factor (zero stiffness pre-gelation)
!=======================================================================
      sigmaE(1) = gel*(C11inf*epsilonE(1)
     &              + C12inf*(epsilonE(2)+epsilonE(3))) + f(1)
      sigmaE(2) = gel*(C12inf*epsilonE(1)
     &              + C22inf*epsilonE(2)
     &              + C23inf*epsilonE(3)) + f(2)
      sigmaE(3) = gel*(C12inf*epsilonE(1)
     &              + C23inf*epsilonE(2)
     &              + C22inf*epsilonE(3)) + f(3)
      sigmaE(4) = gel*C44inf*epsilonE(4) + f(4)
      sigmaE(5) = gel*C55inf*epsilonE(5) + f(5)
      sigmaE(6) = gel*C55inf*epsilonE(6) + f(6)
!
      do K=1,6
          stre(K) = sigmaE(K)
      end do
!     rotate stress back to the global frame
      if (iorien.gt.0) call visc_rotsym(stre,skl,+1)
!
!=======================================================================
!    12. Consistent tangent stiffness --> stiff(21) upper-triangular
!        C_eff(I,J) = C_inf(I,J) * gel + sum_arms{gammai*Ap*dC(I,J)*gel}
!
!        Packing (CalculiX convention, column-major upper triangle):
!          stiff( 1)=c(1,1)
!          stiff( 2)=c(1,2)  stiff( 3)=c(2,2)
!          stiff( 4)=c(1,3)  stiff( 5)=c(2,3)  stiff( 6)=c(3,3)
!          stiff( 7)=c(1,4)  stiff( 8)=c(2,4)  stiff( 9)=c(3,4)
!          stiff(10)=c(4,4)
!          stiff(11..15): column 5
!          stiff(16..21): column 6
!=======================================================================
      if (icmd .ne. 3) then
!
!         Initialise effective stiffness matrix
          do J=1,6
              do I=1,6
                  c_eff(I,J) = 0.0d0
              end do
          end do
!
!         Base: long-time (fully relaxed) stiffness * gel factor
          c_eff(1,1) = gel*C11inf
          c_eff(1,2) = gel*C12inf
          c_eff(2,1) = gel*C12inf
          c_eff(1,3) = gel*C12inf      ! C13 = C12 (transv. isotropy)
          c_eff(3,1) = gel*C12inf
          c_eff(2,2) = gel*C22inf
          c_eff(2,3) = gel*C23inf
          c_eff(3,2) = gel*C23inf
          c_eff(3,3) = gel*C22inf      ! C33 = C22
          c_eff(4,4) = gel*C44inf
          c_eff(5,5) = gel*C55inf
          c_eff(6,6) = gel*C55inf      ! C66 = C55
!
!         Add Prony relaxation contribution
          do I=1,9
              gammai = gammai_arr(I)
              c_eff(1,1)=c_eff(1,1)+gel*A1p(I)*dC11*gammai
              c_eff(1,2)=c_eff(1,2)+gel*A1p(I)*dC12*gammai
              c_eff(2,1)=c_eff(2,1)+gel*A2p(I)*dC12*gammai
              c_eff(1,3)=c_eff(1,3)+gel*A1p(I)*dC12*gammai
              c_eff(3,1)=c_eff(3,1)+gel*A2p(I)*dC12*gammai
              c_eff(2,2)=c_eff(2,2)+gel*A2p(I)*dC22*gammai
              c_eff(2,3)=c_eff(2,3)+gel*A2p(I)*dC23*gammai
              c_eff(3,2)=c_eff(3,2)+gel*A2p(I)*dC23*gammai
              c_eff(3,3)=c_eff(3,3)+gel*A2p(I)*dC22*gammai
              c_eff(4,4)=c_eff(4,4)+gel*A2p(I)*dC44*gammai
              c_eff(5,5)=c_eff(5,5)+gel*A1p(I)*dC55*gammai
              c_eff(6,6)=c_eff(6,6)+gel*A1p(I)*dC55*gammai
          end do
!
!         Rotate stiffness to the global frame: column IQ of the
!         global operator = global stress response to a unit global
!         strain, evaluated through the local operator. Exact w.r.t.
!         the code's own Voigt convention.
          if (iorien.gt.0) then
              do IQ=1,6
                  do K=1,6
                      vloc(K) = 0.0d0
                  end do
                  vloc(IQ) = 1.0d0
                  call visc_rotsym(vloc,skl,-1)
                  do K=1,6
                      col(K) = 0.0d0
                      do J=1,6
                          col(K) = col(K) + c_eff(K,J)*vloc(J)
                      end do
                  end do
                  call visc_rotsym(col,skl,+1)
                  do K=1,6
                      cg(K,IQ) = col(K)
                  end do
              end do
              do IQ=1,6
                  do K=1,6
                      c_eff(K,IQ) = cg(K,IQ)
                  end do
              end do
          end if
!
!         Pack upper-triangular column-major into stiff(21)
          KK = 0
          do J=1,6
              do I=1,J
                  KK = KK + 1
                  stiff(KK) = c_eff(I,J)
              end do
          end do
!
      end if   ! icmd .ne. 3
!
!=======================================================================
!    13. Persist state to xstate
!=======================================================================
      xstate(1 ,iint,iel) = alpha_np1
      xstate(2 ,iint,iel) = Tg_new
      xstate(3 ,iint,iel) = cum_eigen(1)
      xstate(4 ,iint,iel) = cum_eigen(2)
      xstate(5 ,iint,iel) = cum_eigen(3)
      xstate(6 ,iint,iel) = curerate
      xstate(7 ,iint,iel) = 0.0d0         ! reserved
      xstate(8 ,iint,iel) = 0.0d0         ! reserved
      xstate(9 ,iint,iel) = 0.0d0         ! reserved
      xstate(10,iint,iel) = t1l            ! T_prev for next increment
!
!     Prony arm internal variables
      Countnum = 11
      do I=1,9
          do K=1,6
              xstate(Countnum,iint,iel) = q(K,I)
              Countnum = Countnum + 1
          end do
      end do
!
      return
      end

!=======================================================================
! dflux -- CalculiX volumetric heat source (replaces Abaqus HETVAL)
!          Identical interface to dflux.f; reads q_exo from cure_state.
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
      flux(1) = 0.0d0
      iscale  = 0
!
      if (loadtype(1:2).ne.'BF') return
      if (noel.le.0 .or. noel.gt.max_el) return
      if (npt .le.0 .or. npt .gt.max_ip) return
!
      flux(1) = q_exo_arr(npt, noel)
!
      return
      end

!=======================================================================
! visc_frame -- build the local->global rotation matrix skl from a
!               CalculiX *ORIENTATION (rectangular) record orab(1:7):
!               a = orab(1:3), b = orab(4:6);
!               x' = a/|a|, z' = (a x b)/|a x b|, y' = z' x x'.
!               Columns of skl are the local axes in global coords,
!               so  t_global = skl * t_local * skl^T.
!=======================================================================
      subroutine visc_frame(orab7,skl)
      implicit none
      real*8 orab7(7),skl(3,3)
      real*8 a(3),b(3),x(3),y(3),z(3),rn
      integer i
      do i=1,3
          a(i)=orab7(i)
          b(i)=orab7(i+3)
      end do
      rn = dsqrt(a(1)**2+a(2)**2+a(3)**2)
      do i=1,3
          x(i)=a(i)/rn
      end do
      z(1)=a(2)*b(3)-a(3)*b(2)
      z(2)=a(3)*b(1)-a(1)*b(3)
      z(3)=a(1)*b(2)-a(2)*b(1)
      rn = dsqrt(z(1)**2+z(2)**2+z(3)**2)
      do i=1,3
          z(i)=z(i)/rn
      end do
      y(1)=z(2)*x(3)-z(3)*x(2)
      y(2)=z(3)*x(1)-z(1)*x(3)
      y(3)=z(1)*x(2)-z(2)*x(1)
      do i=1,3
          skl(i,1)=x(i)
          skl(i,2)=y(i)
          skl(i,3)=z(i)
      end do
      return
      end

!=======================================================================
! visc_rotsym -- rotate a symmetric tensor stored in CalculiX Voigt
!                order (11,22,33,12,13,23), tensor shear components.
!                idir=+1 : local -> global   B = R A R^T
!                idir=-1 : global -> local   B = R^T A R
!=======================================================================
      subroutine visc_rotsym(t6,R,idir)
      implicit none
      real*8 t6(6),R(3,3)
      integer idir
      real*8 A(3,3),B(3,3),Q(3,3)
      integer i,j,k,l
      A(1,1)=t6(1)
      A(2,2)=t6(2)
      A(3,3)=t6(3)
      A(1,2)=t6(4)
      A(2,1)=t6(4)
      A(1,3)=t6(5)
      A(3,1)=t6(5)
      A(2,3)=t6(6)
      A(3,2)=t6(6)
      if (idir.gt.0) then
          do i=1,3
              do j=1,3
                  Q(i,j)=R(i,j)
              end do
          end do
      else
          do i=1,3
              do j=1,3
                  Q(i,j)=R(j,i)
              end do
          end do
      end if
      do i=1,3
          do j=1,3
              B(i,j)=0.0d0
              do k=1,3
                  do l=1,3
                      B(i,j)=B(i,j)+Q(i,k)*A(k,l)*Q(j,l)
                  end do
              end do
          end do
      end do
      t6(1)=B(1,1)
      t6(2)=B(2,2)
      t6(3)=B(3,3)
      t6(4)=B(1,2)
      t6(5)=B(1,3)
      t6(6)=B(2,3)
      return
      end
