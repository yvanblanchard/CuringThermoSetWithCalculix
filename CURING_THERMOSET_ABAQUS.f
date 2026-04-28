      ! CASH2D: Cure Advancement, Shrinkage & Heat model
      ! User subroutines to simulate all aspects of 
      ! thermoset cure behaviour within a full coupled 
      ! thermal-stress abaqus simulation using 2D section.
      ! The model consists of:
      ! An UMATHT - User-defined material for Heat-Transfer analysis
      !! This subroutine computes the curing, glass transition tem-
      !! rature as well as conductivity, heat capacity in relation to
      !! temperature, degree of cure and based on energy balance, the 
      !! development of exotherm in a composite material.
C
      ! An UEXPAN - User-defined expansion and shrinkage routine
      !! Based on the thermal behaviour from UMATHT the thermal expan-
      !! sion and the cure shrinkage in the composite is computed. The
      !! routine takes into account the load-transferring part of the 
      !! resin expansion and shrinkage, as well as the fibre expansion,
      !! pre- and post-gel in orthotropic material directions. 
C
      ! An UMAT - User-defined material for cure stiffness development
      !! An orthotropic (transverse isotropic) material, using the CHI-
      !! LE model for the resin stiffness build-up during curing is im-
      !! plemented.
C
      ! UMATHT is used for Cure Advanced & Heat
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
C
      ! Define paramters used
      integer i
      real*8 A,E,n,m,D,ac0,acT,R,H
      real*8 eta,TgInf,Tg0,Tg
      real*8 K,Kd,curerate,ak
      real*8 pr,prInf,pr0,pf,pc,vf,wf,sqVf
      real*8 cpr0,cprT,cpf0,cpfT,cpr,cpf,cp,rk
      real*8 kr,klf,klf0,ktf,pi,kr0,krT,kra
      real*8 akr,bkr,ckr,dkr
      real*8 K11,K22,K33,cond(2)
      real*8 DcondDT(2)
      real*8 theta
C
      ! Define Properties
      A = PROPS(1)       ! Pre-exponential factor
      E = PROPS(2)       ! Activation energy
      n = PROPS(3)       ! First Power law constant
      m = PROPS(4)       ! Second power law
      D = PROPS(5)       ! Diffusion constant
      ac0 = PROPS(6)     ! Initial critical DoC
      acT = PROPS(7)     ! Critical DoC Increase
      eta = PROPS(8)     ! Fitting coefficient Dibenedetto
      TgInf = PROPS(9)   ! Ultimate glass transition temperature
      Tg0 = PROPS(10)    ! Initial glass transition temperature
      H = PROPS(11)      ! Total Enthalpy of reaction
      Vf = PROPS(12)     ! FVF
      pr0 = PROPS(13)    ! Density of uncured resin
      prInf = PROPS(14)  ! Density of cured resin
      pf = PROPS(15)     ! Density of fibres
C
      theta= 0.0 !Ply orientation passed from model
C
      ! Specific heat capacity definition
      cprT = 3.0    ! Increase in heat capacity resin
      cpr0 = 1900.  ! Initial Heat capacity resin
      cpf0 = 810.   ! Initial heat capacity fibre
      cpfT = 2.05   ! Increase in heat capacity fibre
      R = 8.3145    ! Universal gas constant
      klf0 = 1.3    ! Initial Longitudinal fibre conductivity
      klfT = 0.0156 ! Fibre conductivity increase - temperature
      kr0 = 0.135   ! Initial resin conductivity
      krT = 0.000343! Resin conductivity increase - temperature
      kra = 0.00607 ! Resin conductivity increase - Degree of cure
      pr0 = 1088.   ! Initial density of resin
C
      !--state variables from previous increment--
      ak=statev(1)
      curerate=statev(11)
      IF (ak<0.000001) THEN
       ak=0.000001
      ENDIF        
C
      !---degree of cure at the end of increment--
      ak=ak+dtime*curerate
C
      !--cure kinetics---------
      Tg=Tg0+(eta*ak*(TgInf-Tg0))/(1.0-(1.0-eta)*ak) ! Dibenedetto eq
      K=A*exp(-1.0*E/(R*(temp+273.15))) ! Reaction model
      Kd=exp(D*(ak-ac0-acT*(temp+273.15))) ! Diffusion model
      curerate=K*(ak**m)*((1.0-ak)**n)/(1.0+Kd) ! Cole model
      ! Check for max degree of cure
      IF (ak>0.999) THEN
       curerate=0.0
       ak=1.0
      ENDIF
C
      !--specific heat capacity---------
      ! Heat capacity of resin func of temperature    
      cpr = cprT*temp+cpr0
C
    ! Heat capacity of fibre
      cpf = cpf0 !cpfT*temp+cpf0
C  
      ! Define density of composite - ROM
      pr = pr0*ak+(1.0-ak)*prInf
      pc = vf*pf+(1.0-vf)*pr
C     
    !Define FWF for composite heat capacity
      wf = vf*pf/pc
      cp = wf*cpf+(1.0-wf)*cpr
C
      ! Compute constituent conductivity
      kr = kr0+krT*TEMP+kra*ak
      klf = klf0!klf0+klfT*TEMP
      rk = kr/klf ! Allocate ratio of klf and kr
      sqVf = Vf**(1.0/2.0) ! Allocate the squareroot of FVF
C
      !Compute direction https://doi.org/10.1520/CTR10326J 
      K11 = vf*klf+(1.0-vf)*kr
      K22 = (1.0-sqVf)*kr+(sqVf*kr)/(1.0-sqVf*(1.0-kr/klf))
      K33 = K22
C
      !--thermal conductivity matrix----
      cond(1)=K11  !x-axis--
      cond(2)=K22  !y-axis--
      ! Thermal conductivity derivatives matrix
      DcondDT(1) = 0 ! Vf*klfT+(1.0-Vf)*krT
      DcondDT(2) = 0 ! (1.0-sqVf)*krT+(sqVf*krT)/(1.0-sqVf*(1.0-rk))
      ! 1 +(Vf*kr*(-krT/klf+(kr*klfT)/(klf)**2.0))
      ! 2 /(1.0-sqVf*(1.0-rk))**2.0
      !DcondDT(3)=DcondDT(2)
C
      !---energy balance-----------
      DUDT=cp
      DU=DUDT*DTEMP
      U=U+DU-(1.0-vf)*H*curerate*dtime*pr/pc
C
      ! Loop through energy balance
      DO i=1,NTGRD
         FLUX(i)=-COND(i)*DTEMDX(i)
         DFDG(i,i)=-cond(i)
         DFDT(i)=-DcondDT(i)*DTEMDX(i)
      END DO
C
      !--store state variables--
      statev(1) = ak
      statev(2) = Tg
      statev(11) = curerate
      statev(23) = NTGRD
C      
      RETURN
      END
C
      ! UMAT is used for cure dependent mechanical properties
      ! For example: https://abaqus-docs.mit.edu/2017/English
      ! /SIMACAESUBRefMap/simasub-c-umat.htm
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C     
      !Declaration of variables for the UMAT
      REAL*8 Em,Em0,EmI,Gm,num,Ef,G12f,G23f,nuf,Vf,aHg,aCg
     1 alphaf,lamb,mu,q, Efc, E, nu, aX0
      REAL*8 ak,ag,aInf,E11,E22,G12,G23,Kf,Km,Kt,v12,v23,v21
      REAL*8 c(NTENS,NTENS), DSTRESS(NTENS) 
      integer k,j
C
      !Material definition
      Em0=PROPS(1)     ! Resin Uncured E
      EmI=PROPS(2)     ! Resin Fully cured E
      ag=PROPS(3)      ! Gel Point / Load transfer Point
      aInf=PROPS(4)    ! DoC achieving the final properties i.e. EmI
      num=PROPS(5)     ! Resin Poissons ratio (assumed constant)
      aHg=PROPS(6)     ! Resin CTE heat-up - glassy
      aCg=PROPS(7)     ! Resin CTE cooldown - glassy
      alphaf=PROPS(8)  ! Fibre CTE - Isotropic
      Ef=PROPS(9)      ! Fibre E
      G12f=PROPS(10)   ! Fibre G12
      G23f=PROPS(11)   ! Fibre G23
      Vf=PROPS(12)     ! Fibre volume fraction
      nuf=PROPS(13)    ! Fibre Poissons ratio
      VgInf = PROPS(14)! Total load-transferring shrinkage
      aX0 = PROPS(15)
      aX1 = PROPS(16)
      aX2 = PROPS(17)
      a1 = PROPS(18)
      a2 = PROPS(19)
      a3 = PROPS(20)
      a4 = PROPS(21)
      Tp1 = PROPS(22)
      Tp2 = PROPS(23)
      Tp3 = PROPS(24)
      Tp4 = PROPS(25)
      Tp5 = PROPS(26)
C
      ! Define basic stiffness parameters and cure deltas
      ak=STATEV(1)
      Tg = STATEV(2)
      am=(ak-ag)/(aInf-ag)
      IF (ak<ag) THEN
        Em=Em0
      ELSEIF (ak.GE.ag) THEN
        Em=Em0*(1.0-am)+EmI*am
      ELSEIF (ak.GE.ag.AND.TEMP(1).GE.Tg) THEN
        Em=Em0*(1.0-am)+0.01*EmI*am
      ELSEIF (ak.GE.aInf.AND.TEMP(1).LT.Tg) THEN
        Em=EmI
      ELSE
        Em = EmI*0.01
      ENDIF
      Gm=0.5*Em/(1.0+num)
      q = -10 ! = -100 - hard, = -10 - very soft
C
      IF (STRAN(1)<0.0 .or. ak<ag) THEN
        Ef = (0.01 + 0.99/(1.0+EXP(q*(ak-ag))))*Ef
        G12f = (0.01 + 0.99/(1.0+EXP(q*(ak-ag))))*G12f
        G23f = (0.01 + 0.99/(1.0+EXP(q*(ak-ag))))*G23f
      ENDIF
C
      ! Determine terms using https://doi.org/10.3390/polym14142903
      Kf=Ef/(2.0*(1.0-nuf-2.0*(nuf**2.0)))
      Km=Em/(2.0*(1.0-num-2.0*(num**2.0)))
      Kt=((Kf+Gm)*Km+(Kf-Km)*Gm*Vf)/(Kf+Gm-(Kf-Km)*Vf)
      G12=Gm*(G12f+Gm+(G12f-Gm)*Vf)/(G12f+Gm-(G12f-Gm)*Vf)
      G23=Gm*(Km*(Gm+G23f)+2.0*Gm*G23f+Km*(G23f-Gm)*Vf)/(Km*(Gm+G23f)
     1 +2.0*Gm*G23f-(Km+2.0*Gm)*(G23f-Gm)*Vf)
      E11=Ef*Vf+Em*(1.0-Vf)+((4.0*(num-nuf**2.0)*Kf*Km*Gm*(1.0-Vf)*Vf)
     1 /((Kf+Gm)*Km+(Kf-Km)*Gm*Vf))
      v12=nuf*Vf+num*(1.0-Vf)+((num-nuf)*(Km-Kf)*Gm*(1.0-Vf)*Vf)
     1 /((Kf+Gm)*Km+(Kf-Km)*Gm*Vf)
      E22=1.0/(0.25/Kt+0.25/G23+(v12**2.0)/E11)
      v23=(2.0*E11*Kt-E11*E22-4.0*(v12**2.0)*Kt*E22)/(2.0*E11*Kt)
      v21=v12*E22/E11
      !Definition of matrix of elasticity (Orthotropic)
!     JACOBIAN MATRIX
      IF (NTENS.EQ.1) THEN ! TRUSS element
        c(1,1) = E11
      ELSEIF (NTENS.EQ.4) THEN !Plane Strain
C
        c(1,1)=E11*(1.0-v23**2.0)/((1.0-v23-v12*v21)*(1.0+v23))
        c(1,2)=E11*(v21+v21*v23)/((1.0-v23-v12*v21)*(1.0+v23))
        c(1,3)=E11*(v21+v21*v23)/((1.0-v23-v12*v21)*(1.0+v23))
        c(1,4)=0.0d0
C
        c(2,1)=E22*(v12+v12*v23)/((1.0-v23-v12*v21)*(1.0+v23))
        c(2,2)=E22*(1.0-v12*v21)/((1.0-v23-v12*v21)*(1.0+v23))
        c(2,3)=E22*(v23+v12*v21)/((1.0-v23-v12*v21)*(1.0+v23))
        c(2,4)=0.0d0
C
        c(3,1)=E22*(v12+v12*v23)/((1.0-v23-v12*v21)*(1.0+v23))
        c(3,2)=E22*(v23+v12*v21)/((1.0-v23-v12*v21)*(1.0+v23))
        c(3,3)=E22*(1.0-v12*v21)/((1.0-v23-v12*v21)*(1.0+v23))
        c(3,4)=0.0d0
C
        c(4,1)=0.0d0
        c(4,2)=0.0d0
        c(4,3)=0.0d0
        c(4,4)=G12
C
      ELSEIF (NTENS.EQ.3) THEN !Plane Stress
        ! WRITE(6,*)'Using plane stress orthotropic jacobian'
        c(1,1)=E11/(1.0-v12**2.0)
        c(1,2)=v12*E11/(1.0-v12**2.0)
        c(1,3)=0.0d0
C
        c(2,1)=v12*E22/(1.0-v12**2.0)
        c(2,2)=E22/(1.0-v12**2.0)
        c(2,3)=0.0d0
C
        c(3,1)=0.0d0
        c(3,2)=0.0d0
        c(3,3)=G12
C
      ELSEIF(NTENS.EQ.6) THEN !3D
        ! WRITE(6,*)'Using 3D orthotropic jacobian'
        c(1,1)=E11*(1.0-v23**2.0)/((1.0-v23-v12*v21)*(1.0+v23))
        c(1,2)=E11*(v21+v21*v23)/((1.0-v23-v12*v21)*(1.0+v23))
        c(1,3)=E11*(v21+v21*v23)/((1.0-v23-v12*v21)*(1.0+v23))
        c(1,4)=0.0d0
        c(1,5)=0.0d0
        c(1,6)=0.0d0
C
        c(2,1)=E22*(v12+v12*v23)/((1.0-v23-v12*v21)*(1.0+v23))
        c(2,2)=E22*(1.0-v12*v21)/((1.0-v23-v12*v21)*(1.0+v23))
        c(2,3)=E22*(v23+v12*v21)/((1.0-v23-v12*v21)*(1.0+v23))
        c(2,4)=0.0d0
        c(2,5)=0.0d0
        c(2,6)=0.0d0
C
        c(3,1)=E22*(v12+v12*v23)/((1.0-v23-v12*v21)*(1.0+v23))
        c(3,2)=E22*(v23+v12*v21)/((1.0-v23-v12*v21)*(1.0+v23))
        c(3,3)=E22*(1.0-v12*v21)/((1.0-v23-v12*v21)*(1.0+v23))
        c(3,4)=0.0d0
        c(3,5)=0.0d0
        c(3,6)=0.0d0
C
        c(4,1)=0.0d0
        c(4,2)=0.0d0
        c(4,3)=0.0d0
        c(4,4)=G23
        c(4,5)=0.0d0
        c(4,6)=0.0d0
C
        c(5,1)=0.0d0
        c(5,2)=0.0d0
        c(5,3)=0.0d0
        c(5,4)=0.0d0
        c(5,5)=G12
        c(5,6)=0.0d0
C
        c(6,1)=0.0d0
        c(6,2)=0.0d0
        c(6,3)=0.0d0
        c(6,4)=0.0d0
        c(6,5)=0.0d0
        c(6,6)=G12
C
      ELSE
        WRITE(6,*)'Case with ntens =',NTENS,' not implemented'
      ENDIF
C     Calculate stress increment
      DO k=1,NTENS
        DSTRESS(k)=0.d0
        DO i=1,NTENS
           DSTRESS(k) = DSTRESS(k) + c(k,i)*DSTRAN(i)
        ENDDO 
      ENDDO
C
C     Update stress components
      DO k=1,NTENS
        STRESS(k)=STRESS(k)+DSTRESS(k)
      ENDDO
C
C     Determine Jacobian
      DO k=1,NTENS
        DO j=1,NTENS
           DDSDDE(k,j)=c(k,j)
        ENDDO
      ENDDO
      ! Save Properties as STATEV
      STATEV(9) = E11       ! Young's Modulus 11
      STATEV(10) = E22      ! Young's Modulus 22
      STATEV(12) = ag       ! DoC at Load transfer point
      STATEV(13) = aInf     ! Assumed DoC at shrinkage stop
      STATEV(14) = VgInf    ! Total load-transferring shrinkage
      STATEV(15) = aHg      ! Resin CTE heat-up - glassy
      STATEV(16) = aCg      ! Resin CTE cooldown - glassy
      STATEV(17) = alphaf   ! Fibre CTE - Isotropic!
      STATEV(18) = Ef       ! Fibre Young's modulus - Isotropic!
      STATEV(19) = nuf      ! Fibre Poisson's ratio - isotropic
      STATEV(20) = num      ! Resin Poisson's ratio - Isotropic!
      STATEV(21) = Em       ! Resin Young's modulus - func of DoC
      STATEV(22) = Vf       ! FVF - IF==0.0 THEN neat resin, IF==1.0 pure fibre
      STATEV(24) = aX0
      STATEV(25) = aX1
      STATEV(26) = aX2
      STATEV(27) = a1
      STATEV(28) = a2
      STATEV(29) = a3
      STATEV(30) = a4
      STATEV(31) = Tp1
      STATEV(32) = Tp2
      STATEV(33) = Tp3
      STATEV(34) = Tp4
      STATEV(35) = Tp5
      RETURN
      END
C
      SUBROUTINE UEXPAN(EXPAN,DEXPANDT,TEMP,TIME,DTIME,PREDEF,
     1 DPRED,STATEV,CMNAME,NSTATV,NOEL)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION EXPAN(*),DEXPANDT(*),TEMP(2),TIME(2),PREDEF(*),
     1 DPRED(*),STATEV(NSTATV)
C
      ! Declare dimensions for all variables to be used.
      REAL*8 ak,curerate,Tg,VgInf,ag,aInf,DVg,am
      REAL*8 dechm,dech1,dech2,dech3,deth1,deth2,deth3
      REAL*8 alpham,alpha1,alpha2,alpha3,alphaf,Vf,Em,num,nuf,Ef
      REAL*8 Ts,TC3,TC2,TC1,TH3,TH2,TH1,Tp1,Tp2,Tp3,Tp4,Tp5
      REAL*8 aCHT,aC2,aC1,aCg,aHHT,aH4,aH3,aH2,aH1,aHg,aX2,aX1,aX0
      REAL*8 a1,a2,a3,a4
C
      ! Variables needed for calculation
      ak = STATEV(1)        ! Degree of Cure
      Tg = STATEV(2)        ! Glass transition temperature
C
      !For "internal" procedures only!
      curerate = STATEV(11) ! Rate of Cure
      ag = STATEV(12)       ! DoC at Load transfer point
      aInf = STATEV(13)     ! Assumed DoC at shrinkage stop
      VgInf = STATEV(14)    ! Total load-transferring shrinkage
      aHg = STATEV(15)      ! Resin CTE heat-up - glassy
      aCg = STATEV(16)      ! Resin CTE cooldown - glassy
      alphaf = STATEV(17)   ! Fibre CTE - Isotropic!
      Ef = STATEV(18)       ! Fibre Young's Modulus - Isotropic
      nuf = STATEV(19)      ! Fibre Poisson's ratio - Isotropic
      num = STATEV(20)      ! Resin Poisson's ratio - Isotropic
      Em = STATEV(21)       ! Resin Young's modulus - Fun of DoC
      Vf = STATEV(22)       ! FVF - IF == 0.0 THEN neat resin
      aX0 = STATEV(24)      !
      aX1 = STATEV(25)
      aX2 = STATEV(26)
      a1 = STATEV(27)
      a2 = STATEV(28)
      a3 = STATEV(29)
      a4 = STATEV(30)
      Tp1 = STATEV(31)
      Tp2 = STATEV(32)
      Tp3 = STATEV(33)
      Tp4 = STATEV(34)
      Tp5 = STATEV(35)
      Ash = 0.0d0
C
      Ts = TEMP(1)-Tg
      ! aC3 = 2.76702108e-06
      ! aC2 = 9.49639559e-07
      ! aC1 = 5.18111403e-07
C
      aHHT = aX2*ak**2+aX1*ak+aX0 ! Only valid from ak=>ag
      aCHT = aHHT
      aH4 = a4
      aH3 = a3
      aH2 = a2
      aH1 = a1
C
      TH1 = Tp1
      TH2 = Tp2
      TH3 = Tp3
      TH4 = Tp4
      TH5 = Tp5
C
      ! TC1 = -56.
      ! TC2 = -29.
      ! TC3 = -16.
      ! TC4 = -3.
      TC1 = TH1
      TC2 = TH2
      TC3 = TH3
      TC4 = TH4
      aC1 = aH1
      aC2 = aH2
      aC3 = aH2
C
      IF (Ef<1.0d0) THEN
        Ef = 100000.0d0
      ENDIF
      IF (Em<1.0d0) THEN
        Em = 100000.0d0
      ENDIF
      dTdt = TEMP(2)/DTIME
C
      ! Define resin CTE w. Tg
      IF (ak.LT.ag) THEN
        alpham=0.0d0
      ELSEIF (ak.GE.ag) THEN
        IF (dTdt.GE.0.0001) THEN
          IF (Ts.GE.TH5) THEN
            alpham = aHHT
          ELSEIF (Ts.GE.TH4.AND.Ts.LT.TH5) THEN
            alpham=aH4*(Ts-TH4)+aH3*(TH4-TH3)+aH2*(TH3-TH2)+aH1*(TH2-TH1)+aHg
          ELSEIF (Ts.GE.TH3.AND.Ts.LT.TH4) THEN
            alpham = aH3*(Ts-TH3)+aH2*(TH3-TH2)+aH1*(TH2-TH1)+aHg
          ELSEIF (Ts.GE.TH2.AND.Ts.LT.TH3) THEN
            alpham = aH2*(Ts-TH2)+aH1*(TH2-TH1)+aHg
          ELSEIF (Ts.GE.TH1.AND.Ts.LT.TH2) THEN
            alpham = aH1*(Ts-TH1)+aHg
          ELSEIF (Ts.LT.TH1) THEN
            alpham = aHg
          ENDIF
        ELSEIF(dTdt.LT.0.0001) THEN
          IF (Ts.GE.TC4) THEN
            alpham = aCHT
          ELSEIF (Ts.GE.TC3.AND.Ts.LT.TC4) THEN
            alpham = aC3*(Ts-TC3)+aC2*(TC3-TC2)+aC1*(TC2-TC1)+aCg
          ELSEIF (Ts.GE.TC2.AND.Ts.LT.TC3) THEN
            alpham = aC2*(Ts-TC2)+aC1*(TC2-TC1)+aCg
          ELSEIF (Ts.GE.TC1.AND.Ts.LT.TC2) THEN
            alpham = aC1*(Ts-TC1)+aCg
          ELSEIF (Ts.LT.TC1) THEN
            alpham = aCg 
          ENDIF
        ELSE
          alpham = 0.0d0
        ENDIF
      ELSE
        alpham = 0.0d0
      ENDIF
C
C     !Compute composite expansions
      alpha1 = (alphaf*Ef*Vf+alpham*Em*(1.0d0-Vf))/(Ef*Vf+Em*(1.0-Vf))
      alpha2 = (alphaf+nuf*alphaf)*Vf+(alpham+num*alpham)*(1.0-Vf) 
     1-(nuf*Vf+num*(1.0-Vf))*alpha1
      alpha3 = alpha2
C
C   !Compute composite expansion components
      deth1 = alpha1*TEMP(2) ! Thermal strains 1,2,3-direction
      deth2 = alpha2*TEMP(2)
      deth3 = alpha3*TEMP(2)
C
C
      ! Find change in cure Shrinkage
      IF (ak<ag) THEN
        DVg=0.0d0
      ELSEIF (ak>ag) THEN
        DVg=curerate*DTIME*(Ash/(aInf-ag)
     1+(2.0d0*(VgInf-Ash)*(ak-ag))/((aInf-ag)**2))
      ELSE ! Make up for any remaining shrinkage
        DVg=VgInf-(statev(6)+1.0d0)**(3)+1.0d0
      ENDIF
C
      ! Compute composite shrinkage components
      dechm = ((1.0d0+DVg)**(1.0d0/3.0d0))-1.0d0
      ! Compute the chemical shrinkage in material directions
      IF (Vf.GT.0.95) THEN
        dech1 = 0.0d0
        dech2 = 0.0d0
        dech3 = 0.0d0
      ELSE
        dech1 = (dechm*Em*(1.0-Vf))/(Ef*Vf+Em*(1.0-Vf))
        dech2 = (dechm+num*dechm)*(1.0-Vf)
     1  -(nuf*Vf+num*(1.0-Vf))*dech1
        dech3 = dech2 ! For Transverse Isotropic conditions
      ENDIF
C
      !Set thermal expansion and derivative w.r.t. temperature
      !Isotropic
      !EXPAN(1)=deth1+dech1
      !DEXPANDT(1)=alpha1
      !Orthotropic
      EXPAN(1) = deth1+dech1 !e11
      EXPAN(2) = deth2+dech2 !e22
      ! EXPAN(3)=deth3+dech3 !e33
      EXPAN(3) = 0.0
      EXPAN(4) = 0.0 !e12
      DEXPANDT(1) = alpha1
      DEXPANDT(2) = alpha2
      ! DEXPANDT(3)=alpha3
      DEXPANDT(3) = 0.0
      DEXPANDT(4) = 0.0
      !Anisotropic
      !EXPAN(1)=deth1+dech1 !e11
      !EXPAN(2)=deth2+dech2 !e22
      !EXPAN(3)=deth3+dech3 !e33
      !EXPAN(4)=0.0 !e12
      !EXPAN(5)=0.0 !e13
      !EXPAN(6)=0.0 !e23
      !DEXPANDT(1)=alpha1
      !DEXPANDT(2)=alpha2
      !DEXPANDT(3)=alpha3
      !DEXPANDT(4)=0.0
      !DEXPANDT(5)=0.0
      !DEXPANDT(6)=0.0
C
      ! Store strain components
      !Halved as the subroutine is called twice per iteration
      !See:https://abaqus-docs.mit.edu/2017/English/
      !SIMACAESUBRefMap/simasub-c-uexpan.htm#simasub-c-uexpan
      STATEV(3) = STATEV(3)+0.5d0*deth1     ! Thermal strain 11
      STATEV(4) = STATEV(4)+0.5d0*deth2     ! Thermal strain 22
      STATEV(5) = STATEV(5)+0.5d0*deth3     ! Thermal strain 33
      STATEV(6) = STATEV(6)+0.5d0*dech1     ! Chemical strain 11
      STATEV(7) = STATEV(7)+0.5d0*dech2     ! Chemical strain 22
      STATEV(8) = STATEV(8)+0.5d0*dech3     ! Chemical strain 33
      STATEV(36) = Ts
      STATEV(37) = dTdt
      RETURN
      END
      ! FILM is supplimentary and used for custom heat transfer.
      SUBROUTINE FILM(H,SINK,TEMP,KSTEP,KINC,TIME,NOEL,NPT,
     1 COORDS,JLTYP,FIELD,NFIELD,SNAME,NODE,AREA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION H(2),TIME(2),COORDS(3),FIELD(NFIELD)
      CHARACTER*80 SNAME
C
      real SINKARRAY(6,2)
      integer ii
C
      SINKARRAY(1,1)=0.0 ! Initiate @ 20 C
      SINKARRAY(1,2)=20.0
      SINKARRAY(2,1)=1200.0 ! Heat-up
      SINKARRAY(2,2)=40.0
      SINKARRAY(3,1)=19200.0 ! 5 hr @ 40 C
      SINKARRAY(3,2)=40.0
      SINKARRAY(4,1)=21600.0 ! Heat-up
      SINKARRAY(4,2)=80.0
      SINKARRAY(5,1)=39600.0 ! 5hr @ 80 C
      SINKARRAY(5,2)=80.0
      SINKARRAY(6,1) = 43200 ! Cool-down
      SINKARRAY(6,2) = 20
C
      ! Find current sink temperature
      SINK=20.0
      do ii = 1,5
        if (TIME(2)==SINKARRAY(ii,1)) then
          SINK=SINKARRAY(ii,2)
        else if (TIME(2)>SINKARRAY(ii,1)) then
          SINK=SINKARRAY(ii,2)+(SINKARRAY(ii+1,2)-
     1 SINKARRAY(ii,2))*(TIME(2)-SINKARRAY(ii,1))/
     2 (SINKARRAY(ii+1,1)-SINKARRAY(ii,1))
        end if
      end do
C
      ! Set film coefficient temperature
      if (SINK>TEMP) then
        H(1)=500.0
        H(2)=0.0
      else
        H(1)=20.0
        H(2)=0.0
      end if
C
      RETURN
      END
C
C      ! SUBROUTINE FRIC(LM,TAU,DDTDDG,DDTDDP,DSLIP,SED,SFD,
C     ! 1 DDTDDT,PNEWDT,STATEV,DGAM,TAULM,PRESS,DPRESS,DDPDDH,SLIP,
C     ! 2 KSTEP,KINC,TIME,DTIME,NOEL,CINAME,SLNAME,MSNAME,NPT,NODE,
C     ! 3 NPATCH,COORDS,RCOORD,DROT,TEMP,PREDEF,NFDIR,MCRD,NPRED,
C     ! 4 NSTATV,CHRLNGTH,PROPS,NPROPS)
C
C      ! INCLUDE 'ABA_PARAM.INC'
C
C      ! CHARACTER*80 CINAME,SLNAME,MSNAME
C
C      ! DIMENSION TAU(NFDIR),DDTDDG(NFDIR,NFDIR),DDTDDP(NFDIR),
C     ! 1 DSLIP(NFDIR),DDTDDT(NFDIR,2),STATEV(*),DGAM(NFDIR),
C     ! 2 TAULM(NFDIR),SLIP(NFDIR),TIME(2),COORDS(MCRD),
C     ! 3 RCOORD(MCRD),DROT(2,2),TEMP(2),PREDEF(2,*),PROPS(NPROPS)
C      ! REAL*8 Gts, Gts0,GtsI, Em, ut, u1p, hs, num 
C     
C      ! hs = PROPS(1)
C      ! num = STATEV(20)
C      ! Em = STATEV(21)
C      ! Gts=0.5*Em/(1.0+num)
C
C      ! IF (ak<ag) THEN
C          ! LM = 2
C      ! ELSEIF (ak.GE.ag) THEN
C         ! LM = 0
C      ! ENDIF
C
C      ! TAU = Gts*(ut-u1p)/(hs) 
C      ! user coding to define LM, TAU, DDTDDG, DDTDDP,
C      ! and, optionally, DSLIP, SED, SFD, DDTDDT, PNEWDT, STATEV
C      ! RETURN
C      ! END