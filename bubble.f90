module parameters
  IMPLICIT NONE

  !-- common parameters and variables ------------------------------
  ! THIS IS ALMOST PROJECT-INDEPENDENT 
  double precision, parameter :: tm32   = 1.d0/(2.d0**32.d0)
  double precision, parameter :: eps    = 1.d-14            ! very small number
  double precision, parameter :: tol    = 0.15d0            ! tolerance for Cor
  logical                     :: prt                        ! flag for write2file
  integer,          parameter :: Mxint  = 2147483647        ! maximum integer
  integer,          parameter :: Mnint  =-2147483647        ! minimum integer
  double precision, parameter :: pi=3.1415926

  !-- Parameters -------------------------------------------------
  integer, parameter :: D=2  
  integer, parameter :: UP=1
  integer, parameter :: DOWN=0
  integer, parameter :: MxL=512     !Max size of the system
  integer          :: PID      ! the ID of this job
  integer :: L     !the actual size of the system
  integer          :: Order
  double precision :: Beta
  double precision :: Mu
  double precision :: rs
  double precision :: EF, kF
  double precision :: Mass2
  double precision :: TotalStep  !total steps of this MC simulation
  double precision :: Step    ! a counter to keep track of the current step number
  integer, parameter :: UpdateNum=2    ! number of updates
  integer                      :: Seed                   ! random-number seed
  double precision, dimension(UpdateNum) :: PropStep
  double precision, dimension(UpdateNum) :: AcceptStep

  !-- Diagram Permutation Table ----------------------------------
  integer, parameter :: PHYSICAL=1
  integer, parameter :: NORMALIZATION=2
  integer, parameter :: MaxOrder=8 ! Max diagram order
  integer, parameter :: MaxDiagNum=1000 ! Max diagram number 
  integer, parameter :: MaxIndepdentG=10000 ! Max indepdent Green function number 
  integer, parameter :: MaxIndepdentVer=10000 ! Max indepdent vertex number
  integer, parameter :: MaxEK=10000 ! Max indepdent vertex number
  integer, parameter :: MaxTau=10000 ! Max indepdent vertex number
  integer :: GNum, VerNum, LoopNum, DiagNum, TauNum, HugenholtzNum !Number of G, Vertex, Loop, diagram, 
  integer :: iGNum, iVerNum !Number of independent G and Vertex
  integer :: numDiagV, DiagNum1H

  !double precision, dimension(0:MaxEK, 0:MaxTau) :: GreenTable  !propagator for different k and tau
  ! all diagram-related arraies have two copies: 1 is for physical diagrams, 2 is normalization diagrams
  integer, allocatable :: GIndex(:,:), GIndexNorm(:,:)  ! Index of Green function
  integer, allocatable :: VerIndex(:,:), VerIndexNorm(:,:) ! Index of vertex
  double precision,allocatable ::SymFactor(:), SymFactorNorm(:)  ! Symmetry Factor (includes diagram sign)
  double precision,allocatable ::SpinFactor(:,:), SpinFactorNorm(:,:)  ! Spin Factor (includes diagram sign)
  double precision,allocatable ::SpinCache(:,:), SpinCacheNorm(:,:) 
  integer, allocatable :: LoopBases(:,:), LoopBasesNorm(:,:) ! Bases for loops
  integer, allocatable :: LoopBasesVer(:,:), LoopBasesVerNorm(:,:) ! Bases for loops including vertex
  integer, allocatable :: TauBases(:,:), TauBasesNorm(:,:) ! Permutation 
  double precision, allocatable :: DiagWeight(:), DiagWeightABS(:)
  double precision, allocatable :: GWeight(:), GWeightNorm(:) !Weight of green function
  double precision, allocatable :: VerWeight(:), VerWeightNorm(:) ! Weight of vertex (interation)
  double precision, allocatable :: NewGWeight(:), NewGWeightNorm(:) !Weight of green function
  double precision, allocatable :: NewVerWeight(:), NewVerWeightNorm(:) ! Weight of vertex (interation)
  double precision, allocatable :: LoopMom(:,:), LoopMomNorm(:,:) ! values to attach to each loop basis
  double precision, allocatable:: LoopSpin(:), LoopSpinNorm(:) ! values to attach to each spin
  double precision, allocatable :: TauTable(:), TauTableNorm(:) ! Time table for each Tau, all Tau are between [0,beta)


  !integer, dimension(MaxDiagNum, 2*MaxOrder) :: GIndex, GIndexNorm  ! Index of Green function
  !integer, dimension(MaxDiagNum, 2*MaxOrder) :: VerIndex,  VerIndexNorm ! Index of vertex
  !double precision,dimension(MaxDiagNum) ::SymFactor, SymFactorNorm  ! Symmetry Factor (includes diagram sign)
  !integer, dimension(MaxOrder+1, MaxIndepdentG) :: LoopBases, LoopBasesNorm ! Bases for loops
  !integer, dimension(MaxOrder+1, MaxIndepdentVer) :: LoopBasesVer, LoopBasesVerNorm ! Bases for loops including vertex
  !integer, dimension(2*MaxOrder, MaxIndepdentG) :: TauBases, TauBasesNorm ! Permutation 
  !double precision, dimension(MaxDiagNum) :: DiagWeight, DiagWeightABS

  !double precision, dimension(MaxIndepdentG) :: GWeight, GWeightNorm !Weight of green function
  !double precision, dimension(MaxIndepdentVer) :: VerWeight, VerWeightNorm ! Weight of vertex (interation)

  !double precision, dimension(MaxIndepdentG) :: NewGWeight, NewGWeightNorm !Weight of green function
  !double precision, dimension(MaxIndepdentVer) :: NewVerWeight, NewVerWeightNorm ! Weight of vertex (interation)

  !!double precision, dimension(MaxIndepdentG) :: FreqTable, FreqTableNorm  !Freq variable for each propagator
  !!integer, dimension(D, MaxIndepdentG) :: MomTable, MomTableNorm  !Momentum variable for each propagator, kx, ky, kz
  !!integer, dimension(MaxIndepdentG) :: SpinTable, SpinTableNorm  !spin variable for each vertex function, 1: spin for in-leg, 2:spin for out-leg

  !double precision, dimension(D, MaxOrder+1) :: LoopMom, LoopMomNorm ! values to attach to each loop basis
  !integer, dimension(MaxOrder+1) :: LoopSpin, LoopSpinNorm ! values to attach to each spin
  !double precision, dimension(2*MaxOrder) :: TauTable, TauTableNorm ! Time table for each Tau, all Tau are between [0,beta)

  integer :: Offset, OffVer, OffDiag


  double precision                   :: OldFreq, ExtFreq
  double precision                   :: OldTau, ExtTau
  double precision, dimension(D)     :: OldMom
  integer                            :: OldSpin, ExtSpin
  double precision                            :: OldWeight
  integer                            :: ExtMomBin

  integer                           :: Sector  !a flag, 1: physical diagram, 2:  normalization diagram
  !the weight of all diagrams in current sector=AbsWeight*Phase

  !-- Measurement  ------------------------------------------------
  double precision :: NormWeight    !the accumulated weight of the normalization diagram
  double precision :: normConstant
  integer, parameter          :: QBinNum=32     !number of q bins
  double precision            :: ExtMomMax
  double precision            :: DeltaQ
  !integer, parameter          :: WBinNum=128     !number of q bins
  double precision, dimension(20) :: F1, F2, F3

  double precision, dimension(QBinNum) :: Polarization  !the accumulated weight of the Spin-zz polarization
  double precision, dimension(MaxDiagNum, QBinNum) :: PolarDiag  !the accumulated weight of the Spin-zz polarization
  double precision, dimension(QBinNum) :: TauPolarization  !the accumulated weight of the Spin-zz polarization
  double precision, dimension(D, QBinNum)   :: ExtMomMesh
  
  logical, dimension(MaxIndepdentVer) :: reducibleTable

end module

INCLUDE "rng.f90"

program main
    use mt19937
    use parameters
    implicit none
    double precision :: x
    integer :: PrintCounter, SaveCounter
  
    print *, 'Beta, rs, Mass2, Order, HugenholtzNum, TotalStep(*1e6), Seed, PID'
    read(*,*)  Beta, rs, Mass2, Order, HugenholtzNum, TotalStep, Seed, PID
    !Beta=10.0
    !rs=1.0
    !Mass2=1.0
    !kF=(9.0*pi/4.0)**(1.0/3.0)/rs !3D
    kF=sqrt(2.0)/rs !2D
    EF=kF*kF
    Mu=EF
    print *,"Inverse Temperature:", Beta
    print *,"rs:", rs
    print *,"Fermi Mom:", kF
    print *,"Fermi Energy:", EF
  
    call sgrnd(Seed) 
  
    print *, "Initializing ..."
    call Initialize()
    print *, "Initialized!"
  !  call ForTest
  !  stop "stop here"
  
    call Test() !call test first to make sure all subroutines work as expected
  
    TotalStep=TotalStep*1.0e6
    print *, "Start simulation..."
    do while (Step<TotalStep)
        Step=Step+1.0
        x=grnd()
        if (x<1.0/UpdateNum) then
          call ChangeTau()
        else if (x<2.0/UpdateNum) then
          call ChangeMom()
        else if (x<3.0/UpdateNum) then
          call NormToPhysical()
        else if (x<4.0/UpdateNum) then
          call PhysicalToNorm()
        !else if (x<5.0/UpdateNum) then
          !!call ChangeSpin()
          !call SwapMom()
        endif
        !if(mod(int(Step), 4)==0) call Measure()
        call Measure()
  
        !call DynamicTest()
  
        PrintCounter=PrintCounter+1
        if (PrintCounter==1e7)  then
          write(*,*) 
          write(*,"(f8.2, A15)") Step/1e6, "million steps"
          write(*,"(A14, A6, f8.3, A5, f8.3)") "Accept Ratio: ", "Freq:", AcceptStep(1)/PropStep(1),&
             & "Mom:", AcceptStep(2)/PropStep(2)
  !        write(*,"(A14, A6, f8.3)") "Accept Ratio: ", "Swap:", AcceptStep(3)/PropStep(3)
          PrintCounter=0
  
          print *, "mom1", norm2(LoopMom(:,1)), LoopMom(:,1)
          print *, "mom2", norm2(LoopMom(:,2)), LoopMom(:,2)
          print *, "Momnorm", norm2(LoopMom(:,2)+LoopMom(:,1))**2-Mu, norm2(LoopMom(:,2))**2-Mu
          !print *, "mom3", norm2(LoopMom(:,3))
          !print *, "mom4", norm2(LoopMom(:,4))
          print *, "tau table", TauTable(1:TauNum)
          print *, "green", GWeight(1), GWeight(2)
          print *, "weight", OldWeight
        endif
        if (PrintCounter==1e8)  then
          call SaveToDisk()
          call SaveToDiskF()
        endif
    end do

    call SaveToDisk()
    call SaveToDiskF()
    
    print *, "End simulation."
  
    CONTAINS


!    subroutine ForTest
!        implicit none
!        integer :: QNum = 96
!        double precision :: ddeltaQ
!
!        TauTable(1:2) = [ 0.0, 1.2, 1.5, 1.8, 2.0, 2.3 ]
!        ddeltaQ = 3.0*kF/QNum 
!
!        LoopMom(:,1) = [  1.5, 0.0  ]
!        LoopMom(:,2) = [  1.0, 2.0  ]
!
!        call ResetWeightTable()
!
!        print *, "momentum=", LoopMom(:,1), LoopMom(:,2)
!        print *, "GWeight=", GWeight(1:6)
!        print *, "GWeight=", VerWeight(1:6)
!
!        OldWeight = CalcWeight(0)
!
!        print *, "OldWeight=", OldWeight 
!    end subroutine


  subroutine Test()
    implicit none
    integer :: i
    double precision :: ratio, old, new
    ! Put all tests here
    !if (cabs(Green(0.d0, -1.d0)+Green(0.d0, Beta-1.d0))>1e-6) then
      !print *, "Green's function is not anti-periodic"
      !stop
    !endif

    old=0.0
    !call GenerateNewFreq(old, new, ratio)
    print *, old, new, ratio
  end subroutine

  subroutine Initialize()
    implicit none
    integer :: i, num

    GNum = 2*Order   	
    VerNum = Order-1 	
    LoopNum = Order+1	
    TauNum = 2*Order	
    DiagNum1H = 2**VerNum

    DiagNum = DiagNum1H * HugenholtzNum  
    iGNum = GNum * HugenholtzNum
    iVerNum = 2 * VerNum  * HugenholtzNum 
    
    
    allocate(GIndex(HugenholtzNum, 2*Order))
    allocate(GIndexNorm(HugenholtzNum, 2*Order))

    allocate(VerIndex(HugenholtzNum, 2*Order), VerIndexNorm(HugenholtzNum, 2*Order)) ! Index of vertex
    allocate(SymFactor(HugenholtzNum), SymFactorNorm(HugenholtzNum))  ! Symmetry Factor (includes diagram sign)
    allocate(SpinFactor(2**VerNum, HugenholtzNum), SpinFactorNorm(2**VerNum, HugenholtzNum))  ! Symmetry Factor (includes diagram sign)
    allocate(SpinCache(2**VerNum, Order-1), SpinCacheNorm(2**VerNum, Order-1)) 

    allocate(LoopBases(Order+1, iGNum), LoopBasesNorm(Order+1, iGNum)) ! Bases for loops
    allocate(LoopBasesVer(Order+1, iVerNum), LoopBasesVerNorm(Order+1, iVerNum)) ! Bases for loops including vertex
    allocate(TauBases(2*Order, iGNum), TauBasesNorm(2*Order, iGNum)) ! Permutation 

    allocate(DiagWeight(HugenholtzNum), DiagWeightABS(HugenholtzNum))
    allocate(GWeight(iGNum), GWeightNorm(iGNum)) !Weight of green function
    allocate(NewGWeight(iGNum), NewGWeightNorm(iGNum)) !Weight of green function
    allocate(VerWeight(iVerNum), VerWeightNorm(iVerNum))
    allocate(NewVerWeight(iVerNum), NewVerWeightNorm(iVerNum))
    allocate(LoopMom(D, Order+1), LoopMomNorm(D, Order+1)) ! values to attach to each loop basis
    allocate(LoopSpin(Order+1), LoopSpinNorm(Order+1)) ! values to attach to each spin
    allocate(TauTable(2*Order), TauTableNorm(2*Order)) ! Time table for each Tau, all Tau are between [0,beta)

    print *, "Reading Diagram ..."
    call ReadDiagram()
    print *, "Read Diagram done!"
    call Reducible
    !call CreateGTable()

    normConstant = 0.3

    PropStep=0.0
    AcceptStep=0.0

    ExtMomMax = kF
    DeltaQ=ExtMomMax/QBinNum

    Polarization=0.0
    PolarDiag=0.0
    F1 = 0.0
    F2 = 0.0
    F3 = 0.0

    ExtMomMesh=0.0
    do i=1, QBinNum
      ExtMomMesh(1, i)=(i-0.5)*DeltaQ
    enddo

    Step=0.0
    PrintCounter=0
    SaveCounter=0

    ExtMomBin=1

    do i=0, Order-1
      num=2*i+1
      if(num==1) then
        TauTable(num)=0.0
        TauTable(num+1)=Beta/2.0
      else
        TauTable(num)=Beta*grnd()
        TauTable(num+1) = TauTable(num)
      endif
    enddo
    !TauTable(2)=-0.0001


    LoopMom(:,1)=ExtMomMesh(:, ExtMomBin)
    LoopMom(:,2:)=kF/sqrt(real(D))
    LoopSpin(1)=0
    LoopSpin(2)=1



    call ResetWeightTable()
    OldWeight = CalcWeight(0)
    call UpdateState()
    !call TestIni(OldWeight)
end subroutine


subroutine TestIni(weight)
    implicit none
    double precision :: weight, x

    if (weight>1.0e-20) then
        return
    else
        do while (weight<1.0e-20)
            x = grnd()
            if (x<1.0/UpdateNum) then
                call ChangeTau()
            else if (x<2.0/UpdateNum) then
                call ChangeMom()
            end if
            weight = OldWeight
        end do
    end if

end subroutine

    subroutine ReadDiagram()
        implicit none
        character(90) :: fname
        character (len=20) :: charc
        integer :: baseNum, i, numDiagV, num
        character( len = 3 ) :: Orderstr
    
        OffDiag = 0
        Offset = 0
        OffVer = 0
    
        write( Orderstr,'(I3)' )  Order 
    	fname = 'DiagPolar'//trim(adjustl(Orderstr))//'.txt'
    
        open(unit=10, file=trim(fname), action='read', status="old")
            do numDiagV=1, HugenholtzNum
                GIndex(numDiagV, 1:GNum) = (/ ((numDiagV-1)*GNum+i, i=1, GNum) /)
                VerIndex(numDiagV, 1:2*VerNum) = (/ ( (numDiagV-1)*2*VerNum + i, i=1, 2*VerNum ) /)
    
                Read(10, *) charc
                Read(10, *) TauBases(2, Offset+1:Offset+GNum)
                TauBases(1, Offset+1:Offset+GNum) = (/ (i, i=1,GNum) /) 
                TauBases(2, Offset+1:Offset+GNum) = TauBases(2, Offset+1:Offset+GNum) + 1
                Read(10, *) charc
                Read(10, *) SymFactor(numDiagV)
                Read(10, *) charc
                do baseNum=1, LoopNum
                    Read(10, *)  LoopBases(baseNum, Offset+1:Offset+GNum)
                end do
                Read(10, *) charc
                do baseNum=1, Order+1
                    Read(10, *)  LoopBasesVer(baseNum, OffVer+1:OffVer+2*VerNum)
                end do
                Read(10, *) charc
                Read(10, *) SpinFactor(:, numDiagV)
    
                OffDiag = OffDiag + DiagNum1H
                Offset = Offset + GNum
                OffVer = OffVer + 2*VerNum
            end do
    
        CLOSE(unit=10)
 
!        SymFactor(2:5) = 0.5 * SymFactor(2:5)

    end subroutine

    !double precision function Green(tau ,Mom, spin)
        !!calculate Green's function
        !implicit none
        !double precision :: tau, k2, s, Ek, x, y, w, r, coshv
        !integer :: spin, i
        !double precision, dimension(D) :: Mom
    !!    print *, tau, Mom, Spin
    !!    stop
        !s=1.0
        !if(tau<0) then
          !tau=beta+tau
          !s=-s
        !endif
        !if(tau>=beta) then
          !tau=tau-beta
          !s=-s
        !endif
        !Ek=sum(Mom**2)   !kinetic energy
        !x=Beta*(Ek-Mu)/2.0
        !y=2.0*tau/Beta-1.0
        !if(x>100.0) then
          !Green=dexp(-x*(y+1.0))
        !else if(x<-100.0) then
          !Green=dexp(x*(1.0-y))
        !else
          !Green=dexp(-x*y)/(2.0*cosh(x))
        !endif
        !!if(spin==1 .or. spin==-1) then
        !Green=s*Green
    
        !if(isnan(Green)) then
          !print *, "Green is too large!", tau, Ek, Green
          !stop
        !endif
        !return
      !end function Green

    double precision function Green(tau ,Mom, spin)
        !calculate Green's function
        implicit none
        double precision :: tau, k2, s, Ek, x, y, w, r, coshv
        integer :: spin, i
        double precision, dimension(D) :: Mom
    !    print *, tau, Mom, Spin
    !    stop
        Ek=sum(Mom**2)   !kinetic energy
        Green=1.0/(Ek+Mass2)
    
        return
      end function Green

      !subroutine CreateGTable()
        !implicit none
        !integer :: iEk, iTau, i
        !double precision :: dEk, dTau, Ek, Tau, x
        !dEk=5.0*EF/MaxEK
        !dTau=Beta/MaxTau
        !do iEk=0, MaxEK
          !do iTau=0, MaxTau-1 !treat the last point, namely tau==beta later
            !Ek=iEk*dEk
            !Tau=iTau*dTau
            !GreenTable(iEk, iTau)=GreenFunc(Tau, Ek, 1)
          !enddo
          !GreenTable(iEk, MaxTau)=GreenFunc(Beta-1.0e-10, Ek, 1) !set the largest table point to tau--->beta^-
        !enddo
      !end subroutine

      !double precision function Green(tau, Mom, spin)
        !implicit none
        !integer :: iEK, iTau, spin
        !double precision :: Ek, s, tau, wEk, wTau, rEk, rTau
        !double precision :: g11, g12, g21, g22
        !double precision :: ratioEk, ratioTau
        !double precision, dimension(D) :: Mom
        !Ek=sum(Mom**2)   !kinetic energy
        !s=1.0
        !if(tau<0) then
          !tau=beta+tau
          !s=-s
        !endif
        !if(tau>=beta) then
          !tau=tau-beta
          !s=-s
        !endif
        !if(Ek>=5.0*EF) then
          !Green=GreenFunc(tau, Ek, spin)
        !else
          !wEk=Ek/5.0/EF*MaxEK
          !wTau=tau/Beta*MaxTau
          !iEK=int(wEk)
          !iTau=int(wTau)
          !rEk=wEk-iEk
          !rTau=wTau-iTau
          !Green=(1-rTau)*((1-rEk)*GreenTable(iEk, iTau)+rEk*GreenTable(iEk+1, iTau))
          !Green=Green+rTau*((1-rEk)*GreenTable(iEk, iTau+1)+rEk*GreenTable(iEk+1, iTau+1))
          !Green=s*Green
          !if(abs(Green-s*GreenFunc(tau, Ek, spin))>1e-4) then
            !print *, Ek, tau
            !print *, iEk, iTau
            !print *, wEk, wtau
            !print *, rEk, rTau
            !print *, GreenTable(iEk, iTau), GreenTable(iEk+1, iTau), GreenTable(iEk, iTau+1), GreenTable(iEk+1, iTau+1)
            !print *, Green, s*GreenFunc(tau, Ek, spin) 
            !stop
          !endif
        !endif
      !end function Green
    
      double precision function Interaction(tau, Mom, spin)
        implicit none
        double precision :: tau,k2
        double precision, dimension(D) :: Mom
        integer :: spin
        Interaction=1.0/(sum(Mom**2)+Mass2)
    !    Interaction=1.0/(Mass2)
        return
      end function Interaction
    
      subroutine DynamicTest()
        implicit none
        integer   :: i,j
        double precision   :: weight
        double precision :: tau,k2, Ek1, Ek2
        !if(Step==19.0) then
          !print *, "Before Reset"
          !print *, "LoopFreq", LoopFreq
          !print *, "LoopMom", LoopMom
          !print *, "GWeight", GWeight(1:GNum), GWeightNorm(1:GNum)
          !print *, "VerWeight", VerWeight(1:VerNum), VerWeightNorm(1:VerNum)
        !endif
    
        call ResetWeightTable()
        weight = CalcWeight(0)
        if(abs(weight-OldWeight)>1e-6) then
          print *, "Weight is wrong at step:", Step, weight, OldWeight
          print *, "TauTable", TauTable(1:TauNum)
          print *, "LoopMom", LoopMom(:, 1:LoopNum)
          print *, "GWeight", GWeight(1:GNum), GWeightNorm(1:GNum)
          print *, "VerWeight", VerWeight(1:VerNum), VerWeightNorm(1:VerNum)
          stop
        endif
    
        do i=1, LoopNum
          do j=1, D
            if(isnan(LoopMom(j, i))) then
              print *, "Mom is NaN",Step, j, i, LoopMom(j, i)
              stop
            endif
          enddo
        enddo
    
        do i=1, iVerNum
          if(isnan(VerWeight(i))) then
            print *, "VerWeight is NaN",Step, i, VerWeight(i)
            stop
          endif
    
          if(isnan(VerWeightNorm(i))) then
            print *, "VerWeightNorm is NaN",Step, i, VerWeightNorm(i)
            stop
          endif
        enddo
    
        if(isnan(OldWeight)) then
          print *, "OldWeight is NaN",Step, i, OldWeight
          stop
        endif
    
        if(OldWeight==0.0) then
          print *, "OldWeight is zero",Step, i, OldWeight
          stop
        endif
    
        if(abs(ExtMomMesh(1, ExtMomBin)-LoopMom(1,1))>1e-6) then
          print *, "ExtMom is wrong!",Step, ExtMomBin, ExtMomMesh(:, ExtMomBin), LoopMom(:,1)
        endif
    end subroutine
   
    subroutine ResetWeightTable()
        implicit none
        integer :: sector, i, j, Spin
        double precision :: Tau
        double precision, dimension(D) :: Mom
        !print *, "iGNum", iGNum
        call NewState()
        do i=1, iGNum 
          Tau = TauTable(TauBases(2, i))-TauTable(TauBases(1, i))
          Spin=sum(LoopBases(:, i)*LoopSpin(:))
          do j=1, D
            Mom(j)=sum(LoopBases(:, i)*LoopMom(j, :))
          enddo
          NewGWeight(i)=Green(Tau, Mom, Spin)
        enddo
    
    
        do i=1, iVerNum 
            Spin=sum(LoopBasesVer(:, i)*LoopSpin(:))
            do j=1, D
                Mom(j)=sum(LoopBasesVer(:, i)*LoopMom(j, :))
            enddo

            if (reducibleTable(i)) then
                NewVerWeight(i) = 0.0
            else
                NewVerWeight(i) = Interaction(0.d0, Mom, Spin)
            endif

        enddo
    
        return
    end subroutine

    subroutine NewState()
      implicit none
      NewGWeight(1:iGNum)=GWeight(1:iGNum)
      NewVerWeight(1:iVerNum)=VerWeight(1:iVerNum)
    end subroutine
    
    subroutine UpdateState()
      implicit none
      GWeight(1:iGNum)=NewGWeight(1:iGNum)
      VerWeight(1:iVerNum)=NewVerWeight(1:iVerNum)
    end subroutine
    
    
    subroutine ChangeWeightTableMom(loopindex)
        implicit none
        integer :: sector, i, j, Spin, loopindex
        double precision :: Tau
        double precision, dimension(D) :: Mom
        call NewState()
        do i=1, iGNum 
            if (LoopBases(loopindex, i)==0) then
                cycle !if the loop with #loopindex do not affect the Gline, then do not recalculate 
            endif
    
            Tau = TauTable(TauBases(2, i))-TauTable(TauBases(1, i))
            !Spin = sum(LoopBases(:, i)*LoopSpin(:))
            Spin=1
            do j=1, D
                Mom(j) = sum(LoopBases(:, i)*LoopMom(j, :))
            enddo
            NewGWeight(i) = Green(Tau, Mom, Spin)
        enddo
    
    
        do i=1, iVerNum 

            if (LoopBasesVer(loopindex, i)==0) then
              cycle !if the loop with #loopindex do not affect the Gline, then do not recalculate 
            endif
      
            if (reducibleTable(i)) then
                NewVerWeight(i) = 0.0
            else
              !Spin=sum(LoopBasesVer(:, i)*LoopSpin(:))
              Spin=1
              do j=1, D
                Mom(j)=sum(LoopBasesVer(:, i)*LoopMom(j, :))
              enddo
              NewVerWeight(i) = Interaction(0.d0, Mom, Spin)
            endif
        enddo
    
        return
    end subroutine
    
    subroutine ChangeWeightTableTau(Tauindex)
        implicit none
        integer :: sector, i, j, Spin, Tauindex
        double precision :: Tau
        double precision, dimension(D) :: Mom
        call NewState()
        do i=1, iGNum 
            if (TauBases(1, i)==Tauindex .or. TauBases(2, i)==Tauindex .or. &
                TauBases(1, i)==Tauindex+1 .or. TauBases(2, i)==Tauindex+1) then
                Tau=TauTable(TauBases(2, i))-TauTable(TauBases(1, i))
                !Spin=sum(LoopBases(:, i)*LoopSpin(:))
                Spin=1
                do j=1, D
                    Mom(j)=sum(LoopBases(:, i)*LoopMom(j, :))
                enddo
                NewGWeight(i)=Green(Tau, Mom, Spin)
            endif
        enddo
        return
    end subroutine
    
    double precision function CalcWeight(SaveWeight)
        !calculate the weight for ALL diagrams in a given sector
        implicit none
        integer :: sector, i, j, k, BlockNum, Shift
        double precision :: Q2, Q
        double precision :: GWeight, TempWeight, VerWeight, AbsVerWeight
        integer :: SaveWeight
        !double precision, dimension(D) :: Mom1, Mom2
        !double precision :: Ek1, Ek2
        CalcWeight=0.0
    
        do i=1, HugenholtzNum
            GWeight = SymFactor(i)  !initialize weight wiht the symmetry factor
            do j=1, GNum
                GWeight = GWeight * NewGWeight(GIndex(i, j))
            enddo

            !TempWeight=1.0
            !do j=1, VerNum
                !TempWeight = TempWeight * ( NewVerWeight(VerIndex(i,2*j-1)) - NewVerWeight(VerIndex(i,2*j)) )
            !enddo

            !!!!!!!!!!!!!! Calculate Feynman diagram weight with binary tree expansion !!!!!!!!!!!!!!!

            SpinCache(1, 1)=NewVerWeight(VerIndex(i,1))
            SpinCache(2, 1)=NewVerWeight(VerIndex(i,2))
            BlockNum=2
            do j=2, VerNum
              !print *, "j:", j
              do k=1, BlockNum
                SpinCache(2*k-1, j)=SpinCache(k, j-1)*NewVerWeight(VerIndex(i,2*j-1))
                SpinCache(2*k, j)=SpinCache(k, j-1)*NewVerWeight(VerIndex(i,2*j))
              enddo
              !print *, "Spincahce", j, SpinCache(:, j)
              BlockNum=BlockNum*2
            enddo

            VerWeight=0.0
            AbsVerWeight=0.0
            do j=1, 2**(Order-1)
              VerWeight=VerWeight+SpinCache(j, VerNum)*SpinFactor(j, i)
              AbsVerWeight=AbsVerWeight+abs(SpinCache(j, VerNum)*SpinFactor(j, i))
            enddo

            !if(abs(TempWeight-VerWeight)>1e-7) then
              !print *, "binary tree fail", TempWeight, VerWeight
              !!print *, SpinFactor(:,i)
              !print *, SpinCache(:, 1)
              !print *, SpinCache(:, 2)
              !print *, NewVerWeight(VerIndex(i, 1))
              !print *, NewVerWeight(VerIndex(i, 2))
              !print *, NewVerWeight(VerIndex(i, 3))
              !print *, NewVerWeight(VerIndex(i, 4))
              !stop
            !endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(SaveWeight==1) then
              DiagWeight(i) = GWeight*VerWeight
              DiagWeightABS(i) = abs(GWeight)*abs(AbsVerWeight)
            endif

            !CalcWeight = CalcWeight + abs(GWeight*VerWeight)
            CalcWeight = CalcWeight + GWeight*VerWeight
        enddo
    
        return
    end function CalcWeight
    
    subroutine Measure()
        implicit none
        integer :: Spin, Num, i, j, TauBin
        double precision :: AbsWeight, dQ, Q, Freq, ReWeight, Weight
        double precision, dimension(D) :: Mom
        double precision :: Phase
        double precision, dimension(4) :: t, u, s, up, sp
    
        AbsWeight=abs(OldWeight)
        Phase=OldWeight/AbsWeight
    
        Num=ExtMomBin
        Q=norm2(ExtMomMesh(:, ExtMomBin))
        !ReWeight=exp(Q*Q/1.d0/kF/kF)
        Reweight=1.0
    
        !if(Num<=QBinNum) then
          !if(D==2) then
            !!Polarization(Num)=Polarization(Num)+Phase/2.0/pi/Q*Reweight
            !Polarization(Num)=Polarization(Num)+Phase*Reweight
          !else
            !!Polarization(Num)=Polarization(Num)+Phase/4.0/pi/Q/Q*Reweight
            !Polarization(Num)=Polarization(Num)+Phase*Reweight
          !endif
        !endif
    
        call NewState()
        Weight=CalcWeight(1)
    
        !------------- measure sign blessing ------------------
   
        do i=1, HugenholtzNum
            F1(i) = F1(i) + DiagWeight(i)/AbsWeight 
            F2(i) = F2(i) + ABS(DiagWeight(i))/AbsWeight
            F3(i) = F3(i) + DiagWeightABS(i)/AbsWeight
        enddo

        F1(HugenholtzNum+1) = F1(HugenholtzNum+1) + SUM(DiagWeight(1:HugenholtzNum))/AbsWeight 
        F2(HugenholtzNum+1) = F2(HugenholtzNum+1) + ABS(SUM(DiagWeight(1:HugenholtzNum)))/AbsWeight
        F3(HugenholtzNum+1) = F3(HugenholtzNum+1) + SUM(ABS(DiagWeightABS(1:HugenholtzNum)))/AbsWeight

        !-------------------------------------------------------
        do i = 1, HugenholtzNum
          PolarDiag(i, Num) = PolarDiag(i, Num) + DiagWeight(i)/AbsWeight
          Polarization(Num) = Polarization(Num) + DiagWeight(i)/AbsWeight
        enddo
        return
    end subroutine
    
    subroutine SaveToDisk()
        implicit none
        integer :: i, ref, j
        double precision :: Obs
        !Save NormWeight and Polarization to disk
        character*10 :: ID
        character*10 :: DiagIndex
        character*20 :: filename
        write(ID, '(i10)') PID
        filename="Data/Diag"//trim(adjustl(ID))//".dat"
        write(*,*) "Save to disk ..."
        open(100, status="replace", file=trim(filename))
        write(100, *) "#", Step
        write(100, *) "#", Polarization(1)
        ref=int(kF/DeltaQ)+1
        do i=1, QBinNum
            Obs = Polarization(i)
            write(100, *) norm2(ExtMomMesh(:, i)), Obs
        enddo
        close(100)
    
        do j=1, HugenholtzNum
            write(DiagIndex, '(i10)') j
            filename="Data/Diag"//trim(adjustl(ID))//"_"//trim(adjustl(DiagIndex))//".dat"
            open(100, status="replace", file=trim(filename))
            write(100, *) "#", Step
            write(100, *) "#", PolarDiag(j,1)
            ref=int(kF/DeltaQ)+1
            do i=1, QBinNum
              !Obs=Polarization(i)/sum(abs(Polarization))
              !Obs=Polarization(i)/(sum(abs(Polarization))/QBinNum)
              Obs=PolarDiag(j, i)
              !Obs=Polarization(i)/real(Polarization(ref))
              write(100, *) norm2(ExtMomMesh(:, i)), Obs
            enddo
            close(100)
        enddo
        return
    end subroutine
    
    subroutine SaveToDiskF()
        implicit none
        character*20 :: filename
        integer :: i
        double precision :: res=0.0
        do i=1, HugenholtzNum
            res = res + F2(i)
        enddo

        filename="F123.dat"
        open(100, position='append', file=trim(filename))
            write(100, *)  "OLD BASES:  q=", LoopMom(1,1) 
            do i=1, HugenholtzNum+1
               write(100, *) "Hugenholtz", i,":", F1(i)/F3(HugenholtzNum+1), F2(i)/F3(HugenholtzNum+1), F3(i)/F3(HugenholtzNum+1)
            enddo
            write(100, *) "F2   SUM:", res/F3(HugenholtzNum+1)
        close(100)
        
    end subroutine
    
    
    double precision function  CalcNorm()
        implicit none
    
        CalcNorm = normConstant
    
    end function
    
    
      
    subroutine NormToPhysical()
      !increase diagram order by one/change normalization diagram to physical diagram
        implicit none
        if (Sector==PHYSICAL) return
        !Weight = CalcNorm()
        return
    end subroutine
    
      subroutine PhysicalToNorm()
      !decrease diagram order by one/change physical diagram to normalization diagram
        implicit none
        if (Sector==NORMALIZATION) return
        !if the current diagrams are already in normalization sector, then return
        return
      end subroutine
    
      !subroutine RollBack()
    
      subroutine ChangeTau()
      !randomly choose a vertex, change the time variable
        implicit none
        double precision :: NewTau, prop, dw, R
        double precision :: Weight
        integer :: Num
    
        !print *, "Change Freq", Step
    
        Num=int(grnd()*Order)*2+1 !1,3,5
    
        PropStep(1)=PropStep(1)+1.0
        OldTau=TauTable(Num+1) !in the case of Num==1, then TauTable(1)/=TauTable(2)
        call GenerateNewTau(OldTau, NewTau, prop)
        TauTable(Num+1)=NewTau
        if(Num/=1) TauTable(Num)=NewTau
        call ChangeWeightTableTau(Num)
        Weight = CalcWeight(0)
        R=prop*abs(Weight)/abs(OldWeight)
        if(grnd()<R) then
          !print *, abs(Weight), abs(OldWeight)
          AcceptStep(1)=AcceptStep(1)+1.0
          OldWeight=Weight
          call UpdateState()
        else
          if(Num/=1) TauTable(Num)=OldTau
          TauTable(Num+1)=OldTau
          !call ChangeWeightTableTau(Num)
        endif
        !GenerateNewFreq
        return
      end subroutine
    
    subroutine ChangeMom()
        !randomly choose a vertex, change the space variable
        implicit none
        double precision, dimension(D) :: NewMom
        double precision :: prop, R
        integer :: Num, i, j, NewExtMomBin
        double precision :: Weight
    
        !Num = int( (LoopNum-1)*grnd() ) + 2
        Num = int( (LoopNum)*grnd() ) + 1
    
        PropStep(2) = PropStep(2) + 1.0
        OldMom = LoopMom(:, Num)
    
        if(Num==1) then
            call GenerateNewExtMom(ExtMomBin, NewExtMomBin, prop)
            NewMom=ExtMomMesh(:, NewExtMomBin)
        else
    
          call GenerateNewMom(OldMom, NewMom, prop)
    
        endif
    
        if(prop<0.0) return
        if(Num==1 .and. norm2(NewMom)>ExtMomMax) return
        LoopMom(:,Num)=NewMom
    
        call ChangeWeightTableMom(Num)
        Weight = CalcWeight(0)
        R = prop*abs(Weight)/abs(OldWeight)
    
        if(grnd()<R) then
            AcceptStep(2) = AcceptStep(2)+1.0
            OldWeight = Weight
            if(Num==1) ExtMomBin=NewExtMomBin
            call UpdateState()
        else
          LoopMom(:,Num)=OldMom
          !call ChangeWeightTableMom(Num)
        endif
    
        return
    end subroutine
    
      subroutine SwapMom()
      !randomly choose a vertex, change the space variable
        implicit none
        double precision, dimension(D) :: TempMom1, TempMom2
        double precision :: prop, R
        integer :: Num1, Num2, i, j, NewExtMomBin
        double precision :: Weight
        !print *, "Change Mom", Step
    
        Num1=int(grnd()*LoopNum)+1
        Num2=int(grnd()*LoopNum)+1
        if(Num1==1 .or. Num2==1) return
        if(Num1==Num2) return
    
    !    PropStep(3)=PropStep(3)+1.0
    
        TempMom1=LoopMom(:, Num1)
        TempMom2=LoopMom(:, Num2)
    
        LoopMom(:,Num1)=TempMom1
        LoopMom(:,Num2)=TempMom2
    
        call ResetWeightTable()
        Weight=CalcWeight(0)
        R=abs(Weight)/abs(OldWeight)
    
        if(grnd()<R) then
    !      AcceptStep(3)=AcceptStep(3)+1.0
          OldWeight=Weight
          call UpdateState()
        else
          LoopMom(:,Num1)=TempMom1
          LoopMom(:,Num2)=TempMom2
          !call ResetWeightTable()
        endif
    
        return
    end subroutine
    
    subroutine ChangeSpin()
        !randomly choose a loop basis, change the spin variable
        implicit none
        return
    end subroutine
    
    subroutine GenerateNewMom(old, new, prop)
        implicit none
        double precision, dimension(D) :: old, new
        integer :: Direction, Num
        double precision, parameter :: STEP=5.0
        double precision :: ratio, prop, k, k_new
        double precision :: lambda, x
        x=grnd()
        if(x<1.0/3.0) then
            new=old
            Num=int(grnd()*D)+1
            new(Num)=new(Num)+sign(STEP/Beta*grnd(), grnd()-0.5)
            prop=1.0
        else if(x<3.0/3.0) then
            k=norm2(old) !sqrt(m_x^2+m_y^2+m_z^2)
            if(k==0.0)then
            prop=-1.0
            return
            endif
            lambda=1.5
            k_new=k/lambda+grnd()*(lambda-1.0/lambda)*k
            ratio=k_new/k
            new=old*ratio
            if(D==2) then
            prop=1.0 ! prop=k_old/k_new
            else
            prop=k_new/k ! prop=k_old/k_new
            endif
        else
            new=-old
            prop=1.0
        endif
    end subroutine
    
    
    subroutine GenerateNewExtMom(old, new, prop)
        implicit none
        double precision :: prop, x
        integer :: old, new
    
        new=int(grnd()*QBinNum)+1
        prop=1.0
    end subroutine
    
    
    subroutine GenerateNewTau(OldTau, NewTau, Prop)
        implicit none
        double precision, intent(in) :: OldTau
        double precision, intent(out) :: NewTau
        double precision, intent(out) :: Prop
        double precision :: DeltaT
        if(grnd()<1.0/3.0) then
          DeltaT=Beta/3.0
          NewTau=OldTau+DeltaT*(grnd()-0.5)
        else if(grnd()<2.0/3.0) then
          NewTau=-OldTau
        else
          NewTau=grnd()*Beta
        endif
        if (NewTau<0) then
          NewTau=NewTau+Beta
        else if (NewTau>Beta) then
          NewTau=NewTau-Beta
        endif
        Prop=1.0
        return
    
    end subroutine GenerateNewTau

    subroutine Reducible
        implicit none
        integer :: k1, k2, flag 
        integer, dimension(MaxOrder+1) :: OutLine=0

        OutLine(1) = 1
        do k1=1, iVerNum 

            flag = 1
            do k2=1, LoopNum
                if ( ABS(LoopBasesVer(k2, k1)) .ne. OutLine(k2) ) then
                    flag = 0
                    exit
                end if
            end do 

            if (flag==1) then 
                reducibleTable(k1) = .true.
            else if (flag==0) then
                reducibleTable(k1) = .false.
            end if 
        end do

    end subroutine


end program main
