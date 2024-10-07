!  Routh.f90 
!
!****************************************************************************
!
!  PROGRAM: Routh
!
!****************************************************************************

    program Routh

    implicit none
        
    ! Variables
    real::x,y,z
    integer::deg,syc,rows,cols,k,changecount
    logical::sc1=0,sc2=0
    real, dimension (:), allocatable :: coef 
    real, dimension (:,:), allocatable :: rt 
    logical, dimension (:,:), allocatable :: signs 
    integer,dimension(:),allocatable::frow    
    integer,dimension(:),allocatable::srow
    
    ! Body of Routh
    print *, 'Enter the degree of characteristics equation of the system'
    read*,deg
    rows=deg+1
    cols=deg/2+1
    allocate (coef(rows))
    print *, 'Enter the coefficients of characteristics equation'
    read*,coef
    allocate (rt(rows,cols))   
    allocate (signs(rows,cols)) 
    rt=0;
    if (MOD(deg,2)==0) then
        allocate(frow(cols))
        allocate(srow(cols-1))
        frow=0
        srow=0
        do syc=1,cols
            frow(syc)=syc*2-1
        end do
        do syc=1,deg/2
            srow(syc)=syc*2
        end do
        rt(1,:)=coef(frow)
        rt(2,srow/2)=coef(srow)
    else
        allocate(frow(cols))
        allocate(srow(cols))
        frow=0
        srow=0
        do syc=1,cols
            frow(syc)=syc*2-1
            srow(syc)=syc*2
        end do
        rt(1,:)=coef(frow)
        rt(2,:)=coef(srow)
    end if
    do syc=3,rows
        do k=1,cols-1
            rt(syc,k)=-1/rt(syc-1,1)*(rt(syc-2,1)*rt(syc-1,k+1)-rt(syc-1,1)*rt(syc-2,k+1))
        end do
        if ((rt(syc,1)==0) .AND. (SUM(rt(syc,:))/=0)) then
            sc1=1
            print*,"Special Case-1. Zero at first column of",syc,". row."
            print*,"This value has changed to 0.0001 in order to complete the table."
            rt(syc,1)=0.0001
        end if
        if (SUM(rt(syc,:))==0) then
            print'("Special Case-2. Entire row zero. Row: s^",I1)',rows-syc 
            print*,"This case shows the system is unstable or marginally stable."
            print*,"The zero row elements have been changed with the auxiliary polynom coefficients obtained from the previous row."
            do k=1,(deg-syc+1)/2+1
                rt(syc,k)=rt(syc-1,k)*(deg+1-(syc-1)-(k-1)*2)
            end do
            sc2=.true.
        end if
    end do   
    signs=rt>0
    print*,"Routh Table:"
    do syc=1,rows
        print'("s^",I1,"     ",(*(F13.5)))',rows-syc,rt(syc,:)        
    end do    
    do syc=2,rows
        if (signs(syc,1)/=signs(syc-1,1)) then
            changecount=changecount+1
        end if
    end do     
    if((sc2==.true.) .OR. (sc1==.true.)) then
        if ((changecount==0) .AND. (sc2==.true.)) then
            print*,"No sign changes at first column, the system is marginally stable."
        else
            print'(I2," sign changes at first column, the system is unstable.")',changecount
        end if
    else
        print*,"No special cases."
        if (changecount==0) then
            print*,"No sign changes at first column, the system is stable."
        else
            print'(I2," sign changes at first column, the system is unstable.")',changecount
        end if
    end if
    end program Routh

