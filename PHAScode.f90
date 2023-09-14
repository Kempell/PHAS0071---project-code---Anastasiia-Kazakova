program PHAScode 

  implicit none

  integer :: nx,ny,istep,nstep,nn,ibc,intv
  double precision, dimension(:), allocatable :: h,hp,x,y,kf,kd,b,u,chi,edot,exhum_low,exhum_deep,depth
  real :: time_in,time_out
  double precision :: kfsed,m,n,kdsed,g1,g2,expp
  double precision xl,yl,dt,pi,vex,ratio,dista,f,v,g,fa,ttstep,timestp

  integer i,j,itime

  character cs*5,ct*5
  
  double precision :: dd,ee,ccx,ccy,dist,eea,tt1,tt2,dda,ccy2,ccx2


!>>>>>>>>>>>>>>>>>>>>>>>> Input variables here <<<<<<<<<<<<<<<<<<<<<<<<<<


!>>>>>>>>>> For 1st crater
!crater diamter (m) for crater A
ee=25e3;

!crater depth (km*10)
dd=8.e-2;

!crater x position
ccx=35.d3;

!crater y positiom
ccy=100.d3;

!ratio of crater depth to crater rim height with 1 = no rim; 2 = rim as high as
!the crater itself
f=1.7;

!crater spawn timestep
tt1=1500;




!>>>>>>>>>> For 2nd crater
!crater diamter (m) for crater B
eea=16e3;

!crater depth (km*10)
dda=8.e-2;

!crater x position
ccx2=35.d3;

!crater y positiom
ccy2=50.d3;

!ratio of crater depth to crater rim height with 1 = no rim; 2 = rim as high as
!the crater itself
fa=1.7;

!crater spawn timestep
tt2=2000;



!>>>>>>>>> other parameters
!number of timesteps
ttstep=1.e5;

!set time step (years*100)
timestp=1e4;

!interval at which output files with x y z data are saved (istep) 
intv=50;


!>>>>>>>>>>>>>>>>>>>>>>>>>End variable input <<<<<<<<<<<<<<<<<<<<<<<<<<<<




! set model resolution (200 x 200 m)
  nx = 371
  ny = 1001
  nn = nx*ny

  pi=atan(1.d0)*4.d0

! initialize FastScape
  call FastScape_Init ()
  call FastScape_Set_NX_NY (nx,ny)
  call FastScape_Setup ()

! wset model dimensions
  xl=75.d3
  yl=200.d3
  call FastScape_Set_XL_YL (xl,yl)

! construct nodal coordinate arrays x and y
  allocate (x(nx*ny),y(nx*ny),u(nx*ny))
  allocate (chi(nx*ny),edot(nx*ny),exhum_low(nx*ny),exhum_deep(nx*ny))
  x = (/((xl*float(i-1)/(nx-1), i=1,nx),j=1,ny)/)
  y = (/((yl*float(j-1)/(ny-1), i=1,nx),j=1,ny)/)

! set time step (coresponds to 100 ky with the current erosion and sediemnt transport rates)
  dt=timestp
  call FastScape_Set_DT (dt)

  


!>>>>> erosion and other stream power model parameters
  allocate (kf(nn),kd(nn))
  kf=3.d-7 !bedrock eordibility 
  kfsed=3.d-7 !sediment erodibility
  m=0.5d0
  n=1.d0 ! need to be >=1
  kd=1.5d-3 !bedrock transport coefficient
  kdsed=1.5d-2 !sediement transport coefficient
  g1=1.d0
  g2=1.d0
  expp=1.d0
  call FastScape_Set_Erosional_Parameters (kf,kfsed,m,n,kd,kdsed,g1,g2,expp)





! set uplift rate (uniform while keeping boundaries at base level)
 u = 2.d-5
 do i=1,nx*ny 
   if (y(i)<10) then
   u(i)=0
   endif
 enddo
 
call FastScape_Set_U (u)


! bottom side is fixed only
  ibc=1000
call FastScape_Set_BC (ibc)

!initial topography is a 2000 m high plateau
  allocate (h(nn),b(nn),hp(nn))
  call random_number (h)
  where (y.gt.yl/3d0) h=h+2000

  
  do i=1,nx*ny
    h(i)=h(i)+(100.*(y(i)/yl))+(1d-4*abs(x(i)-(xl/2)));
  enddo

call FastScape_Init_H (h)

! set number of time steps
  nstep = ttstep 

! echo model setup
  call FastScape_View ()

! initializes time step
  call FastScape_Get_Step (istep)

! set vertical exaggeration
  vex = 1.d0

! start of time loop
  call cpu_time (time_in)


!>>>> run the model <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
do while (istep.lt.nstep)
        u = 2.d-5
        do i=1,nx*ny 
           if (y(i)<10) then
           u(i)=0
           endif
        enddo


     !>>>>>> first crater
        if (istep==tt1) then
        do i=1,nx*ny
        dist=sqrt((ccx-x(i))**2+(ccy-y(i))**2);
                if ((dist)<(ee/f)) then
                u(i)=-u(i)-((-30/(1-dist)/ee))-dd;
                endif
        
                if (dist>(ee/f) .and. dist<ee) then
                u(i)=u(i)+(((ee/dist)-1)*dd*f);
                endif

                if (y(i)<10) then
                u(i)=0
                endif
        enddo
        endif

    !>>>> second crater
        if (istep==tt2) then
        do i=1,nx*ny
        dist=sqrt((ccx2-x(i))**2+(ccy2-y(i))**2);
                if ((dist)<(eea/fa)) then
                u(i)=-u(i)-((-30/(1-dist)/eea))-dda;        
                endif       

                if (dist>(eea/fa) .and. dist<eea) then
                u(i)=u(i)+(((eea/dist)-1)*dda*fa);
                endif

                if (y(i)<10) then
                u(i)=0
                endif
        enddo
        endif

    call FastScape_Set_U (u)
    ! execute FastScape step
    call FastScape_Execute_Step ()
    call FastScape_Get_Step (istep)
  
        if (mod(istep,intv)==0) then
        call FastScape_Copy_H (h)
        call FastScape_Copy_Basement (b)
        call FastScape_VTK (h-b, vex)
        call FastScape_Copy_Erosion_Rate (edot)

        itime = istep
        write(cs,'(i5)') itime
        if (itime.lt.10) cs(1:4)='0000'
        if (itime.lt.100) cs(1:3)='000'
        if (itime.lt.1000) cs(1:2)='00'
        if (itime.lt.10000) cs(1:1)='0'
        open(345,file="step__"//cs//".txt",status="unknown")
                do i=1,nx*ny  
                write(345,*) x(i),y(i),h(i),edot(i)
                enddo
        close(345)
        endif    

enddo

! display timing information
  call FastScape_Debug()
  call cpu_time (time_out)
  print*,'Total run time',time_out-time_in

! exits FastScape
  call FastScape_Destroy ()

! deallocate memory
  deallocate (h,x,y,kf,kd,b,hp)

end program PHAScode
