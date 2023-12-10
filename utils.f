c--------------------------------------------------------------------------------------------------
c Performs 2-point (i.e. first order) extrapolation stencil using values 
c Tp1 at p1 with coordinate xp1 and Tp2 at p2 with coordinate xp2 
c (p1 is the closest point to the one we want to extrapolate, p2 is the furthest one)
c--------------------------------------------------------------------------------------------------
        real*8 function firstord_extrap(Tp1,Tp2,xp1,xp2,x_ex)
        implicit none
        real*8 Tp1,Tp2,xp1,xp2,x_ex

        !--------------------------------------------------------------

        firstord_extrap=Tp1*(x_ex-xp2)/(xp1-xp2)
     &                 +Tp2*(x_ex-xp1)/(xp2-xp1)

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c Performs bilinear interpolation on a plane using values 
c Tp2a at (xp2a,yp2a), Tp2b at (xp2b,yp2b), Tp2c at (xp2c,yp2c), Tp2d at (xp2d,yp2d)  
c to obtain the value at (x_interp,y_interp) 
c----------------------------------------------------------------------
        real*8 function bilinear_interp(Tp2a,Tp2b,Tp2c,Tp2d,
     &                                  xp2a,xp2b,xp2c,xp2d,
     &                                  yp2a,yp2b,yp2c,yp2d,
     &                                  x_interp,y_interp)
        implicit none
        real*8 Tp2a,Tp2b,Tp2c,Tp2d
        real*8 xp2a,xp2b,xp2c,xp2d
        real*8 yp2a,yp2b,yp2c,yp2d
        real*8 x_interp,y_interp

        !--------------------------------------------------------------

        bilinear_interp=(1/((xp2c-xp2a)*(yp2c-yp2a)))*
     &    (Tp2a*(xp2c-x_interp)*(yp2c-y_interp)+
     &     Tp2b*(xp2c-x_interp)*(y_interp-yp2a)+
     &     Tp2c*(x_interp-xp2a)*(y_interp-yp2a)+
     &     Tp2d*(x_interp-xp2a)*(yp2c-y_interp))

        return
        end
c--------------------------------------------------------------------------------------


c----------------------------------------------------------------------
c Converts Cartesian coordinates x0,y0,z0 to spherical coordinates rho0,chi0,xi0
c----------------------------------------------------------------------
        subroutine cart_to_sph(rho0,chi0,xi0,x0,y0,z0)
        implicit none
        real*8 x0,y0,z0
        real*8 rho0,chi0,xi0

        real*8 PI
        parameter (PI=3.141592653589793d0)


        !--------------------------------------------------------------

         rho0=sqrt(x0**2+y0**2+z0**2)
         chi0=(1/PI)*acos(x0/rho0)
         if (z0.lt.0) then
             xi0=(1/(2*PI))*(atan2(z0,y0)+2*PI)
         else
             xi0=(1/(2*PI))*atan2(z0,y0)
         end if

        return
        end
c--------------------------------------------------------------------------------------


c----------------------------------------------------------------------
c Converts spherical coordinates rho0,chi0,xi0 to Cartesian coordinates x0,y0,z0
c----------------------------------------------------------------------
        subroutine sph_to_cart(x0,y0,z0,rho0,chi0,xi0)
        implicit none
        real*8 x0,y0,z0
        real*8 rho0,chi0,xi0

        real*8 PI
        parameter (PI=3.141592653589793d0)


        !--------------------------------------------------------------

         x0=rho0*cos(PI*chi0)
         y0=rho0*sin(PI*chi0)*cos(2*PI*xi0)
         z0=rho0*sin(PI*chi0)*sin(2*PI*xi0)

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c Prints parameters to screen
c----------------------------------------------------------------------
        subroutine print_params_toscreen(rho_ex,Nx,Ny,Nz,
     &                          f_id,nametag,n_procs)

        implicit none
        real*8 rho_ex
        integer Nx,Ny,Nz
        integer f_id
        integer n_procs

        character*40 nametag

        !--------------------------------------------------------------

        write (*,*) "Working on ",n_procs," processes"
        write (*,*) "Perform radial extrapolation at rho_ex=",rho_ex
        write (*,*) "with resolution Nx,Ny,Nz=",Nx,Ny,Nz
        if (f_id.eq.0) then
            write(*,*) "Extrapolating f=rho:=sqrt(x**2+y**2+z**2)"
          else
            write (*,*) "ERROR: f_id>0 is not implemented."
            stop
        end if
        write (*,*) "Output file tag is ",nametag

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c Prints separation line
c----------------------------------------------------------------------
        subroutine separation_line()

        !--------------------------------------------------------------

        
        write (*,*) "================================================"

        end
c--------------------------------------------------------------------------------------

c--------------------------------------------------------------------------------------
c For a given process, determines the number of points along each Cartesian direction 
c and the initial value of the corresponding Cartesian coordinate
c--------------------------------------------------------------------------------------
        subroutine setup_local_grid(
     &       Nx,Ny,Nz,
     &       dx,dy,dz,
     &       xinit,yinit,zinit,
     &       Nx_all,Ny_all,Nz_all,n_procs,procid)

        implicit none

        integer n_procs,procid
        integer Nx_all,Ny_all,Nz_all
        real*8 dx,dy,dz
        integer Nx,Ny,Nz
        real*8 xinit,yinit,zinit

!----------------------------------------------------------------------

        dx=2.0d0/(Nx_all-1)
        dy=2.0d0/(Ny_all-1)
        dz=2.0d0/(Nz_all-1)

        Nx=floor(real(Nx_all)/n_procs)
        if (procid.eq.(n_procs-1)) Nx=Nx_all-(n_procs-1)*Nx
        xinit=-1+procid*Nx*dx
        if (procid.eq.(n_procs-1)) xinit=1-(Nx-1)*dx

        Ny=floor(real(Ny_all)/n_procs)
        if (procid.eq.(n_procs-1)) Ny=Ny_all-(n_procs-1)*Ny
        yinit=-1+procid*Ny*dy
        if (procid.eq.(n_procs-1)) yinit=1-(Ny-1)*dy

        Nz=floor(real(Nz_all)/n_procs)
        if (procid.eq.(n_procs-1)) Nz=Nz_all-(n_procs-1)*Nz
        zinit=-1+procid*Nz*dz
        if (procid.eq.(n_procs-1)) zinit=1-(Nz-1)*dz

        return
        end
c--------------------------------------------------------------------------------------

c--------------------------------------------------------------------------------------
c Prints message depending on whether condition is true or false
c--------------------------------------------------------------------------------------

        subroutine print_message(msg_true, msg_false, condition)

            implicit none

            character(len=*) msg_true, msg_false
            logical condition
!----------------------------------------------------------------------
        
            if (condition) then
                write (*,*) msg_true
            else
                write (*,*) msg_false
            endif
        
        end
c--------------------------------------------------------------------------------------
