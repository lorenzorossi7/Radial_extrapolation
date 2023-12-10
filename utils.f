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
c Prints parameters
c----------------------------------------------------------------------
        subroutine print_params(rho_ex,Nx,Ny,Nz,filename)
        implicit none
        real*8 rho_ex
        integer Nx,Ny,Nz

        character*40 filename

        !--------------------------------------------------------------

        write (*,*) "Perform radial extrapolation at rho_ex=",rho_ex
        write (*,*) "with resolution Nx,Ny,Nz=",Nx,Ny,Nz
        write (*,*) "and grid spacing dx,dy,dz=",
     &         2.0d0/(Nx-1),2.0d0/(Ny-1),2.0d0/(Nz-1)
        write (*,*) "Output file is ",filename

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c Prints separation line
c----------------------------------------------------------------------
        subroutine separation_line()

        !--------------------------------------------------------------

        
        write (*,*) "======================================="

        end
c--------------------------------------------------------------------------------------