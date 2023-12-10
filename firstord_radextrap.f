c---------------------------------------------------------------------------------------------
c Given the grid point p1=i,j,k, the routine returns the coordinates and function value 
c at the point p_ex along the radial direction through p1 on a sphere at rho_ex
c
c Strategy
c 1. Gets coordinates and function value at p1
c 2. Gets coordinates of p_ex along radial direction through p1
c FROM HERE FORMULAE DEPEND ON LOCATION OF p1 ON THE GRID, SO CONSIDER VARIOUS CASES
c 3. Gets coordinates of 4 grid points a,b,c,d at the vertices of a square,
c     and the coordinates of p2 along the radial direction away from the sphere
c 4. Determines value of function at p2 via bilinear interpolation from 4 Cartesiang grid points a,b,c,d
c 5. Determines value of f at sphere point p_ex extrapolated using p1 and p2 (first order radial extrapolation)
c---------------------------------------------------------------------------------------------

        subroutine sphere_coords_func(
     &                      x_ex,y_ex,z_ex,
     &                      f_ex,  
     &                      i,j,k,             
     &                      f,
     &                      rho_ex,
     &                      x,y,z,Nx,Ny,Nz)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        real*8 x(Nx),y(Ny),z(Nz)
        integer N_ex
        real*8 x_ex,y_ex,z_ex

        integer i,j,k
        integer ip2a,jp2a,kp2a
        integer ip2b,jp2b,kp2b
        integer ip2c,jp2c,kp2c
        integer ip2d,jp2d,kp2d
        integer a,b,c,d

        real*8 x0,y0,z0,rho0,q,chi0,xi0
        real*8 xp2a,xp2b,xp2c,xp2d
        real*8 yp2a,yp2b,yp2c,yp2d
        real*8 zp2a,zp2b,zp2c,zp2d

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 f(Nx,Ny,Nz)
        real*8 f_p1
        real*8 f_ex

        real*8 f_p2a
        real*8 f_p2b
        real*8 f_p2c
        real*8 f_p2d
        real*8 f_p2

        real*8 xp1,yp1,zp1
        real*8 xp2,yp2,zp2
        real*8 rho_ex
        real*8 rhop1,chip1,xip1
        real*8 rhop2,chip2,xip2

        real*8 bilinear_interp
        real*8 firstord_extrap
!----------------------------------------------------------------------

!              call MPI_comm_rank(MPI_COMM_WORLD, rank, ierr)
!              call MPI_comm_size(MPI_COMM_WORLD, n_ranks, ierr)
!              if (rank.eq.1) then
!              do m=1,Nx
!                write (*,*) "rank,i,x(i),x(i-1)=",
!     &             rank,i,x(i),x(i-1)
!              end do
!              end if

              !coords and function at p1
              xp1=x(i)
              yp1=y(j)
              zp1=z(k)
              call cart_to_sph(rhop1,chip1,xip1,xp1,yp1,zp1)
              f_p1=f(i,j,k)

              !coords of p_ex
              call sph_to_cart(x_ex,y_ex,z_ex,rho_ex,chip1,xip1)

              !find coords of a,b,c,d and point p2 to interpolate
              !the formulae for these depend on the location of p1 on the Cartesian grid
              chip2=chip1
              xip2=xip1


              if ((abs(xp1).ge.abs(yp1)).and.
     &            (abs(xp1).ge.abs(zp1))) then !(i.e., |xp1|>=|yp1|,|zp1|, so xp1 cannot be 0)
            
                if (xp1.gt.0) then !(i.e., we are in the upper part of 
                                    !the spatial grid, called "a" sector)

                  ip2a=i-1
                  ip2b=i-1
                  ip2c=i-1
                  ip2d=i-1

                  xp2a=x(ip2a)
                  xp2b=x(ip2b)
                  xp2c=x(ip2c)
                  xp2d=x(ip2d)
                  xp2=x(ip2a)

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jp2a=j
                      jp2b=j
                      jp2c=j+1
                      jp2d=j+1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kp2a=k
                      kp2b=k
                      kp2c=k-1
                      kp2d=k-1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)


                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., yp1<=0: quadrant IIIa)

                        if ((abs(yp1).lt.10.0d-10).and.
     &                      (abs(zp1).lt.10.0d-10)) then
                          yp2a=0.0d0
                          yp2b=0.0d0
                          yp2c=0.0d0
                          yp2d=0.0d0

                          yp2=0.0d0
                          zp2=0.0d0
                          rhop2=abs(xp2)

                          f_p2=f(i-1,j,k)


                        else

                          jp2a=j
                          jp2b=j+1
                          jp2c=j+1
                          jp2d=j

                          yp2a=y(jp2a)
                          yp2b=y(jp2b)
                          yp2c=y(jp2c)
                          yp2d=y(jp2d)


                          yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                          zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                          rhop2=abs(xp2/cos(PI*chip2))

                          f_p2a=f(ip2a,jp2a,kp2a)
                          f_p2b=f(ip2b,jp2b,kp2b)
                          f_p2c=f(ip2c,jp2c,kp2c)
                          f_p2d=f(ip2d,jp2d,kp2d)

                          f_p2=
     &                        bilinear_interp(
     &                             f_p2a,
     &                             f_p2b,
     &                             f_p2c,
     &                             f_p2d,
     &                             yp2a,yp2b,yp2c,yp2d,
     &                             zp2a,zp2b,zp2c,zp2d,
     &                             yp2,
     &                             zp2)

                        end if

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                else !(i.e., xp1<0, i.e., we are in the lower part of 
                     !the spatial grid, called "b" sector)

                  ip2a=i+1
                  ip2b=i+1
                  ip2c=i+1
                  ip2d=i+1

                  xp2a=x(ip2a)
                  xp2b=x(ip2b)
                  xp2c=x(ip2c)
                  xp2d=x(ip2d)
                  xp2=x(ip2a)

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)
                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)
                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jp2a=j
                      jp2b=j
                      jp2c=j+1
                      jp2d=j+1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)
                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kp2a=k
                      kp2b=k
                      kp2c=k-1
                      kp2d=k-1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)


                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., yp1<=0: quadrant IIIa)

                        if ((abs(yp1).lt.10.0d-10).and.
     &                      (abs(zp1).lt.10.0d-10)) then
                          yp2a=0.0d0
                          yp2b=0.0d0
                          yp2c=0.0d0
                          yp2d=0.0d0

                          yp2=0.0d0
                          zp2=0.0d0
                          rhop2=abs(xp2)

                          f_p2=
     &                        f(i+1,j,k)


                        else

                          jp2a=j
                          jp2b=j+1
                          jp2c=j+1
                          jp2d=j

                          yp2a=y(jp2a)
                          yp2b=y(jp2b)
                          yp2c=y(jp2c)
                          yp2d=y(jp2d)


                          yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                          zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                          rhop2=abs(xp2/cos(PI*chip2))

                          f_p2a=f(ip2a,jp2a,kp2a)
                          f_p2b=f(ip2b,jp2b,kp2b)
                          f_p2c=f(ip2c,jp2c,kp2c)
                          f_p2d=f(ip2d,jp2d,kp2d)

                          f_p2=
     &                        bilinear_interp(
     &                             f_p2a,
     &                             f_p2b,
     &                             f_p2c,
     &                             f_p2d,
     &                             yp2a,yp2b,yp2c,yp2d,
     &                             zp2a,zp2b,zp2c,zp2d,
     &                             yp2,
     &                             zp2)

                        end if

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                end if !closes condition on (xp1.gt.0)


















              else if ((abs(yp1).ge.abs(zp1)).and.
     &                 (abs(yp1).ge.abs(xp1))) then !(i.e., |yp1|>=|zp1|,|xp1|, so yp1 cannot be 0)


                if (yp1.gt.0) then !(i.e., we are in either Q.Ia or Q.Ib or Q.IVa or Q.IVb)

                  jp2a=j-1
                  jp2b=j-1
                  jp2c=j-1
                  jp2d=j-1

                  yp2a=y(jp2a)
                  yp2b=y(jp2b)
                  yp2c=y(jp2c)
                  yp2d=y(jp2d)
                  yp2=y(jp2a)

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      kp2a=k
                      kp2b=k
                      kp2c=k-1
                      kp2d=k-1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., xp1<=0: quadrant Ib)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., zp1<=0, so Q.IVa or IVb)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., xp1<=0: quadrant IVb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1| - xp1=zp1=0 is not an issue)

                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IVb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      if (zp1.gt.0) then !(i.e., quadrant Ib)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., zp1<=0: quadrant IVb)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                else !(i.e., yp1<0, so we are in either Q.IIa or Q.IIb or Q.IIIa or Q.IIIb)

                  jp2a=j+1
                  jp2b=j+1
                  jp2c=j+1
                  jp2d=j+1

                  yp2a=y(jp2a)
                  yp2b=y(jp2b)
                  yp2c=y(jp2c)
                  yp2d=y(jp2d)
                  yp2=y(jp2a)

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant IIa or IIb)
                      kp2a=k
                      kp2b=k
                      kp2c=k-1
                      kp2d=k-1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)
                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., xp1<=0: quadrant IIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IIIb)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)
                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., xp1<=0: quadrant IIIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1| - the case xp1=zp1=0 is not an issue)

                    if (xp1.gt.0) then !(i.e., either quadrant IIa or IIIa)

                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., xp1<=0, so Q.IIb or IIIb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      if (zp1.gt.0) then !(i.e., quadrant IIb)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., zp1<=0: quadrant IIIb)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)


                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                end if !closes condition on (yp1.gt.0)











              else !(i.e., |zp1|>=|xp1|,|yp1|, so zp1 cannot be 0)

                if (zp1.gt.0) then !(i.e., we are in either Q.Ia or Q.Ib or Q.IIa or Q.IIb)

                  kp2a=k-1
                  kp2b=k-1
                  kp2c=k-1
                  kp2d=k-1

                  zp2a=z(kp2a)
                  zp2b=z(kp2b)
                  zp2c=z(kp2c)
                  zp2d=z(kp2d)
                  zp2=z(kp2a)

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IIb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      if (yp1.gt.0) then !(i.e., quadrant Ib)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., yp1<=0: quadrant IIb)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1| - yp1=xp1=0 is not an issue)

                    if (yp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., xp1<=0: quadrant Ib)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIb)

                      jp2a=j
                      jp2b=j
                      jp2c=j+1
                      jp2d=j+1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., xp1<=0: quadrant IIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))

                else !(i.e., zp1<0, so we are in either Q.IIIa or Q.IIIb or Q.IVa or Q.IVb)

                  kp2a=k+1
                  kp2b=k+1
                  kp2c=k+1
                  kp2d=k+1

                  zp2a=z(kp2a)
                  zp2b=z(kp2b)
                  zp2c=z(kp2c)
                  zp2d=z(kp2d)
                  zp2=z(kp2a)

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant IIIa or IVa)
                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)
                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., yp1<=0: quadrant IIIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., xp1<=0, so Q.IIIb or IVb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      if (yp1.gt.0) then !(i.e., quadrant IVb)
                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      else !(i.e., yp1<=0: quadrant IIIb)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1| - yp1=xp1=0 is not an issue)

                    if (yp1.gt.0) then !(i.e., either quadrant IVa or IVb)

                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., xp1<=0: quadrant IVb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)

                      end if

                    else !(i.e., yp1<=0, so Q.IIIa or IIIb)

                      jp2a=j
                      jp2b=j
                      jp2c=j+1
                      jp2d=j+1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        f_p2a=f(ip2a,jp2a,kp2a)
                        f_p2b=f(ip2b,jp2b,kp2b)
                        f_p2c=f(ip2c,jp2c,kp2c)
                        f_p2d=f(ip2d,jp2d,kp2d)

                        f_p2=
     &                      bilinear_interp(
     &                           f_p2a,
     &                           f_p2b,
     &                           f_p2c,
     &                           f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        f_ex=
     &                      firstord_extrap(
     &                           f_p1,
     &                           f_p2,
     &                           rhop1,rhop2,rho_ex)


                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))

                end if !closes condition on (zp1.gt.0)



              end if !closes condition on ((abs(xp1).gt.abs(yp1)).and.(abs(xp1).gt.abs(zp1)))

!                if (abs(rho_ex-f_ex).ge.0.03125d0) then
!                    write(*,*) "TEST FAILED"
!                    write(*,*) "at xp1,yp1,zp1,
!     &                 xp2a,yp2a,zp2a,
!     &                 xp2b,yp2b,zp2b,
!     &                 xp2c,yp2c,zp2c,
!     &                 xp2d,yp2d,zp2d,
!     &                 xp2,yp2,zp2,               
!     &                 x_ex,y_ex,z_ex,f_ex=",
!     &                xp1,yp1,zp1,
!     &                 xp2a,yp2a,zp2a,
!     &                 xp2b,yp2b,zp2b,
!     &                 xp2c,yp2c,zp2c,
!     &                 xp2d,yp2d,zp2d,
!     &                 xp2,yp2,zp2,
!     &                x_ex,y_ex,z_ex,
!     &               f_ex
!                end if

        return
        end
c--------------------------------------------------------------------------------------