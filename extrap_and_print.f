c--------------------------------------------------------------------------
c Calls radial extrapolation routine and print to file routine
c Runs a test of the extrapolated values 
c--------------------------------------------------------------------------
        subroutine radextrap_printtofile_test(f,num_f,x,y,z,
     &              rho_ex,mask,tag,Nx,Ny,Nz,N_ex,filename)

        implicit none

        integer Nx,Ny,Nz
        integer i,j,k
        integer N_ex
        integer ind_ex
        real*8 x(Nx),y(Ny),z(Nz)
        integer mask(Nx,Ny,Nz),tag
        real*8 f(Nx,Ny,Nz)
        real*8 x_ex(N_ex),y_ex(N_ex),z_ex(N_ex)
        real*8 dx
        real*8 f_ex(N_ex)
        real*8 rho_ex

        character*40 filename

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer num_f

!--------------------------------------------------------------------

        write (*,*) "Extrapolating..."
        call radextrap(f,x,y,z,rho_ex,mask,tag,Nx,Ny,Nz,N_ex
     &                ,x_ex,y_ex,z_ex,f_ex)
        call separation_line()
        call print_to_file(x_ex,y_ex,z_ex,f_ex,N_ex,filename)
        call separation_line()


        ! implement test only for the case num_f=0, i.e., f=rho:=sqrt(x**2+y**2+z**2)
        if (num_f.eq.0) then
            write(*,*) "TEST for f=rho:=sqrt(x**2+y**2+z**2)"
            write(*,*) "Sphere values must be close to rho_ex=",rho_ex
            dx=2.0d0/(Nx-1)
            do ind_ex=1,N_ex
                if (abs(rho_ex-f_ex(ind_ex)).ge.dx) then
                    write(*,*) "TEST FAILED"
                    write(*,*) "at ind_ex=",ind_ex
                    stop
                end if
            end do
            write(*,*) "TEST PASSED"
        else
            write(*,*) "No test for the chosen test function"
        end if

        return
        end
c--------------------------------------------------------------------------------------



c------------------------------------------------------------------------------------------------------------------------------------
c Loop over grid and, for each grid point (p1) such that mask(i,j,k)=tag, 
c get coordinates and function value of sphere point p_ex along the radial direction through p1 via radial extrapolation
c------------------------------------------------------------------------------------------------------------------------------------
        subroutine radextrap(f,x,y,z,rho_ex,mask,tag,Nx,Ny,Nz,N_ex,
     &                      x_ex,y_ex,z_ex,f_ex)

        implicit none

        integer Nx,Ny,Nz
        integer i,j,k
        integer N_ex
        integer ind_ex
        real*8 x(Nx),y(Ny),z(Nz)
        integer mask(Nx,Ny,Nz),tag
        real*8 f(Nx,Ny,Nz)
        real*8 x_ex(N_ex),y_ex(N_ex),z_ex(N_ex)
        real*8 f_ex(N_ex)
        real*8 rho_ex

        real*8 PI
        parameter (PI=3.141592653589793d0)

!--------------------------------------------------------------------

        ind_ex=0
        do i=1,Nx
            do j=1,Ny
                do k=1,Nz
                    if (mask(i,j,k).eq.tag) then
                        !index identifying sphere points
                        ind_ex=ind_ex+1

                        !get coordinates and function value of p_ex along radial direction through p1
                        call sphere_coords_func(
     &                      x_ex(ind_ex),y_ex(ind_ex),z_ex(ind_ex),
     &                      f_ex(ind_ex),    
     &                      i,j,k,           
     &                      f,
     &                      rho_ex,
     &                      x,y,z,Nx,Ny,Nz)

                    end if
                end do
            end do
        end do

        return
        end
c--------------------------------------------------------------------------------------


c--------------------------------------------------------------------------
c Print to file
c--------------------------------------------------------------------------
        subroutine print_to_file(x,y,z,f,N,filename)

        implicit none

        integer N
        integer ind
        real*8 x(N),y(N),z(N)
        real*8 f(N)

        character*40 filename

        integer unitnum
        parameter (unitnum=1)

        integer ierr

!--------------------------------------------------------------------

        write(*,*) "Printing to ", filename
        open(unitnum, file=filename, status="UNKNOWN", iostat=ierr)
        if (ierr.ne.0) then
            write(*,*) "Error opening ", filename
            stop
        end if
        do ind=1,N
            write(unitnum, *)    ind, 
     &                       " ",x(ind),
     &                       " ",y(ind),
     &                       " ",z(ind),
     &                       " ",f(ind)
        end do
        close(unitnum)

        return
        end
c--------------------------------------------------------------------------------------