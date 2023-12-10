c--------------------------------------------------------------------------
c Calls radial extrapolation routine and print to file routine
c Runs a test of the extrapolated values 
c--------------------------------------------------------------------------
        subroutine radextrap_printtofile_test(f,f_id,x,y,z,test_passed,
     &         printed,rho_ex,mask,tag,Nx,Ny,Nz,dx,dy,dz,N_ex,filename)

        implicit none

        integer Nx,Ny,Nz
        integer i,j,k
        integer N_ex
        integer ind_ex
        real*8 x(Nx),y(Ny),z(Nz)
        integer mask(Nx,Ny,Nz),tag
        real*8 f(Nx,Ny,Nz)
        real*8 x_ex(N_ex),y_ex(N_ex),z_ex(N_ex)
        real*8 dx,dy,dz
        real*8 f_ex(N_ex)
        real*8 rho_ex

        character*40 filename

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer f_id
        logical test_passed,test_feqrho
        logical printed, print_to_file

!--------------------------------------------------------------------

        !perform radial extrapolation
        call radextrap(x_ex,y_ex,z_ex,f_ex
     &          ,f,x,y,z,rho_ex,mask,tag,Nx,Ny,Nz,N_ex)

        !run test
        if (f_id.eq.0)
     &      test_passed=
     &       test_feqrho(x_ex,y_ex,z_ex,f_ex,N_ex,dx,dy,dz,rho_ex)

        !print to file
        printed=print_to_file(x_ex,y_ex,z_ex,f_ex,N_ex,filename)

        return
        end
c--------------------------------------------------------------------------------------



c------------------------------------------------------------------------------------------------------------------------------------
c Loop over grid and, for each grid point (p1) such that mask(i,j,k)=tag, 
c get coordinates and function value of sphere point p_ex along the radial direction through p1 via radial extrapolation
c------------------------------------------------------------------------------------------------------------------------------------
        subroutine radextrap(x_ex,y_ex,z_ex,f_ex
     &                ,f,x,y,z,rho_ex,mask,tag,Nx,Ny,Nz,N_ex)

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

        real*8 xp1,yp1,zp1
        real*8 rhop1,chip1,xip1
        real*8 xex,yex,zex

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
c Run a test - implemented only for the case f_id=0, i.e., f=rho:=sqrt(x**2+y**2+z**2)
c--------------------------------------------------------------------------
        logical function test_feqrho(x,y,z,f,N,dx,dy,dz,rho_ex)

        implicit none

        integer f_id
        integer proc
        integer N
        integer ind
        real*8 x(N),y(N),z(N)
        real*8 f(N)
        real*8 dx,dy,dz
        real*8 rho_ex

!--------------------------------------------------------------------
        test_feqrho=.true.
        do ind=1,N
            if (abs(rho_ex-f(ind)).ge.dx) then
                test_feqrho=.false.
                return
            end if
        end do

        return
        end
c--------------------------------------------------------------------------------------

c--------------------------------------------------------------------------
c Print to file
c--------------------------------------------------------------------------
        logical function print_to_file(x,y,z,f,N,filename)

        use mpi
        implicit none

        integer rank, n_ranks

        integer N
        integer ind
        real*8 x(N),y(N),z(N)
        real*8 f(N)

        character*40 filename

        integer unitnum
        parameter (unitnum=1)

        integer ierr

!--------------------------------------------------------------------

        call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
        call MPI_comm_size(MPI_COMM_WORLD,n_ranks,ierr)
        open(unitnum, file=filename, status="UNKNOWN", iostat=ierr)
        if (ierr.ne.0) then
            print_to_file=.false.
            return
        end if
        do ind=1,N
            write(unitnum, *)    x(ind),
     &                       " ",y(ind),
     &                       " ",z(ind),
     &                       " ",f(ind)
        end do
        close(unitnum)
        print_to_file=.true.

        return
        end
c--------------------------------------------------------------------------------------