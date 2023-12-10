c---------------------------------------------------------------------------------------------
c On a 3-dimensional Cartesian grid in Cartesian coordinates x,y,z in [-1,1]
c with Nx/Ny/Nz points along the x/y/z-direction 
c
c Perform first order radial extrapolation to a sphere with radius rho_ex, i.e.,
c 1. identify the points closest to the sphere of radius rho_ex via a mask,
c 2. for each of these points p1, identify the value of the function f_p1
c and the value of the function f_p2 at the point p2 
c along the radial direction away from the sphere (done via bilinear interpolation)
c 3. extrapolate the value of f, called f_ex, at the sphere point p_ex along the radial direction using p1 and p2
c---------------------------------------------------------------------------------------------

c--------------------------------------------------------------------------------------
    	program main

!----------------------------------------------------------------------

        use mpi
        implicit none

        integer rank, n_ranks, ierr

        integer Nx_all,Ny_all,Nz_all
        parameter (Nx_all=17,Ny_all=17,Nz_all=17)
        real*8 rho_ex
        parameter (rho_ex=1.0d0)
        integer tag
        parameter (tag=1)
        character*40 output_nametag
        parameter (output_nametag="output_index_xex_yex_zex_fex_rank")
        character*50 output_filename

        !specifies the function to extrapolate
        integer f_id
        parameter (f_id=0)

        real*8 dx,dy,dz
        integer Nx,Ny,Nz
        real*8 xinit,yinit,zinit
        logical test_passed
        logical printed

!----------------------------------------------------------------------

                call MPI_init(ierr)
                call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
                call MPI_comm_size(MPI_COMM_WORLD,n_ranks,ierr)

                !print chosen parameters
                if (rank.eq.0) call print_params_toscreen(rho_ex,
     &             Nx_all,Ny_all,Nz_all,f_id,output_nametag,n_ranks)

                write(output_filename, '(A, I0, A)') 
     &             trim(output_nametag), rank, ".txt"

                call setup_local_grid(
     &                  Nx,Ny,Nz,
     &                  dx,dy,dz,
     &                  xinit,yinit,zinit,
     &                  Nx_all,Ny_all,Nz_all,
     &                  n_ranks,rank)

                call local_main(
     &             Nx,Ny,Nz,dx,dy,dz,
     &             xinit,yinit,zinit,
     &             test_passed,printed,
     &             rho_ex,tag,output_filename,f_id)

                !check result of all local tests
                call MPI_ALLREDUCE(test_passed,test_passed,1,
     &             MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
                if (rank.eq.0) then
                 if (f_id.eq.0) then
                  call print_message(
     &             "Test passed - sphere values are close to rho_ex",
     &             "Test failed", test_passed)
                 else
                  write(*,*) 
     &             "No test for the chosen test function"
                 end if
                end if

                !check that all processes printed
                call MPI_ALLREDUCE(printed,printed,1,
     &             MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)

                if (rank.eq.0) then
                  call print_message("Printing successful",
     &               "Printing failed", printed)
                end if

                call MPI_finalize(ierr)

    	end program main
c--------------------------------------------------------------------------------------

c--------------------------------------------------------------------------------------
c Calls the main routines for a single process
c--------------------------------------------------------------------------------------
        subroutine local_main(Nx,Ny,Nz,dx,dy,dz,
     &       xinit,yinit,zinit,
     &       test_passed,printed,
     &       rho_ex,tag,output_filename,f_id)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        real*8 dx,dy,dz

        real*8 rho_ex
        integer tag
        character*40 output_filename
        integer f_id

        real*8 xinit,yinit,zinit
        real*8 x(Nx),y(Ny),z(Nz)
        integer mask(Nx,Ny,Nz)
        real*8 f(Nx,Ny,Nz)
        integer N_ex
        logical test_passed
        logical printed

!----------------------------------------------------------------------

        !initialise Cartesian coordinates, function to extrapolate, and mask that determines points p1
        call initialise_coords_func_mask(f,f_id,x,y,z,
     &         xinit,yinit,zinit,
     &         rho_ex,mask,tag,Nx,Ny,Nz,dx,dy,dz,N_ex)


        !perform radial extrapolation, print to file and test
        call radextrap_printtofile_test(f,f_id,x,y,z,test_passed,
     &   printed,rho_ex,mask,tag,Nx,Ny,Nz,dx,dy,dz,N_ex,output_filename)


        return
        end
c--------------------------------------------------------------------------------------
