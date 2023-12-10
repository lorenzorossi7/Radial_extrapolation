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
    	program main

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        parameter (Nx=257,Ny=257,Nz=257)
        real*8 rho_ex
        parameter (rho_ex=1.0d0)
        integer tag
        parameter (tag=1)
        character*40 output_filename
        parameter (output_filename="output_index_xex_yex_zex_fex.txt")

        !specifies the function to extrapolate
        integer num_f
        parameter (num_f=0)

        real*8 x(Nx),y(Ny),z(Nz)
        integer mask(Nx,Ny,Nz)
        real*8 f(Nx,Ny,Nz)
        integer N_ex

!----------------------------------------------------------------------

                write (*,*) "Starting..."
                call separation_line()

                !print chosen parameters
                call print_params(rho_ex,Nx,Ny,Nz,output_filename)

                call separation_line()
		!initialise Cartesian coordinates, function to extrapolate, and mask that determines points p1
		call initialise_grid_func_mask(f,num_f,x,y,z,
     &                            rho_ex,mask,tag,Nx,Ny,Nz,N_ex)

                call separation_line()
		!perform radial extrapolation, print to file and test
		call radextrap_printtofile_test(f,num_f,x,y,z,rho_ex,mask,tag
     &				,Nx,Ny,Nz,N_ex,output_filename)
                call separation_line()

                write (*,*) "Finished"


    	end
c--------------------------------------------------------------------------------------
