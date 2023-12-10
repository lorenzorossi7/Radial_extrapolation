c--------------------------------------------------------------------------
c Initialises
c 1. values of Cartesian coordinates x,y,z on the grid
c 2. points p1 closest to the sphere to use for extrapolation (using mask)
c 3. number of extrapolated points on the sphere, N_ex
c--------------------------------------------------------------------------
         subroutine initialise_grid_func_mask(f,num_f,x,y,z,
     &                         rho_ex,mask,tag,
     &                         Nx,Ny,Nz,N_ex)
 
         implicit none
 
         integer Nx,Ny,Nz
         integer i,j,k
         integer N_ex
         integer ind_ex
         real*8 x(Nx),y(Ny),z(Nz)
         integer mask(Nx,Ny,Nz),tag
         real*8 rho,chi,xi
         real*8 dx,dy,dz
         real*8 rho_ex
         real*8 f(Nx,Ny,Nz)
         integer num_f

!--------------------------------------------------------------------

         write(*,*) "Initialising Cartesian coordinates and"
         if (num_f.eq.0) then
            write(*,*) "Initialising f=rho:=sqrt(x**2+y**2+z**2) and"
         else
            write (*,*) "ERROR: num_f>0 is not implemented."
            stop
         end if
         write(*,*) "Initialising mask function"

         dx=2.0d0/(Nx-1)
         dy=2.0d0/(Ny-1)
         dz=2.0d0/(Nz-1)
         N_ex=0
         do i=1,Nx
            x(i)=-1+(i-1)*dx
            do j=1,Ny
               y(j)=-1+(j-1)*dy
               do k=1,Nz
                  z(k)=-1+(k-1)*dz

                  !get spherical coordinate values of point (i,j,k)
                  call cart_to_sph(rho,chi,xi,x(i),y(j),z(k))
 
                  !specify function to extrapolate
                  f(i,j,k)=rho
 
                  !set mask to tag closest points to the sphere of radius rho_ex
                  if (rho.lt.rho_ex.and.rho.gt.(rho_ex-3.0d0*dx/2)) then
                     mask(i,j,k)=tag
                     N_ex=N_ex+1
                  else
                     mask(i,j,k)=1-tag
                  end if
               end do
            end do
         end do
 
         return
         end
c--------------------------------------------------------------------------------------