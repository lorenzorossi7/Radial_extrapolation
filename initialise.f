c--------------------------------------------------------------------------
c Initialises
c 1. values of Cartesian coordinates x,y,z on the grid
c 2. points p1 closest to the sphere to use for extrapolation (using mask)
c 3. number of extrapolated points on the sphere, N_ex
c--------------------------------------------------------------------------
         subroutine initialise_coords_func_mask(f,f_id,x,y,z,
     &                         xinit,yinit,zinit,    
     &                         rho_ex,mask,tag,
     &                         Nx,Ny,Nz,dx,dy,dz,N_ex)

         implicit none
 
         integer Nx,Ny,Nz
         real*8 dx,dy,dz
         integer i,j,k
         integer N_ex

         real*8 xinit,yinit,zinit
         integer ind_ex
         real*8 x(Nx),y(Ny),z(Nz)
         integer mask(Nx,Ny,Nz),tag
         real*8 rho,chi,xi
         real*8 rho_ex
         real*8 f(Nx,Ny,Nz)
         integer f_id

!--------------------------------------------------------------------

         N_ex=0
         do i=1,Nx
            x(i)=xinit+(i-1)*dx
            do j=1,Ny
               y(j)=yinit+(j-1)*dy
               do k=1,Nz
                  z(k)=zinit+(k-1)*dz

                  !get spherical coordinate values of point (i,j,k)
                  call cart_to_sph(rho,chi,xi,x(i),y(j),z(k))
 
                  !specify function to extrapolate
                  f(i,j,k)=rho
 
                  !set mask to tag closest points to the sphere of radius rho_ex
                  if ((rho.lt.rho_ex.and.rho.gt.(rho_ex-3.0d0*dx/2))
     &                .and.((i-1).ge.1).and.((i+1).le.Nx).and.  
     &                ((j-1).ge.1).and.((j+1).le.Ny).and.  
     &                ((k-1).ge.1).and.((k+1).le.Nz)
     &                   ) then
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