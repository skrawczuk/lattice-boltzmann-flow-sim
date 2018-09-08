subroutine step(f, height, width, omega, u0, w, rho, ux, uy, u, sim_object)
    ! Calculates the new number density, density, and velocities for one 
    ! timestep of the lattice-Boltzmann simulation. Subroutine is wrapped  
    ! in python using f2py. 
    
    implicit none
    
    ! input parameters
    integer                           :: height, width
    integer, dimension(height, width) :: sim_object
    real, dimension(height, width, 9) :: f
    real                              :: omega, u0
    real, dimension(9)                :: w
    
    ! local parameters 
    integer                           :: i, j
    real, dimension(height, width)    :: rho, ux, uy, u
    real, dimension(9, height, width) :: eu 
    real, dimension(height, width, 9) :: f_copy, f_eq
    
    ! f2py parameters
    !f2py integer, intent(in     )    :: height, width
    !f2py real,    intent(in     )    :: omega, u0
    !f2py real(9), intent(in     )    :: w
    !f2py real,    intent(in     ), dimension(height, width), depend(height, width)    :: sim_object
    !f2py real,    intent(in, out), dimension(height, width, 9), depend(height, width) :: f
    !f2py real,    intent(    out), dimension(height, width), depend(height, width)    :: rho, ux, uy, u
    
    
    ! ==========================================================================
    ! collision step: calculating macroscopic quantities
    ! ==========================================================================

    ! updating density 
    rho = sum(f, dim=3)
    
    ! updating velocities
    ux = (f(:,:,3) + f(:,:,6) + f(:,:,9) - (f(:,:,1) + f(:,:,4) + f(:,:,7))) / rho
    uy = (f(:,:,1) + f(:,:,2) + f(:,:,3) - (f(:,:,7) + f(:,:,8) + f(:,:,9))) / rho
    u = sqrt(ux**2 + uy**2)
    
    ! updating velocity vectors
    eu(1,:,:) = uy-ux 
    eu(2,:,:) = uy
    eu(3,:,:) = ux+uy
    eu(4,:,:) = -ux
    eu(5,:,:) = 0
    eu(6,:,:) = ux
    eu(7,:,:) = -uy-ux
    eu(8,:,:) = -uy
    eu(9,:,:) = -uy+ux     

    ! updating equilibrium number density 
    !$omp parallel do
    do i=1, 9 
        f_eq(:,:,i) = rho * w(i) * (1 + 3*eu(i,:,:) + 4.5*eu(i,:,:)**2 - 1.5*u**2) 
    end do
    !$omp end parallel do
    
    f = f + omega * (f_eq - f)
    
    ! setting steady flow to the left 
    f(:,width,3) = w(3) * (1. + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    f(:,width,6) = w(6) * (1. + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    f(:,width,9) = w(9) * (1. + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    f(:,width,1) = w(1) * (1. - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    f(:,width,4) = w(4) * (1. - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    f(:,width,7) = w(7) * (1. - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    
    f_copy = f
    
    
    ! ==========================================================================
    ! streaming step: moving particles to neighboring lattice sites
    ! ==========================================================================
    
    !$omp parallel do
    do i=1,height
        ! all non-boundary streaming 
        do j=2,width-1
            f(i,j,1) = f_copy(min(i+1, height), j-1, 1)
            f(i,j,2) = f_copy(min(i+1, height), j, 2)
            f(i,j,3) = f_copy(min(i+1, height), j+1, 3)
            f(i,j,4) = f_copy(i, j-1, 4)
            f(i,j,5) = f_copy(i, j, 5)
            f(i,j,6) = f_copy(i, j+1, 6)
            f(i,j,7) = f_copy(max(i-1, 1), j-1, 7)
            f(i,j,8) = f_copy(max(i-1, 1), j, 8)
            f(i,j,9) = f_copy(max(i-1, 1), j+1, 9)               
        end do
        
        ! left and right boundary conditions
        f(i,1,1) = f_copy(min(i+1, height), width, 1)
        f(i,1,2) = f_copy(min(i+1, height), 1, 2)
        f(i,1,3) = f_copy(min(i+1, height), 2, 3)
        f(i,1,4) = f_copy(i, width, 4)
        f(i,1,5) = f_copy(i, 1, 5)
        f(i,1,6) = f_copy(i, 2, 6)
        f(i,1,7) = f_copy(max(i-1, 1), width, 7)
        f(i,1,8) = f_copy(max(i-1, 1), 1, 8)
        f(i,1,9) = f_copy(max(i-1, 1), 2, 9)
    
        f(i,width,1) = f_copy(min(i+1, height), width-1, 1)
        f(i,width,2) = f_copy(min(i+1, height), width, 2)
        f(i,width,3) = f_copy(min(i+1, height), 1, 3)
        f(i,width,4) = f_copy(i, width-1, 4)
        f(i,width,5) = f_copy(i, width, 5)
        f(i,width,6) = f_copy(i, 1, 6)
        f(i,width,7) = f_copy(max(i-1, 1), width-1, 7)
        f(i,width,8) = f_copy(max(i-1, 1), width, 8)
        f(i,width,9) = f_copy(max(i-1, 1), 1, 9)
    end do
    !$omp end parallel do
    
    f_copy = f
    
    ! top and bottom boundary conditions
    f(1,:,7) = f(1,:,1)
    f(1,:,8) = f(1,:,2)
    f(1,:,9) = f(1,:,3)
    
    f(height,:,3) = f_copy(height,:,9)
    f(height,:,2) = f_copy(height,:,8)
    f(height,:,1) = f_copy(height,:,7)
    
    ! object boundary conditions
    !$omp parallel do
    do i=1, height  
        do j=1, width
            if (sim_object(i,j) == 1) then
                f(i,j,1) = f(i,j,9)
                f(i,j,4) = f(i,j,6)
                f(i,j,7) = f(i,j,3)
                f(i,j,2) = f(i,j,8)
                
                f(i,j,3) = f_copy(i,j,7)
                f(i,j,6) = f_copy(i,j,4)
                f(i,j,9) = f_copy(i,j,1)
                f(i,j,8) = f_copy(i,j,2)
            end if
        end do
    end do
    !$omp end parallel do

end subroutine step