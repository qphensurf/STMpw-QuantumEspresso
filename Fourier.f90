  module Fourier
!
! In order to make these subroutines work
! please link with lfftw
!
CONTAINS

!  subroutine matching (CW, A_G, C_G, ngx, ngy, ngz, z_s, z_t,bg,xp,NPL,enmax,temp)  
   subroutine matching (CW, A_G, C_G, ngx, ngy, ngz, z_s, z_t,temp)


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine performs the matching of the coefficients
! A_G at z_s of the FFT(CW)
! C_G at z_t of the FFT(CW)
! A_G and C_G are returned in reciprocal space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Use volumen
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer, parameter :: qs = SELECTED_REAL_KIND(5)
   integer, intent (in):: ngx, ngy, ngz
   integer ::  NPL
   integer :: rangex (ngx), rangey (ngy), rangez (ngz)
   real (q), intent (in) :: z_s, z_t
   real (q) :: g1, g2, g3, gix, giy, giz, etotal, enmax
   real (q):: xp (3), bg(3,3)
   complex (q), intent (in) :: CW(ngx*ngy*ngz)
   complex (q), intent (out) :: A_G(0:ngx-1,0:ngy-1)
   complex (q), intent (out) :: C_G(0:ngx-1,0:ngy-1)
   complex (q) :: temp(0:ngx-1,0:ngy-1,0:ngz-1)
   integer :: ix, iy, iz ,l, i_s, i_t, test


   temp = (0,0)
   A_G = (0,0)
   C_G = (0,0)

! Cutoff sphere
!     call inilpc (ngx, ngy, ngz, rangex, rangey, rangez)


!      l=0

!      do iz = 1, ngz
!       do iy = 1, ngy
!        do ix = 1, ngx
!         g1 = rangex(ix) + xp (1)
!         g2 = rangey(iy) + xp (2)
!         g3 = rangez(iz) + xp (3)

!         gix = g1*bg(1,1) + g2*bg(2,1) + g3*bg(3,1) 
!         giy = g1*bg(1,2) + g2*bg(2,2) + g3*bg(3,2) 
!         giz = g1*bg(1,3) + g2*bg(2,3) + g3*bg(3,3) 

!         etotal = (gix**2 + giy**2 + giz**2)/2

!         if(etotal < enmax/2) then
!            l=l+1
!         temp (ix-1,iy-1,iz-1) = CW(l)
!         endif
!        enddo
!       enddo
!      enddo

!  if (l /= NPL) then
!        print *,'PROBLEM with cutoff, wf cutoff does not agree with the read one!'
!        print *, 'NPL=',NPL,'Lmax=', l
!  endif

   l=0
   do iz = 1, ngz
   do iy = 1, ngy
   do ix = 1, ngx
   l=l+1
   temp (ix-1,iy-1,iz-1) = CW(l)
   enddo
   enddo
   enddo

   call fft1d_r (temp, ngx, ngy, ngz)



! find index for the matching

   i_s = mod (int(z_s*ngz), ngz)
   i_t = mod (int(z_t*ngz), ngz)


! matching at the corresponding plane

   A_G(0:ngx-1,0:ngy-1) = temp (0:ngx-1,0:ngy-1,i_s-1) ! sample
   C_G(0:ngx-1,0:ngy-1) = temp (0:ngx-1,0:ngy-1,i_t-1) ! tip


   return
   end subroutine matching 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1-D fft BACKWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft1d_r (temp, ngx, ngy, ngz)
   use declfft
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer i_x, i_y,ngx, ngy, ngz
   complex (q) matrix (0:ngz-1)
   complex (q) temp (0:ngx-1,0:ngy-1,0:ngz-1)

   do i_x = 0, ngx-1
   do i_y = 0, ngy-1
   matrix(:) = temp (i_x,i_y,:)
   call dfftw_plan_dft_1d (plan,ngz,matrix,matrix,FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   temp (i_x,i_y,:) = matrix(:)

   enddo
   enddo

   end subroutine fft1d_r


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3-D fft FORWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft3d (matrix, ngx, ngy, ngz)
   use declfft
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer ngx, ngy, ngz
   complex (q) matrix (0:ngx-1,0:ngy-1,0:ngz-1)

   call dfftw_plan_dft_3d (plan,ngx,ngy,ngz,matrix,matrix,FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   end subroutine fft3d 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3-D fft BACKWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft3d_r (matrix, ngx, ngy, ngz)
   use declfft
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer ngx, ngy, ngz
   complex (q) matrix (0:ngx-1,0:ngy-1,0:ngz-1)

   call dfftw_plan_dft_3d (plan,ngx,ngy,ngz,matrix,matrix,FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   end subroutine fft3d_r
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2-D fft FORWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft2d (matrix, ngx, ngy)
   use declfft
   implicit none
   integer ngx, ngy
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   complex (q) matrix (0:ngx-1,0:ngy-1)

   call dfftw_plan_dft_2d (plan,ngx,ngy,matrix,matrix,FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   end subroutine fft2d 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2-D fft BACKWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft2d_r (matrix, ngx, ngy)
   use declfft
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer ngx, ngy
   complex (q) matrix (0:ngx-1,0:ngy-1)

   call dfftw_plan_dft_2d (plan,ngx,ngy,matrix,matrix,FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)
   end subroutine fft2d_r 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine density (density1,CW,ngx, ngy, ngz, IGX,IGY,IGZ,NPL)
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer, parameter :: qs = SELECTED_REAL_KIND(5)
   integer, intent (in):: ngx, ngy, ngz, NPL
   integer :: IGX (NPL), IGY (NPL), IGZ (NPL)
   complex (qs), intent (in) :: CW(ngx*ngy*ngz)
   real (q), intent (out) :: density1(ngz)
   complex (q) :: temp(0:ngx-1,0:ngy-1,0:ngz-1)
   integer :: ix, iy, iz ,l, i_s, i_t,test

   temp = (0,0)
   density1 = 0
   do l = 1, NPL
   temp (IGX(l)-1,IGY(l)-1,IGZ(l)-1) = CW(l)*1.0_q
   enddo


   call fft3d_r (temp, ngx, ngy, ngz)

   do iz =1, ngz
   do ix = 1, ngx
   do iy = 1, ngy
   density1(iz) = density1(iz) + real(temp (ix-1,iy-1,iz-1)*conjg(temp (ix-1,iy-1,iz-1)))
   enddo
   enddo
   enddo

   density1 = density1/ (ngx*ngy)
   end subroutine density
!
! assign ranges
!
      subroutine inilpc (ngx,ngy,ngz,rangex,rangey,rangez)
      implicit none
      integer ngx, ngy, ngz
      integer rangeX(ngx),rangeY(ngy),rangez(ngz)
      integer ix, iy, iz

      do ix = 1,(ngx/2)+1
        rangex(ix) = ix-1
      enddo

      do ix = (ngx/2)+2,ngx
        rangex(ix) = ix-1-ngx
      enddo

      do iy = 1,(ngy/2)+1
        rangey(iy) = iy-1
      enddo

      do iy = (ngy/2)+2,ngy
        rangey(iy) = iy-1-ngy
      enddo

      do iz = 1,(ngz/2)+1
        rangez(iz) = iz-1
      enddo

      do iz = (ngz/2)+2,ngz
        rangez(iz) = iz-1-ngz
      enddo

      return

      end subroutine inilpc


  end module Fourier
