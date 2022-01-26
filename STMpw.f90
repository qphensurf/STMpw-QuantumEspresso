program STMpw
!
! Version for QUANTUM ESPRESSO 6.6
!
! Simple implementation for Bardeen's transfer hamiltonian
! and tunneling extension for the STM
!
! Gnu licence, check out copy of the LICENCE in this directory.
!
!
! Coded by N. Lorente and R. Robles
!                    San Sebastian, year 2020
! Obtainable from GitHub
!
!STMpw modules
  use declaration
  use volumen
  use Fourier
  use currentBardeen
!QE modules

  USE kinds, ONLY : DP
  USE io_files,  ONLY : prefix, tmp_dir, diropn, restart_dir
  USE wvfct,     ONLY : nbnd, npwx, et, wg
  USE klist,     ONLY : xk, nks, ngk, igk_k, wk
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE mp_global, ONLY : mp_startup
  USE mp_images, ONLY : intra_image_comm
  USE mp_pools,  ONLY : npool
  USE wavefunctions, ONLY : psic, evc
  USE io_files,  ONLY : nwordwfc, iunwfc
  USE gvect,     ONLY : ngm, g
  USE noncollin_module, ONLY : npol, noncolin
  USE environment, ONLY : environment_start, environment_end
  USE lsda_mod,  ONLY : nspin
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE scatter_mod, ONLY : gather_grid
  USE ener,      ONLY: Efermi => ef
  USE gvecw,     ONLY : ecutwfc
  USE pw_restart_new, ONLY : read_collected_wfc
  USE cell_base,  ONLY : at, alat, tpiba, bg



  implicit none
!QE
  CHARACTER (len=256) :: outdir
  CHARACTER(LEN=256), external :: trimcheck
  character(len=256) :: filename
  LOGICAL            :: needwf= .TRUE., exst
  INTEGER :: lgk, mgk, first_k, last_k, first_band, last_band
! real (q) :: xx (3)

  NAMELIST / inputpp / outdir, prefix, first_k, last_k, first_band, last_band




! test if input.STMpw exists, if not stop
     inquire (file = 'input.STMpw', exist = fichero)
     if (fichero) then
! reading input.STMpw
     open (8, file = 'input.STMpw', status='old')
     read (8,*, IOSTAT=ios, ERR=100,END=200) phi
     phi = phi / hartree ! all calculations in atomic units
     read (8,*, IOSTAT=ios, ERR=100,END=200) nV
     allocate (V(nV))
     read (8,*, IOSTAT=ios, ERR=100,END=200) (V(jv),jv=1,nV)
     write(*,'(A,30F7.3)') "Voltages to calculate (in V): ",(V(jv),jv=1,nV)
     do jv=1,nV
      V(jv) = V(jv) / hartree ! all calculations in atomic units
     enddo 
     read (8,*, IOSTAT=ios, ERR=100,END=200) N_sampling_z
     read (8,*, IOSTAT=ios, ERR=100,END=200) Zmax
     Zmax = Zmax / bohr
     read (8,*, IOSTAT=ios, ERR=100,END=200) Zsurf
     read (8,*, IOSTAT=ios, ERR=100,END=200) z_s
     read (8,*, IOSTAT=ios, ERR=100,END=200) Bardeen
     if (Bardeen) then
       write(*,'(A)') "Bardeen calculation. NOT finished for the QE version. Pls do not run."
       read (8,*, IOSTAT=ios, ERR=100,END=200) Ztip
     else
       write(*,'(A)') "Tersoff-Hamman calculation"
       Ztip = 15._q
     endif  
!! matching distance for the tip is hardwired
!! z_t = Ztip - 2.0 Angstroems

! dIdV curve
     read (8,*, IOSTAT=ios, ERR=100,END=200) LDIDV !.true. or .false.
     if (LDIDV) then
       read (8,*, IOSTAT=ios, ERR=100,END=200) Vmin, Vmax
       read (8,*, IOSTAT=ios, ERR=100,END=200) ndiv
       read (8,*, IOSTAT=ios, ERR=100,END=200) npts
       allocate(cdIdV(npts,3))
       do i=1,npts
        read (8,*, IOSTAT=ios, ERR=100,END=200) (cdIdV(i,j),j=1,3)
       enddo 
       if(npts.gt.1000) then
         write(*,'(A)') "We can not calculate more than 1000 dIdV curves."
         npts = 1000
       endif
       cdIdV = cdIdV / bohr
       allocate(V1(ndiv+1))
       allocate(IV(npts,ndiv+1))
       allocate(ngp(npts,3))
       write(*,'(A,F7.3,A,F7.3,A)') "Calculation of IV curves between ",Vmin," eV and ",Vmax," eV"
     endif
     read (8,*, IOSTAT=ios, ERR=100,END=200) NamePOSCAR
     read (8,*, IOSTAT=ios, ERR=100,END=200) NameWF
     read (8,*, IOSTAT=ios, ERR=100,END=200) MAPfile !.true. or .false.
     if (MAPfile) then
      write(*,'(A)') "No MappingsCAR in the Quantum Espresso version." 
      write(*,'(A)') "Please, put the entry to .false. and remove next line," 
      write(*,'(A)') "in the input.STMpw file. We stop." 
      stop
     endif
     read (8,*, IOSTAT=ios, ERR=100,END=200)
     read (8,*, IOSTAT=ios, ERR=100,END=200) LGAMMA !.true. or .false.
     read (8,*, IOSTAT=ios, ERR=100,END=200) wsxm   !.true. or .false.
     if(wsxm) read (8,*, IOSTAT=ios, ERR=100,END=200) factor
     read (8,*, IOSTAT=ios, ERR=100,END=200) dat    !.true. or .false.
     read (8,*, IOSTAT=ios, ERR=100,END=200) cube   !.true. or .false.
 
     close (8)
     else
       write(*,*) ""
       write(*,*) "------------------------------------"
       write(*,*) " Please, provide an input.STMpw file "
       write(*,*) "------------------------------------"
       write(*,*) ""
       stop
     endif

! end of reading input.STMpw

     if(.not.wsxm.AND..not.dat.AND..not.cube.AND..not.LDIDV) then
       write(*,*) ""
       write(*,*) "------------------------------------"
       write(*,'(A)') "No output format has been requested."
       write(*,'(A)') "Please, write an output format. Stopping here."
       write(*,*) "------------------------------------"
       write(*,*) ""
       stop
     endif  


            call disk_units  !set the numbers of the units to write output
                             !please read sbroutine below to see assignments

! We call QE routines to initialize all wf quantities
! following the WFCK2R.f90 code
! this means we need an input pp file:
! &inputpp
!   prefix='MgB2',
!   outdir='./tmp',
! /
! with "prefix" and "outdir" as in the scf/nscf/band calculation.

  CALL mp_startup ( )

  CALL environment_start ( 'WFCK2R' )
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'

  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  IF ( ionode )  THEN
     !
     ! set defaults:
     first_k = 0
     last_k = 0
     first_band = 0
     last_band = 0
     !
     CALL input_from_file ( )
     ! 
     READ (5, inputpp, err = 900, iostat = ios)
900  CALL errore ('WFCK2R', 'reading inputpp namelist', ABS (ios) )
     ! 
     tmp_dir = trimcheck (outdir)
     ! 
  END IF
 !
  ! ... Broadcast variables
  !

  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )
  CALL mp_bcast( first_k, ionode_id, intra_image_comm )
  CALL mp_bcast( last_k, ionode_id, intra_image_comm )
  CALL mp_bcast( first_band, ionode_id, intra_image_comm )
  CALL mp_bcast( last_band, ionode_id, intra_image_comm )


  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file_new ( needwf )

  if (first_k <= 0) first_k = 1
  if (last_k <= 0) last_k = nks
  if (first_band <= 0) first_band = 1
  if (last_band <= 0) last_band = nbnd

    Total_Bands = nbnd
    Number_of_KPOINTS = nks
    Number_of_Spin = nspin

    ngx = dfftp%nr1x
    ngy = dfftp%nr2x
    ngz = dfftp%nr3x

    print *, 'DIMENSIONS for this run::'
    print *, 'Number of bands and Kpoints and spins', nbnd, nks, nspin
    print *, 'ngx, ngy, ngz=', ngx, ngy, ngz
    print *, 'Efermi (eV)=', Efermi*13.605
  


! transformation into atomic units
! of the unit cell vectors
           A = alat * at 
           cell = at ! for compatibility reasons
           lat_par = alat * bohr ! Lattice parameter in Angstroems
           bg = bg *tpiba !reciprocal unit cell

           print *, 'CELL from QE::'
           print *, '         A (:,1)', A (:,1)
           print *, '         A (:,2)', A (:,2)
           print *, '         A (:,3)', A (:,3)
           print *, bg (:,1), sum (A(:,1)*bg(:,1))/(2*pi_d)
           print *, bg (:,2), sum (A(:,2)*bg(:,2))/(2*pi_d)
           print *, bg (:,3), sum (A(:,3)*bg(:,3))/(2*pi_d)

! third dimension in 90 degrees with respect to surface
        if (A (1,3) > 0.1 .or. A(2,3) > 0.1 ) then

              print *, 'The calculation has to be repeated'
              print *, 'with the 3rd axis perpendicular to the surface plane!'
              stop

        endif
      call volume (A(:,1),A(:,2),A(:,3),vol)
      call surface (A(:,1),A(:,2),surf)

!NL all quantities must be defined before this point!!!
! general allocation
       call allocation

! spin factor for QE
      if (nspin == 1 .or. nspin ==4) then
             spin_factor = 1
      else
             spin_factor = 2
      endif



      Zmin = (z_s - Zsurf) * A(3,3)

      stepX = sqrt(sum(A(:,1)**2))/(Ngx)
      stepY = sqrt(sum(A(:,2)**2))/(Ngy)
! matching point to tip distance
!     distance = (Ztip-z_t) * A(3,3)
      distance = 2.0/bohr !Hardwired 2.0 Ang
      z_t=Ztip-distance/A(3,3)

! z is the difference
! between matching points in this program
! in other words: when z=0 the distance
! between tip and sample is the offset
      offset = Zmin + distance
      stepZ = Zmax / (N_sampling_z-1)
! Hence the initial distance is "offset"
! the maximum distance is "Zmax+offset"

! Initialize current and conductance

       Tersoff_c =  0
       Tersoff_s =  0

      if (LDIDV) then
          call setting_dIdV
          intensity2 = 0
          Tersoff_s2 = 0
      endif  

     if(Bardeen) then
       intensity = 0
       dintensity_dV = 0
       Tersoff_t =  0
     endif

! reciprocal vectors
! Predefining FFT ordering for gx and gy

     i=0
     do i1 = 0, ngx/2
     i = i+1
     ga (i) = i1
     enddo
     do i1 = 1, ngx/2-1
     i = i+1
     ga (i) = (-ngx/2+i1)
     enddo
     i=0
     do i2 = 0, ngy/2
     i = i+1
     gb (i) = i2
     enddo
     do i2 = 1, ngy/2-1
     i = i+1
     gb (i) = (-ngy/2+i2)
     enddo


! QE first k-point then spin
! Loop on spin and k-points, but actually
! spin is included in the k-point loop

! We start reading the WF file

   KPOINTS: do kp=1, Number_of_KPOINTS

          write(*,'(A,i4,A)') "  Calculating kpoint ", kp, " ..."

        NPL=ngk(kp)

! reciprocal vectors

      gx = 0._q; gy = 0._q ! padding all values to zero

      i=0
      do iy = 1,ngy
        do ix = 1,ngx

      i=i+1

          g1 = ga(ix) + xk(1,kp)
          g2 = gb(iy) + xk(2,kp)

          gx (i) = ( g1*bg(1,1) + g2*bg(2,1) ) 
          gy (i) = ( g1*bg(1,2) + g2*bg(2,2) ) 

        enddo
      enddo


!       xx (1)=xk(1,kp); xx (2)=xk(2,kp);xx (3)=xk(3,kp);
      print *, 'gx (1:3)=', gx (1:3)

     write(*,*) "   Reading and matching wavefunctions ..."

     t_wf = 0.0_q
     t_ma = 0.0_q

     call cpu_time(start)

!read in bands
     CALL read_collected_wfc ( restart_dir(), kp, evc )

         EIG = zero
     BANDS: do iband = 1, Total_Bands

! EIG in Hartree from the Fermi energy
         EIG (iband) = cmplx (et (iband,kp) /2 - Efermi/2,0._q)

       psic (:) =(0.0,0.0)
       CW (:) =(0.0,0.0)
       call cpu_time(st_tmp)
       do l=1,NPL
       psic (dfftp%nl(igk_k(l,kp))) = evc (l, iband)
       enddo
! Transform as in QE and then transform back using fftw
         call invfft ('Rho', psic, dfftp) ! real space
         call fft3d (psic, ngx, ngy, ngz) ! Now it is in usual fftw format
              CW = psic/(ngx*ngy*ngz)
  
       call cpu_time(fi_tmp)
       t_wf = t_wf + fi_tmp - st_tmp


! obtain 2-D coefficients for the surface and the tip
! by performing the exponential matching

   call cpu_time(st_tmp)
!  call matching (CW, A_G, C_G, ngx, ngy, ngz, z_s, z_t,bg,xx,NPL,ecutwfc,temp)
   call matching (CW, A_S, C_T, ngx, ngy, ngz, z_s,z_t,temp)


   call cpu_time(fi_tmp)
   t_ma = t_ma + fi_tmp - st_tmp

     A_G(:,:,iband)=A_S(:,:)
     C_G(:,:,iband)=C_T(:,:)
! close Loop on bands

           enddo BANDS
     call cpu_time(finish)
     write (*,'(A, F8.1, A)') "     Time to read and match wavefunctions = ", finish-start, " s"
     write (*,'(A, F8.1, A)') "       Time to read wavefunctions = ", t_wf, " s"
     write (*,'(A, F8.1, A)') "       Time to match wavefunctions = ", t_ma, " s"

! Compute the current calculation for all different tip
! positions in the cell and from Z_min to Z_max.
! Loop on tip-surface distances:

       call cpu_time(start)
       write(*,*) "   Computing current ..."

! Loop on voltages
   
   VOLTAGES: do jv = 1, nV
      
       write(*,'(A,F7.3)') "      Computing V = ", V(jv)*hartree

       Heights: do iz=1,N_sampling_z

       z=stepZ*(iz-1)


       if (Bardeen) then
! Loop on tip bands
         do tband=1,Total_Bands


! Loop on substrate bands
         do sband=1,Total_Bands
  
! Check voltage sign (tip to mass, then positive
! means empty substrate states)
   
    if ( V(jv) < 0 ) then

          Fermi_t = 1-Fermi(dreal(EIG(tband)))
          Fermi_s = Fermi (dreal(EIG(sband)))

     constant_I= - 2 * pi_d /( sqrpi * sigma)
     constant_dIdV=  4 * pi_d /( sqrpi * sigma**3)

    else

          Fermi_t = Fermi(dreal(EIG(tband)))
          Fermi_s = 1-Fermi (dreal(EIG(sband)))

     constant_I=  2 * pi_d /( sqrpi * sigma)
     constant_dIdV= - 4 * pi_d /( sqrpi * sigma**3)

    endif
     constant_I=spin_factor*constant_I*surf**2/vol**2
     constant_dIdV=spin_factor*constant_dIdV*surf**2/vol**2

! before calculating we check for empty states ...
! and energy conservation

             if (Fermi_t*Fermi_s > 0.1) then

             delta_E = ((dreal(EIG(tband)-EIG(sband)) + V(jv))/sigma) **2

             if (delta_E <3.) then

! compute Bardeen's matrix element: currentSQ which is already
! the squared modulus

     A_S(:,:)=A_G(:,:,sband)
     C_T(:,:)=C_G(:,:,tband)
 
   call matrix_element (A_S,C_T,ngx,ngy,phi,gx,gy,z,currentSQ)

! compute current

      if (.not.LGAMMA) then
   intensity(:,:,iz,jv) = intensity (:,:,iz,jv) + constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E)/                    &
    & (Number_of_KPOINTS)
      elseif (kp == 1) then
   intensity(:,:,iz,jv) = intensity (:,:,iz,jv) + constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E)/                    &
    & (2*Number_of_KPOINTS-1)
      else
   intensity(:,:,iz,jv) = intensity (:,:,iz,jv) + 2.0*constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E)/                    &
    & (2*Number_of_KPOINTS-1)
      end if



! compute conductance

      if (.not.LGAMMA) then
dintensity_dV(:,:,iz,jv) = dintensity_dV(:,:,iz,jv) + constant_dIdV * Fermi_t * &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E) *  &
    & (dreal(EIG(tband)-EIG(sband)) + V(jv))/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
dintensity_dV(:,:,iz,jv) = dintensity_dV(:,:,iz,jv) + constant_dIdV * Fermi_t * &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E) *  &
    & (dreal(EIG(tband)-EIG(sband)) + V(jv))/ &
       & (2*Number_of_KPOINTS-1)
      else
dintensity_dV(:,:,iz,jv) = dintensity_dV(:,:,iz,jv) + 2.0*constant_dIdV * Fermi_t * &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E) *  &
    & (dreal(EIG(tband)-EIG(sband)) + V(jv))/ &
       & (2*Number_of_KPOINTS-1)
      end if

             endif
             endif

          enddo
          enddo

      endif ! Bardeen part

! Redo calculation to
! calculate Tersoff-Hamman picture

    if ( V(jv) < 0 ) then
     constant_I= - 2 * pi_d /( sqrpi * sigma)
    else
     constant_I=  2 * pi_d /( sqrpi * sigma)
    endif
     constant_I=spin_factor*constant_I*surf**2/vol**2
     
! Loop on bands
         do tband=1,Total_Bands

         if ( V(jv) > 0) then
         if (dreal(EIG(tband)) > 0 .and. dreal(EIG(tband)) < V(jv) ) then 
      
     A_S(:,:)=A_G(:,:,tband)
     C_T(:,:)=C_G(:,:,tband)

   call matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z,TH_t,TH_s)

      if (.not.LGAMMA) then
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (Number_of_KPOINTS) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I 
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + TH_t (:,:)/vol/ &
       & (Number_of_KPOINTS)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I 
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + TH_t (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
      else
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + 2.0*TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I 
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + 2.0*TH_t (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + 2.0*TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
      end if

        endif

        else

        if (dreal(EIG(tband)) < 0 .and. dreal(EIG(tband)) > V(jv) ) then 
      
     A_S(:,:)=A_G(:,:,tband)
     C_T(:,:)=C_G(:,:,tband)

   call matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z,TH_t,TH_s)

      if (.not.LGAMMA) then
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (Number_of_KPOINTS) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + TH_t (:,:)/vol/ &
       & (Number_of_KPOINTS)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + TH_t (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
      else
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + 2.0*TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + 2.0*TH_t (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + 2.0*TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
      end if

          endif
          endif

! Loop on bands
         enddo

 
           enddo Heights


         enddo VOLTAGES  
     call cpu_time(finish)
     write (*,'(A, F8.1, A)') "    Time for current = ", finish-start, " s"

! Calculation on dIdV curves

    if(LDIDV) then 
       call cpu_time(start)
       write(*,*) "   Computing dIdV curves ..."

  npts_dIdV: do ip = 1,npts
   dIdV: do jv = 1, ndiv
      
       iz=ngp(ip,3)

       z=stepZ*(iz-1)


       if (Bardeen) then
! Loop on tip bands
         do tband=1,Total_Bands


! Loop on substrate bands
         do sband=1,Total_Bands
  
! Check voltage sign (tip to mass, then positive
! means empty substrate states)
   
    if ( V1(jv) < 0 ) then

          Fermi_t = 1-Fermi(dreal(EIG(tband)))
          Fermi_s = Fermi (dreal(EIG(sband)))

     constant_I= - 2 * pi_d /( sqrpi * sigma)
     constant_dIdV=  4 * pi_d /( sqrpi * sigma**3)

    else

          Fermi_t = Fermi(dreal(EIG(tband)))
          Fermi_s = 1-Fermi (dreal(EIG(sband)))

     constant_I=  2 * pi_d /( sqrpi * sigma)
     constant_dIdV= - 4 * pi_d /( sqrpi * sigma**3)

    endif
     constant_I=spin_factor*constant_I*surf**2/vol**2
     constant_dIdV=spin_factor*constant_dIdV*surf**2/vol**2

! before calculating we check for empty states ...
! and energy conservation

             if (Fermi_t*Fermi_s > 0.1) then

             delta_E = ((dreal(EIG(tband)-EIG(sband)) + V1(jv))/sigma) **2

             if (delta_E <3.) then

! compute Bardeen's matrix element: currentSQ which is already
! the squared modulus

     A_S(:,:)=A_G(:,:,sband)
     C_T(:,:)=C_G(:,:,tband)
 
   call matrix_element (A_S,C_T,ngx,ngy,phi,gx,gy,z,currentSQ)

! compute current

      if (.not.LGAMMA) then
   intensity2(ip,jv) = intensity2 (ip,jv) + constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(ngp(ip,1),ngp(ip,2)) * exp (-delta_E)/                    &
    & (Number_of_KPOINTS)
      elseif (kp == 1) then
   intensity2(ip,jv) = intensity2 (ip,jv) + constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(ngp(ip,1),ngp(ip,2)) * exp (-delta_E)/                    &
    & (2*Number_of_KPOINTS-1)
      else
   intensity2(ip,jv) = intensity2 (ip,jv) + 2.0*constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(ngp(ip,1),ngp(ip,2)) * exp (-delta_E)/                    &
    & (2*Number_of_KPOINTS-1)
      end if

             endif
             endif

          enddo
          enddo

      endif ! Bardeen part

! Redo calculation to
! calculate Tersoff-Hamman picture

    if ( V1(jv) < 0 ) then
     constant_I= - 2 * pi_d /( sqrpi * sigma)
    else
     constant_I=  2 * pi_d /( sqrpi * sigma)
    endif
     constant_I=spin_factor*constant_I*surf**2/vol**2
     
! Loop on bands
         do tband=1,Total_Bands

         if ( V1(jv) > 0) then
         if (dreal(EIG(tband)) > 0 .and. dreal(EIG(tband)) < V1(jv) ) then 
      
     A_S(:,:)=A_G(:,:,tband)
     C_T(:,:)=C_G(:,:,tband)

   call matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z,TH_t,TH_s)

      if (.not.LGAMMA) then
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (2*Number_of_KPOINTS-1)
      else
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + 2.0*TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (2*Number_of_KPOINTS-1)
      end if


         endif

         else

         if (dreal(EIG(tband)) < 0 .and. dreal(EIG(tband)) > V1(jv) ) then 
      
     A_S(:,:)=A_G(:,:,tband)
     C_T(:,:)=C_G(:,:,tband)

   call matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z,TH_t,TH_s)

      if (.not.LGAMMA) then
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (2*Number_of_KPOINTS-1)
      else
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + 2.0*TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (2*Number_of_KPOINTS-1)
      end if

          endif
          endif

! Loop on bands
         enddo

         enddo dIdV      
        enddo npts_dIdV
     call cpu_time(finish)
     write (*,'(A, F8.1, A)') "    Time for computing dIdV curves = ", finish-start, " s"
    endif 


         enddo KPOINTS

! output in WSxM format
! *.siesta
! output in ASCII to be converted into STM.dato
! and read by imagen.f into Mathematica, gnuplot, ...

        write(*,*) "Writing files ..."
        call cpu_time(start)
        if(wsxm) write(*,'(A,F8.0)') "Output in WSxM format. Output will be multiplied by ",factor
        if(dat) write(*,'(A)') "Output for gnuplot"
        if(cube) write(*,'(A)') "Output for TH in cube format"
      do jv=1,nV 
        write(volt,'(F6.3)') V(jv)*hartree
        cmd = 'mkdir V_'//trim(adjustl(volt))
        call system(cmd) 
        UCELL = A
        UCELL (3,3) = Zmax

        if(wsxm) then
         name_file = 'V_'//trim(adjustl(volt))//'/Bardeen_V_'//trim(adjustl(volt))//'.siesta'
         if (Bardeen) open (unitI,file=name_file, form = 'unformatted')
         name_file = 'V_'//trim(adjustl(volt))//'/TH_V_'//trim(adjustl(volt))//'.siesta'
         open (unitTH,file=name_file, form = 'unformatted')
         name_file = 'V_'//trim(adjustl(volt))//'/TH_tip_V_'//trim(adjustl(volt))//'.siesta'
         if (Bardeen) open (unitTIP,file=name_file, form = 'unformatted')
         if (Bardeen) write (unitI) UCELL
         write (unitTH) UCELL
         if (Bardeen) write (unitTIP) UCELL
         if (Bardeen) write (unitI) ngx, ngy, N_sampling_z,1
         write (unitTH) ngx, ngy, N_sampling_z,1
         if (Bardeen) write (unitTIP) ngx, ngy, N_sampling_z,1
         do iz= 1, N_sampling_z
          do iy = 0, ngy-1
           if (Bardeen) write (unitI) (real(abs(intensity (ix,iy,iz,jv))*factor), ix=0,ngx-1)
           write (unitTH) (real(Tersoff_s(ix,iy,iz,jv)*factor), ix=0,ngx-1)
           if (Bardeen) write (unitTIP) (real(Tersoff_t(ix,iy,iz,jv)*factor), ix=0,ngx-1)
          enddo
         enddo
         if (Bardeen) close (unitI)
         close (unitTH)
         if (Bardeen) close (unitTIP)

        endif

        if(dat) then
          name_file = 'V_'//trim(adjustl(volt))//'/Bardeen.dat'
          if (Bardeen) open (unitIdat,file=name_file)
          name_file = 'V_'//trim(adjustl(volt))//'/TH.dat'
          open (unitTHdat,file=name_file)
          name_file = 'V_'//trim(adjustl(volt))//'/TH_tip.dat'
          if (Bardeen) open (unitTIPdat,file=name_file)
          name_file = 'V_'//trim(adjustl(volt))//'/dIdV_TH.dat'
          open (unitdITHdat,file=name_file)
          name_file = 'V_'//trim(adjustl(volt))//'/dIdV_Bardeen.dat'
          if (Bardeen) open (unitdIdat,file=name_file)
          write(unitdITHdat,*) N_sampling_z, ngy, ngx
          write(unitTHdat,*) N_sampling_z, ngy, ngx
          if (Bardeen) write(unitTIPdat,*) N_sampling_z, ngy, ngx
          if (Bardeen) write(unitIdat,*) N_sampling_z, ngy, ngx
          if (Bardeen) write(unitdIdat,*) N_sampling_z, ngy, ngx
          do iz= 1, N_sampling_z
           do iy = 0, ngy-1
            do ix = 0, ngx-1
             write (unitdITHdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)),bohr*(iy)*stepY, bohr*(ix)*stepX, Tersoff_c(ix,iy,iz,jv)
             write (unitTHdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)),bohr*(iy)*stepY, bohr*(ix)*stepX, Tersoff_s(ix,iy,iz,jv)
             if (Bardeen) write (unitTIPdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)), bohr*(iy)*stepY, bohr*(ix)*stepX, Tersoff_t(ix,iy,iz,jv) 
             if (Bardeen) write (unitIdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)),bohr*(iy)*stepY, bohr*(ix)*stepX, intensity (ix,iy,iz,jv)
             if (Bardeen) write (unitdIdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)),bohr*(iy)*stepY, bohr*(ix)*stepX, dintensity_dV (ix,iy,iz,jv)
            enddo
           enddo
          enddo

          close(unitTHdat)
          if (Bardeen) close(unitTIPdat)
          if (Bardeen) close(unitIdat)
          if (Bardeen) close(unitdIdat)
          close(unitdITHdat)
        endif

        if (cube) then
          name_file = 'V_'//trim(adjustl(volt))//'/TH_V_'//trim(adjustl(volt))//'.cube'
          open (unitTHcub,file=name_file)
          write(unitTHcub,'(A)') "STM image in TH approximation."
          write(unitTHcub,'(3A,F6.3,A)') "V = ",volt," V. z = 0 is at",Zmin*bohr," angs over the surface."
          write(unitTHcub,'(i5,3F12.6)') nat, 0.0, 0.0, z_s*cell(3,3)*lat_par/bohr
          write(unitTHcub,'(i5,3F12.6)') ngx, (UCELL(j,1)/ngx,j=1,3)
          write(unitTHcub,'(i5,3F12.6)') ngy, (UCELL(j,2)/ngy,j=1,3)
          write(unitTHcub,'(i5,3F12.6)') N_sampling_z, (UCELL(j,3)/N_sampling_z,j=1,3)
          do i=1,nat
            write(unitTHcub,'(i5,4F12.6)') zat(i),0.000,(coord(i,j)/bohr,j=1,3)
          enddo  
          do ix = 0, ngx-1
            do iy = 0, ngy-1
              write(unitTHcub, '(6E13.5)') (Tersoff_s(ix,iy,iz,jv), iz=1,N_sampling_z)
            enddo
          enddo  
          close(unitTHcub)
! dIdV          
          name_file = 'V_'//trim(adjustl(volt))//'/dIdV_TH_V_'//trim(adjustl(volt))//'.cube'
          open (unitdITHcub,file=name_file)
          write(unitdITHcub,'(A)') "dIdV in TH approximation."
          write(unitdITHcub,'(3A,F6.3,A)') "V = ",volt," V. z = 0 is at",Zmin*bohr," angs over the surface."
          write(unitdITHcub,'(i5,3F12.6)') nat, 0.0, 0.0, z_s*cell(3,3)*lat_par/bohr
          write(unitdITHcub,'(i5,3F12.6)') ngx, (UCELL(j,1)/ngx,j=1,3)
          write(unitdITHcub,'(i5,3F12.6)') ngy, (UCELL(j,2)/ngy,j=1,3)
          write(unitdITHcub,'(i5,3F12.6)') N_sampling_z, (UCELL(j,3)/N_sampling_z,j=1,3)
          do i=1,nat
            write(unitdITHcub,'(i5,4F12.6)') zat(i),0.000,(coord(i,j)/bohr,j=1,3)
          enddo  
          do ix = 0, ngx-1
            do iy = 0, ngy-1
              write(unitdITHcub, '(6E13.5)') (abs(Tersoff_c(ix,iy,iz,jv)), iz=1,N_sampling_z)
            enddo
          enddo  
          close(unitdITHcub)
        endif 

     enddo ! Voltages

     if (LDIDV) then
      cmd = 'mkdir dIdV_curves'
      call system(cmd)
      do ip=1,npts
        unitIV=1000+ip
        write(cnpt,'(I6)') ip
        name_file = 'dIdV_curves/IV_TH_npt_'//trim(adjustl(cnpt))//'.dat'
        open (unitIV,file=name_file)
!        write(unitIV,'(A,3F7.3,A)') "# IV for point (",(cdIdV(ip,j)*bohr,j=1,3),")"
        write(unitIV,'(A,3F7.3,A)') "# Actual point (", &
                cngx(ngp(ip,1),ngp(ip,2)),cngy(ngp(ip,1),ngp(ip,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(ip,3)-1)),")"
        do jv=1,ndiv-1
         if(V1(jv).lt.0) Tersoff_s2(ip,jv) = -1.0 * Tersoff_s2(ip,jv)
         write(unitIV,*) V1(jv)*hartree,Tersoff_s2(ip,jv)
        enddo
        close (unitIV)
      enddo ! points for dIdV
      do ip=1,npts
        unitIV=2000+ip
        write(cnpt,'(I6)') ip
        name_file = 'dIdV_curves/dIdV_TH_npt_'//trim(adjustl(cnpt))//'.dat'
        open (unitIV,file=name_file)
!        write(unitIV,'(A,3F7.3,A)') "# dIdV for point (",(cdIdV(ip,j)*bohr,j=1,3),")"
        write(unitIV,'(A,3F7.3,A)') "#   Actual point (", &
                cngx(ngp(ip,1),ngp(ip,2)),cngy(ngp(ip,1),ngp(ip,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(ip,3)-1)),")"
        do jv=1,ndiv-1
         write(unitIV,*) V1(jv)*hartree,(Tersoff_s2(ip,jv+1)-Tersoff_s2(ip,jv))/(V1(jv+1)-V1(jv))/hartree
        enddo
        close (unitIV)
      enddo ! points for dIdV
!
     if (Bardeen) then
      do ip=1,npts
        unitIV=3000+ip
        write(cnpt,'(I6)') ip
        name_file = 'dIdV_curves/IV_Bardeen_npt_'//trim(adjustl(cnpt))//'.dat'
        open (unitIV,file=name_file)
!        write(unitIV,'(A,3F7.3,A)') "# IV for point (",(cdIdV(ip,j)*bohr,j=1,3),")"
        write(unitIV,'(A,3F7.3,A)') "# Actual point (", &
                cngx(ngp(ip,1),ngp(ip,2)),cngy(ngp(ip,1),ngp(ip,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(ip,3)-1)),")"
        do jv=1,ndiv-1
         if(V1(jv).lt.0) intensity2(ip,jv) = -1.0 * intensity2(ip,jv)
         write(unitIV,*) V1(jv)*hartree,intensity2(ip,jv)
        enddo
        close (unitIV)
      enddo ! points for dIdV
      do ip=1,npts
        unitIV=4000+ip
        write(cnpt,'(I6)') ip
        name_file = 'dIdV_curves/dIdV_Bardeen_npt_'//trim(adjustl(cnpt))//'.dat'
        open (unitIV,file=name_file)
!        write(unitIV,'(A,3F7.3,A)') "# dIdV for point (",(cdIdV(ip,j)*bohr,j=1,3),")"
        write(unitIV,'(A,3F7.3,A)') "#   Actual point (", &
                cngx(ngp(ip,1),ngp(ip,2)),cngy(ngp(ip,1),ngp(ip,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(ip,3)-1)),")"
        do jv=1,ndiv-1
         write(unitIV,*) V1(jv)*hartree,(intensity2(ip,jv+1)-intensity2(ip,jv))/(V1(jv+1)-V1(jv))/hartree
        enddo
        close (unitIV)
      enddo ! points for dIdV
     endif ! Bardeen
    endif ! dIdV curves

     call cpu_time(finish)
     write (*,'(A, F8.1, A)') " Time to write files = ", finish-start, " s"

100  if (ios>0) then
     print *,'ERROR reading input.STMpw (wrong format or corrupt file), IOS=',ios
     endif
200  if (ios<0) then
     print *,'End of input.STMpw (empty file or not enough data), IOS=',ios
     endif

CONTAINS

! allocation of memory

       subroutine allocation

       allocate (EIG(Total_bands))
       allocate (A_G(0:ngx-1,0:ngy-1,Total_bands))
       allocate (C_G(0:ngx-1,0:ngy-1,Total_bands))
       allocate (A_S(0:ngx-1,0:ngy-1))
       allocate (C_T(0:ngx-1,0:ngy-1))
       allocate (cngx(0:ngx,0:ngy))
       allocate (cngy(0:ngx,0:ngy))
       allocate (CW(ngx*ngy*ngz))
       allocate (density_z(ngz))
       allocate (density1(ngz))
       allocate (temp(0:ngx-1,0:ngy-1,0:ngz-1))

       if(Bardeen) then
         allocate (intensity(0:ngx-1,0:ngy-1,N_sampling_z,nV))
         allocate (Tersoff_t(0:ngx-1,0:ngy-1,N_sampling_z,nV))
         allocate (dintensity_dV(0:ngx-1,0:ngy-1,N_sampling_z,nV))
         if(LDIDV) allocate (intensity2(npts,ndiv))
       endif  
       allocate (Tersoff_s(0:ngx-1,0:ngy-1,N_sampling_z,nV))
       allocate (Tersoff_c(0:ngx-1,0:ngy-1,N_sampling_z,nV))
       if(LDIDV) allocate (Tersoff_s2(npts,ndiv))
       allocate (currentSQ(0:ngx-1,0:ngy-1))
       allocate (TH_s(0:ngx-1,0:ngy-1))
       allocate (TH_t(0:ngx-1,0:ngy-1))
       allocate (ga(NGX))
       allocate (gb(NGY))
       allocate (gx(NGX*NGY))
       allocate (gy(NGX*NGY))

       return
       end subroutine allocation

   subroutine disk_units
! Reciprocal vector and wave index file
      unitMAP = 8
! output current units
      unitI = 9
      unitIdat = 19
! output conductance unit
      unitdI = 10
! output conductance unit
      unitdIdat = 25
! output current units Tersoff-Hamman
      unitTH = 11
      unitTHdat = 21
! output conductance unit TH
      unitdITHdat = 20
! output current units Tersoff-Hamman for the TIP
      unitTIP =12
      unitTIPdat =22
! WF input file unit
      unitWF = 13
! output current unit Tersoff-Hamman in cube format
      unitTHcub = 14
! output current unit Tersoff-Hamman in cube format
      unitdITHcub = 15
     return
     end subroutine disk_units


     subroutine setting_dIdV
        do jv=1,ndiv+1
          V1(jv) = Vmin + (Vmax-Vmin)/(ndiv-1)*(jv-1) 
          V1(jv) = V1(jv) / hartree
!          write(*,*) "V(",jv,")=",V1(jv)
        enddo
        ndiv = ndiv + 1

        do iy = 0, ngy
          do ix = 0, ngx
            cngx(ix,iy)=cell(1,1)*ix/(ngx+1)+cell(2,1)*iy/(ngy+1)
            cngy(ix,iy)=cell(1,2)*ix/(ngx+1)+cell(2,2)*iy/(ngy+1)
            cngx(ix,iy)=cngx(ix,iy)*lat_par
            cngy(ix,iy)=cngy(ix,iy)*lat_par
          enddo
        enddo
!
        write(*,'(A)') "Points to calculate dIdV curves (in angs):"
        do i=1,npts
         ngp(i,1)=ngx+100
         do iy = 0, ngy
         do ix = 0, ngx
          if(abs(cdIdV(i,1)-cngx(ix,iy)/bohr).lt.stepX/2.AND.abs(cdIdV(i,2)-cngy(ix,iy)/bohr).lt.stepY/2) then
           ngp(i,1)=ix
           ngp(i,2)=iy
          endif
         enddo
         enddo
          if(ngp(i,1).eq.ngx+100) then
           write(*,'(A,i4,A)') "We do not calculate de dIdV curves because point",i," is out of range:"
           write(*,'(A,2F8.3,A)') "(",(cdIdV(i,j)*bohr,j=1,2),") not in the unit cell."
           LDIDV = .false.
          endif
         ngp(i,3)=N_sampling_z+100
         do iz=1,N_sampling_z
          if(abs(cdIdV(i,3)-(z_s*A(3,3)+stepZ*(iz-1))).lt.stepZ/2) then
           ngp(i,3)=iz
          endif
         enddo
          if(ngp(i,3).eq.N_sampling_z+100) then
           write(*,'(A,i4,A)') "We do not calculate de dIdV curves because point",i," is out of range:"
           write(*,'(F8.3,A,F8.3,A,F8.3,A)') & 
                   cdIdV(i,3)*bohr," not in [", z_s*A(3,3)*bohr,",",(z_s*A(3,3)+Zmax)*bohr,"]"
           LDIDV = .false.
          endif
        write(*,'(A,i3,A,3F8.3,A)') " Requested point",i," = (",(cdIdV(i,j)*bohr,j=1,3),")"
        write(*,'(A,i3,A,3F8.3,A)') "    Actual point",i," = (",&
           cngx(ngp(i,1),ngp(i,2)),cngy(ngp(i,1),ngp(i,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(i,3)-1)),")"
        enddo
        return
     end subroutine setting_dIdV
end program STMpw
