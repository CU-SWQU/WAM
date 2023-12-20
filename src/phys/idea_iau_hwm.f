      module idea_iau_hwm
      USE machine,       ONLY: kind_io8     
      implicit NONE
      integer, parameter  :: nxa_hwm = 768, nya_hwm =384, nza_hwm=40, nta_hwm=24  ! fixed dimension for T254-nudging HWM

      real                :: lat_hwm(nya_hwm)
      real                :: lon_hwm(nxa_hwm)
      real                :: lev_hwm(nza_hwm)

      integer             :: irec_anl, Last_IAU_DATE8
!
      real, dimension(nza_hwm):: rtau
      real, dimension(nza_hwm):: wm_nudge_hwm, wa_nudge_hwm

      real, parameter     ::  pnudge= 200.0  ! 100 Pa or 1mb  48 km ; 200 -43 km
      integer             ::  knudge         ! last v-level of nudging
      integer             ::  dknudge= 10    ! 10  # of layers in the transition zone
      integer             ::  ktrans
      integer             ::  ktrans_hwm
      integer             ::  knudge_hwm


      integer, dimension(6) :: ivar_nudge !/1,1,1,0,0,1/   !T-U-V-Q-Ts-Ps/             

!
      real(kind=kind_io8), dimension(nxa_hwm, nya_hwm, nza_hwm) :: Uhwm, Vhwm


!      character(len=132)  ::  
!     & Dir_list_hwm='/glade/scratch/garimam/hwm_nudge/ncdf_hr'
      character(len=200)  ::  
     & Dir_list_hwm='/glade/scratch/garimam/hwm_nudge/'//
     & 'hwm_at_wam_levs_2010/small_ncdf_hr_24'
      character(len=25)  ::  FILE_PREF_hwm='hr_hwm_24_wam_levs'
      character(len=3)    :: Fend_hwm ='.nc'
      character(len=300)  :: File_anl_hwm   !=trim(Dir_list//File_day//Fend)

      real(kind=kind_io8), allocatable, dimension(:, :, :) :: gpy_U_hwm, gpy_V_hwm

      data ivar_nudge/1,1,1,0,0,1/
      contains

      SUBROUTINE IAU_UPDATE_DATE8(Idate, Fhour, CurDATE8)
      use wam_date_calendar, only : CURRENT_NCEP_JDAT
      IMPLICIT NONE
      integer , intent(in)     :: idate(4)
      integer, intent(inout)   :: CurDATE8
      real , intent(in)        :: Fhour
!locals
      integer                  :: IAU_DATE_8 
      integer                  :: ndays
      integer                  :: y4, m2, d2, d365
!
      integer, parameter       :: ndi=4
      integer, dimension(ndi)  :: jdat, idat
      real                     :: hr_jdat
      IAU_DATE_8 = Idate(4)*10000+Idate(2)*100+Idate(3)
      d365 = 365

      if (mod(Idate(4), 4).eq. 0) d365=366
      ndays = mod(nint(fhour), 24)

      if (ndays <= d365 ) then
          y4=Idate(4)
      else 
          y4=Idate(4)+1
      endif
      idat(1) =idate(4)
      idat(2) =idate(2)
      idat(3) =idate(3)
      idat(4) = 0
      CALL CURRENT_NCEP_JDAT(idat, Fhour,jdat, hr_jdat)
      CurDATE8=Jdat(1)*10000+Jdat(2)*100+Jdat(3)
       
      RETURN
      END SUBROUTINE IAU_UPDATE_DATE8
!
      SUBROUTINE UPDATE_DATE10(Idate, Fhour, CurDATE10)
      use wam_date_calendar, only : CURRENT_NCEP_JDAT
      IMPLICIT NONE
      integer , intent(in)     :: idate(4)
      integer, intent(inout)   :: CurDATE10
      real , intent(in)        :: Fhour
!locals
      integer                  :: IAU_DATE_8 
      integer                  :: ndays
      integer                  :: y4, m2, d2, d365
!
      integer, parameter       :: ndi=4
      integer, dimension(ndi)  :: jdat, idat
      real                     :: hr_jdat
      real                     :: hr_day
      IAU_DATE_8 = Idate(4)*10000+Idate(2)*100+Idate(3)
      d365 = 365

      if (mod(Idate(4), 4).eq. 0) d365=366
      ndays = mod(nint(fhour), 24)
      hr_day = fhour -24.*ndays
      if (ndays <= d365 ) then
          y4=Idate(4)
      else 
          y4=Idate(4)+1
      endif
      idat(1) =idate(4)
      idat(2) =idate(2)
      idat(3) =idate(3)
      idat(4) = 0
      CALL CURRENT_NCEP_JDAT(idat, Fhour,jdat, hr_jdat)
      CurDATE10=Jdat(1)*1000000+Jdat(2)*10000+Jdat(3)*100+nint(hr_jdat)
       
      RETURN
      END SUBROUTINE UPDATE_DATE10
!
      subroutine handle_err( status, message )
      implicit none

      integer,   intent(in)        :: status
      character(len=*), intent(in) :: message
      print*,'** handle_err status = ',status, trim(message)
      stop

      end subroutine handle_err
! 
! TO DO Adapt   IDEA_IAU_SPLIT  intlon_phys(iord,imon,imsk,m1,m2,k1,f1,f2)
!               interpred_phys(1,kmsk0,buffo,buff2, nyg, lonsperlar)
!
        SUBROUTINE IDEA_INIT_IAU_HWM(DATE_8, Irec, me)
    
!       read IAU- coordinates and check dimensions nxa = 192, nya =94, nza=82, nta=4
!                       start = (/ 1, 1, numTimes /),     &
!                           count = (/ numLats, numLons, 1 /))
          use netcdf
          implicit NONE   
          integer, intent(in) ::   DATE_8, Irec, me
          integer    :: nx, ny, nz, nt
          integer    :: status, ncid
          integer    :: iLevDimID, LevDimID, LonDimID, LatDimID
          integer    :: DateDimID, RecordDimID
          integer    :: date_id, var_id
          integer    :: start2(3), count2(3)
          integer    :: start3(4), count3(4)
          integer    :: date_hwm(nta_hwm), date_sec_hwm(nta_hwm)
          integer    :: DATE_4  
          character(len=132)  :: Ftrial
          character(len=8)    :: S8_DATE
          character(len=4)    :: SY_DATE
          integer :: k, kk, kz
          real    :: pwam_pa, tau_6hr, tau_step, decay
          tau_6hr = 6.*3600.   
          !tau_step =180.
          tau_step=1800*2
       
          !print*, 'Gar init func called', DATE_8, Irec, me
          !print *, 'Garr init : ', nta_hwm
          write(S8_DATE, fmt='(I8.8)') DATE_8
          
          !print *, S8_DATE, ' STR-date '

          DATE_4=DATE_8/10000
          !write(SY_DATE, fmt='(I4.4)') DATE_4
          !File_anl_hwm=trim(Dir_list_hwm)//'/'
          File_anl_hwm=trim(Dir_list_hwm)//'/'//trim(FILE_PREF_hwm)//trim(S8_DATE)//trim(Fend_hwm)
          if (me==0) then
              print*, 'Gar init func called', Date_8, Irec, me
              print*, 'Gar File_anl ', File_anl_hwm
          endif
          !print *,  ' Gar File_anl ', File_anl_hwm

          status = nf90_open(File_anl_hwm, nf90_nowrite, ncid)
          if (status /= nf90_noerr) call handle_err(status, ' nf90_open')
          status = nf90_inq_dimid(ncid, 'ntime', RecordDimID)
          if (status /= nf90_noerr) then
            call handle_err(status, 'RecordDimID')
          endif
          status = nf90_inq_dimid(ncid, "nlevs", LevDimID)
          if (status /= nf90_noerr) call handle_err(status, ' LevDimID')

           status = nf90_inq_dimid(ncid, "nlats", LatDimID)
           if (status /= nf90_noerr) call handle_err(status, ' LatDimID')

          status = nf90_inq_dimid(ncid, "nlons", LonDimID)
          if (status /= nf90_noerr) call handle_err(status, ' LonDimID')


            status = nf90_inquire_dimension(ncid, LevDimID, len = nz)
            status = nf90_inquire_dimension(ncid, LatDimID, len = ny)
            status = nf90_inquire_dimension(ncid, LonDimID, len = nx)
            status = nf90_inquire_dimension(ncid, RecordDimID, len = nt)
            
            !print *, 'Gar Read dims',nz,ny,nx,nt

            if (nt.ne.nta_hwm .or. nx.ne.nxa_hwm .or. nz.ne.nza_hwm .or. ny.ne.nya_hwm) then
              print *, ' allocated nx-nt       ', nxa_hwm, nya_hwm,nza_hwm, nta_hwm
              print *, ' written in file nx-nt ', nx, ny, nz, nt
              print *, ' error in Analysis File ', me
              call mpi_quit(23999)
            endif
    
!          dates & coordinates
    
            !status = nf90_inq_varid(ncid, 'Date', date_id)
            !status = nf90_get_var(ncid, date_id, date_hwm)
            !status = nf90_inq_varid(ncid, 'Time', var_id)
            !status = nf90_get_var(ncid, var_id, date_sec_hwm)
    
            !status = nf90_inq_varid(ncid, 'Levels', var_id)
            !status = nf90_get_var(ncid, var_id, lev_hwm)

            !status = nf90_inq_varid(ncid, 'Latitude', var_id)
            !status = nf90_get_var(ncid, var_id, lat_hwm)

            !status = nf90_inq_varid(ncid, 'Longitude', var_id)
            !status = nf90_get_var(ncid, var_id, lon_hwm)
    
!           2D-arrays    count4(:4) = (/ plon, plat, plev, 1 /)
    
            start2(1:3) = [1, 1, Irec]
            start3(1:4) = [1, 1, 1, Irec]
            count2(1:3) = [nx, ny, 1]
            count3(1:4) = [nx, ny, nz, 1]

!           3D-arrays

            status = nf90_inq_varid(ncid, 'Zonal Wind', var_id)
            status = nf90_get_var(ncid, var_id, Uhwm, start=start3,
     &    count=count3)

            status = nf90_inq_varid(ncid, 'Meridional Wind', var_id)
            status = nf90_get_var(ncid, var_id, Vhwm, start=start3, 
     &       count=count3)

            !print*,'Gar close input analysis file on PE=', me
            status = nf90_close(ncid)
    
!           define nudging params to solve dT/dt = -(T-Ta)/tau
!           Tnudge = Tmodel + Ax*dt  Ax=(Ta-Tmodel)/Tau
!           Tnudge = Tmodel*(1.-dt/tau) + Ta*dt/tau = Wm*Tm + Wa*Ta
    
             wm_nudge_hwm(1:nza_hwm)=1.
             wa_nudge_hwm(1:nza_hwm)=0.
             rtau(1:nza_hwm) = 0.0

             knudge_hwm=30
             do k=10, knudge_hwm
                rtau(k) = 1./tau_6hr
                wa_nudge_hwm(k) = rtau(k)*tau_step
                wm_nudge_hwm(k) = 1.- wa_nudge_hwm(k)           
             enddo
             ktrans_hwm=40
             do k=1+knudge_hwm, ktrans_hwm
               !decay=exp(-abs(k-knudge)*0.2)
               decay = abs(ktrans_hwm-k)/float(dknudge)
               rtau(k) = rtau(knudge_hwm)*decay
               wa_nudge_hwm(k) = wa_nudge_hwm(knudge_hwm)*decay
               wm_nudge_hwm(k) = 1.- wa_nudge_hwm(k) 
             enddo
             !print *, 'Gar1 : ',nx,ny,nz,nt,knudge_hwm,dknudge,ktrans_hwm
             !print *, 'Gar2 : ',lev_hwm
             !print *, 'Gar3 : ',rtau
             !print *, 'Gar4 : ', wa_nudge_hwm 
             !print *, 'Gar5 : ', wm_nudge_hwm     
             !print *, 'Gar6 : ',Uhwm
             !print *, 'Gar7 : ',Vhwm
             END SUBROUTINE IDEA_INIT_IAU_HWM



!
         subroutine IDEA_IAU_XZY_DECOMP(A, B, nx, ny, nz)
         implicit none
         integer :: nx, ny,  nz
         real    :: A(nx, ny, nz)
         real    :: B(nx*nz, ny)
         integer :: i, k, j, ixz
         do j=1, ny
            do k=1, nz
               do  i=1,nx
!
!repack 3D=> 2D array
!
              B((k-1)*nx+i,j) = A(i,j,k)
               enddo
            enddo
         enddo
         END subroutine IDEA_IAU_XZY_DECOMP
!
       END MODULE idea_iau_hwm
!
       SUBROUTINE read_analysis(CDATE8, IREC, me)
       use netcdf
       use idea_iau_hwm, only : handle_err      
       use idea_iau_hwm, only : nxa_hwm, nya_hwm, nza_hwm, nta_hwm
       use idea_iau_hwm, only : Uhwm, Vhwm    
       use idea_iau_hwm, only : lev_hwm, lat_hwm, lon_hwm
       use idea_iau_hwm, only : Dir_list_hwm, Fend_hwm, File_anl_hwm, File_pref_hwm


       implicit none
       integer, intent(in) :: CDATE8, IREC, me

     
      integer    :: nx, ny, nz, nt
      integer    :: status, ncid
      integer    :: iLevDimID, LevDimID, LonDimID, LatDimID, RecordDimID
      integer    :: DateDimID
      integer    :: date_id, var_id
      integer    :: start2(3), count2(3)
      integer    :: start3(4), count3(4)
      integer    :: date_hwm(nta_hwm), date_sec_hwm(nta_hwm)
      integer    :: DATE_4  
 
      character(len=132)  :: Ftrial
      character(len=8)    :: S8_DATE
      character(len=4)    :: SY_DATE
      write(S8_DATE, fmt='(I8.8)') CDATE8
!
      !print *, S8_DATE, ' STR-date '
      
      DATE_4=CDATE8/10000
      write(SY_DATE, fmt='(I4.4)') DATE_4

      !File_anl_hwm=trim(Dir_list_hwm)//'/'//trim(FILE_PREF_hwm)//trim(S8_DATE)//trim(Fend_hwm)
!!      File_anl=trim(Dir_list)//'/'//trim(SY_DATE)//'/MERRA_new_Zres'//trim(S8_DATE)//trim(Fend)
!!!!!      File_anl=trim(Dir_list)//trim(S8_DATE)//trim(Fend)
      File_anl_hwm=trim(Dir_list_hwm)//'/'//trim(FILE_PREF_hwm)//trim(S8_DATE)//trim(Fend_hwm)
      if (me==0) then
        print *, Irec, ' analysis-IRec '
        print *, File_anl_hwm, ' analysis-Garima File_anl '
      endif
      status = nf90_open(File_anl_hwm, nf90_nowrite, ncid)

      if (status /= nf90_noerr) call handle_err(status, ' nf90_open')
      status = nf90_inq_dimid(ncid, 'ntime',RecordDimID)
      if (status /= nf90_noerr) call handle_err(status, ' RecordDimID')
      status = nf90_inq_dimid(ncid, "nlevs", LevDimID)
      if (status /= nf90_noerr) call handle_err(status, ' LevDimID')

       status = nf90_inq_dimid(ncid, "nlats", LatDimID)
       if (status /= nf90_noerr) call handle_err(status, ' LatDimID')

      status = nf90_inq_dimid(ncid, "nlons", LonDimID)
      if (status /= nf90_noerr) call handle_err(status, ' LonDimID')


        status = nf90_inquire_dimension(ncid, LevDimID, len = nz)
        status = nf90_inquire_dimension(ncid, LatDimID, len = ny)
        status = nf90_inquire_dimension(ncid, LonDimID, len = nx)
        status = nf90_inquire_dimension(ncid, RecordDimID, len = nt)
        if (me == 0) then
        if (nt.ne.nta_hwm .or. nx.ne.nxa_hwm .or. nz.ne.nza_hwm .or. ny.ne.nya_hwm) then
          print *, ' allocated nx-nt       ', nxa_hwm, nya_hwm, nza_hwm, nta_hwm
          print *, ' written in file nx-nt ', nx, ny, nz, nt
          print *, ' error in Analysis File '
         endif
        endif
       
!
! dates & coordinates
!
        !status = nf90_inq_varid(ncid, 'Date', date_id)
        !status = nf90_get_var(ncid, date_id, date_hwm)
        !status = nf90_inq_varid(ncid, 'Time', var_id)
        !status = nf90_get_var(ncid, var_id, date_sec_hwm)
!
        !status = nf90_inq_varid(ncid, 'Levels', var_id)
        !status = nf90_get_var(ncid, var_id, lev_hwm)

        !status = nf90_inq_varid(ncid, 'Latitude', var_id)
        !status = nf90_get_var(ncid, var_id, lat_hwm)

        !status = nf90_inq_varid(ncid, 'Longitude', var_id)
        !status = nf90_get_var(ncid, var_id, lon_hwm)
!
! 2D-arrays
!
        start2(1:3) = [1, 1, Irec]
        start3(1:4) = [1, 1, 1, Irec]
        count2(1:3) = [nx, ny, 1]
        count3(1:4) = [nx, ny, nz, 1]

! 3D-arrays

        status = nf90_inq_varid(ncid, 'Zonal Wind', var_id)
        status = nf90_get_var(ncid, var_id, Uhwm, start=start3, count=count3)

        status = nf90_inq_varid(ncid, 'Meridional Wind', var_id)
        status = nf90_get_var(ncid, var_id, Vhwm, start=start3, count=count3)

        !print*,'close input analysis file'
        status = nf90_close(ncid)


        !print *
        !print *, Irec
        !print *

!         print *, maxval(Tanl), minval(tanl), ' Garima-Tanl'
!         print *, maxval(Uhwm), minval(Uhwm), ' Garima-Uanl'
!         print *, maxval(Vhwm), minval(Vhwm), ' Garima-Vanl'
!         print *, maxval(Psanl), minval(Psanl), ' VAY-PSanl'
!         print *, maxval(Tsanl), minval(Tsanl), ' VAY-TSanl'
    
       END  SUBROUTINE read_analysis
!========================================================================
! layout1.f
!      integer           nodes, nodes_comp,nodes_io,
!     x                  me,
!     x                  ls_dim,
!     x                  ls_max_node,
!     x                  lats_dim_r,
!     x                  lats_dim_ext,
!     x                  lats_node_r,
!     x                  lats_node_r_max,
!     x                  lats_node_ext,
!     x                  ipt_lats_node_r,
!     x                  ipt_lats_node_ext,
!     x			 lonf, latg,
!     x                  len_trio_ls, len_trie_ls,
!     x                  me_l_0
! Where and how add analysis increments:
!
! gfs_physics_output.f  : Create "ADD_BUNDLE_TO_WRITE" as 2D-fields "SFC"
!  gfs_phy_states_mod.f : Import Phys/State Vars 2D-3D =>   GOCART
!<<<<<<
! 
! gfs_physics_initialize_mod.f: call fix_fields(LONSPERLAR, GLOBAL_LATS_R,ORO/HPRIME/PL-chem
!
! gfs_physics_run_mod.f       :      call do_physics_one_step
!
! do_physics_one_step.f       : call gcycle; CALL GLOOPR ( grid_fld, g3d_fld, aoi_fld, lats_nodes_r
!                             : call gloopb (grid_fld,     g3d_fld,       sfc_fld,
!grep -i call gloopb.f        :        CALL GW_unified_init(levs, dtphys_gw, me, ak5, bk5)
!                                      call read_wam_f107_kp_txt
!                                      call get_gmao_iau_forcing
!                                      call apply_gmao_iau_forcing
! Check read NC-files on   me=0 =>  MPI_BCAST(6-vars)
!                          time-variable read records/days
!
!========================================================================
