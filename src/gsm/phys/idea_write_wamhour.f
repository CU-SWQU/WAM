
          Module idea_ncout_phys_1hr
!
! fixed 3D fields
!         use layout1                            !lats_node_r, lats_node_r_max
          use resol_def, only : lonr, latr, levs !lonf, latg
!
! accumulation id performed at lonr, lats_node_r  with grid_fld%ps( lonr, lats_node_rmax)
!
          real,allocatable, dimension(:, :)    ::  wps, wphis
          real,allocatable, dimension(:, :, :) :: wU, wV, wT
          real,allocatable, dimension(:, :, :) :: wQ, wO3, wOP, wO2
          real,allocatable, dimension(:, :, :) :: Wp3d, Wz3d, Wdp

!
          real,allocatable, dimension(:, :, :)    ::  w3d
          real,allocatable, dimension(:, :)       ::  w2d
!


!===============================================================================
! Total # of fields: 7-daily
! 
!================================================================================
          contains
!
          subroutine idea_init_nc_1hr
          use resol_def, only : lonr, latr, levs
          use layout1,only : n1 => lats_node_r, nmax=>lats_node_r_max
          implicit none   
        
          allocate(wps(lonr, n1), wphis(lonr,n1), w2d(lonr,n1) )
!   
          allocate(wu(lonr, n1, levs), wv(lonr,n1,levs))
          allocate(w3d(lonr, n1, levs),  wT(lonr,n1,levs))
          allocate(wQ(lonr, n1, levs), wO3(lonr,n1,levs))
          allocate(wOP(lonr, n1, levs), wO2(lonr,n1,levs))
          allocate(wdP(lonr, n1, levs), wP3d(lonr,n1,levs))   ! pressures
          allocate(wZ3d(lonr, n1, levs))                      ! geo-height in KM


          end subroutine idea_init_nc_1hr

          end Module idea_ncout_phys_1hr
!
       subroutine idea_copy_inst(grid_fld, ind_lats, ind_lons)
!
!23456
!
          use layout1, only : lats_node_r, lats_node_r_max,  me
          use resol_def, only : lonr, latr, levs                  !lonf, latg
          use gfs_physics_gridgr_mod, only : Grid_Var_Data
!
          use  idea_ncout_phys_1hr
          use  wam_pass_diag_types,    only : Gis_wam     
     
          implicit none
!input
          TYPE(Grid_Var_Data)     :: grid_fld     ! dims: lonr, lats_node_r_max, levs, gis_phy%ntrac in gfs_physics_initialize_mod.f
!
          integer, dimension(lats_node_r) :: ind_lats, ind_lons

!locals
          integer :: n1, n2 
!
          integer ::  i,j,k         
          n1 = lonr
          n2 = lats_node_r
      
!
!
   
!$omp parallel do private(i,j)
      do j=1,n2
        do i=1,lonr
        wps(i,j)= grid_fld%ps(i, j)
        wphis(i,j)=grid_fld%z(i, j)
        enddo
      enddo
!$omp parallel do private(i,j,k)
      do k=1, levs
!
      do j=1,n2
        do i=1,lonr
        wu(i,j,k)=grid_fld%U(i, j,k)
        wv(i,j,k)=grid_fld%V(i, j,k)
        wt(i,j,k)=grid_fld%T(i, j,k)
        wQ(i,j,k)=grid_fld%tracers(1)%flds(i,j,k)
        wO3(i,j,k)=grid_fld%tracers(2)%flds(i,j,k)
        wOP(i,j,k)=grid_fld%tracers(4)%flds(i,j,k)
        wO2(i,j,k)=grid_fld%tracers(5)%flds(i,j,k)
        wdp(i,j,k)=grid_fld%dp(i, j,k)
        wp3d(i,j,k)=grid_fld%p(i, j,k)
        wz3d(i,j,k)=gis_wam%zgkm(i, j,k)
        enddo
      enddo
      ENDDO
!
       END subroutine idea_copy_inst
!
!
       subroutine idea_write_wamhour( global_lats_r,lonsperlar,
     &        fhour,idate,Curr_NC_WAMDAY, Curr_NC_WAMDAYHR)
!
! Prototype from wam_nc_output16.f
!
        use resol_def,   ONLY: latr, levs, levp1, lonr
        use layout1,     ONLY: me, nodes, lats_node_r
        use mpi_def,     ONLY: liope, info, mpi_comm_all, 
     &                                   mc_comp, mpi_comm_null
        USE machine,   ONLY: kind_io4, kind_io8
!
        use  idea_ncout_phys_1hr
!
        use gg_def,         only : colrad_r      ! latr
        use coordinate_def, only : ak5,bk5       ! levs+1
        use netcdf
        use idea_lat_gaus, only : lat_gaus_t62
!
	use idea_gdas_calendar, only : hist_gdas, Dirnc_case    !hourly files
	
          implicit none

          real        :: fhour
          integer     :: idate(4)
          integer     :: Curr_NC_WAMDAYHR            !, Wam_daysteps
          integer     :: Curr_NC_WAMDAY    

!
          integer, dimension(latr)   :: global_lats_r, lonsperlar
!
          integer         :: IOPROC, node, lat, ierr
!
      real, allocatable   :: tmp3d(:,:,:), tmp2d(:,:)
!=======================================================================
! need to be passed: hyam, hybm, hyai, hybi, lons_deg, lats_deg, levi
!=======================================================================
      character(len=10)    :: S8_DATE
      character(len=8)    :: file_wam='hist_hr_'
!123456789012345678901
!hist_hr_2016012500.nc
      character(len=21)   :: file_hist               
      character(len=10)   :: File_day='2016012500' 
      character(len=3)    :: Fend ='.nc'
      character(len=121)  :: file_hist_das
      character(len=100)  :: dir_case='/scratch3/NCEPDEV/swpc/scrub/Valery.Yudin/NC_DAS/'
!
!
! Coordinate arrays
!
!
      real, dimension(levs)   :: pmbm, hyam, hybm 
      real, dimension(levs+1) :: pmbi, hyai, hybi 
      real                    :: lonwam(lonr)
      real                    :: latwam(latr)
      real                    :: Pref
!
      integer    :: status, ncid, iernc    
      real       ::  time  
  
      integer :: NxDimID, NyDimID, NzDimID, NzIDimID, NtDimID
      integer :: scalDimID, Vid_sp0
      integer :: VidLat, VidLon, Vidlev, VidIlev
      integer :: VidHyam, VidHybm, VidHyai, VidHybi
      integer :: VidTime, VidDate, VidDsec
!
      integer :: VidHs, VidPs
      integer :: VidT, VidQ, VidV, VidU, VidZgkm
      integer :: Vidqo2, Vidqo3, Vidqo  

      integer :: Vid

      integer :: Datesec, Ymd

         integer :: start1(1), count1(1)
         integer :: start(4), count(4)
         integer :: nz, ny, nx
!
         integer :: start3(4), count3(4)
         integer :: ihr
!
         integer :: i, j,k
!
             ihr =1                        ! daily-record
          YMD = Curr_NC_WAMDAY
!
! define IOPROC
!
          IOPROC =nodes-1                  !like in SFC-RSTSR  wrtout_physics.f  
!
        call mpi_barrier(mpi_comm_all,ierr)
!      t3  = rtc()
!      call mpi_barrier(mpi_comm_all,ierr)
!      t4  = rtc()
!      tba = t4-t3


          nz =levs
          ny = latr
          nx = lonr

          start=(/1,1,1, ihr/)
          count=(/nz, ny, nx,1/)

          start1(1) =1
          count1(1) =1


      allocate ( tmp2d (lonr,latr), tmp3d (lonr,latr, levs) )

      if(me == ioproc) then   
!
! Create name of the HIST-file
!
    
      write(S8_DATE, fmt='(I10.10)') Curr_NC_WAMDAYHR

      print *, S8_DATE, ' NC_STR-date '
      File_hist=trim(file_wam)//trim(S8_DATE)//trim(Fend)
      print *, File_hist, ' VAY File_hist NC-file '
!      File_hist_das=trim(dir_case)//trim(file_wam)//trim(S8_DATE)//trim(Fend)
      File_hist_das=trim(dirnc_case)//trim(hist_gdas)//trim(S8_DATE)//trim(Fend)      
!
! open-File_hist
!23456
      ierNC=NF90_OPEN(trim(File_hist_das), nf90_write, ncid)
      status = nf90_create(trim(File_hist_das), nf90_clobber, ncid)
!
! Dimensions 
!
      status = nf90_def_dim(ncid, 'lon',  lonr,          NxDimID)
      status = nf90_def_dim(ncid, 'lat',  latr,          NyDimID)
      status = nf90_def_dim(ncid, 'time', NF90_UNLIMITED,NtDimID)
      status = nf90_def_dim(ncid, 'scalar', 1,         scalDimID)

      status = nf90_def_dim(ncid, 'lev',  levs,          NzDimID)
      status = nf90_def_dim(ncid, 'ilev', levs+1,        NzIDimID)
!
! Create vars and attributes .........
! Coordinates
!
      status = nf90_def_var(ncid, 'time',nf90_float, (/ NtDimID /), VidTime)
          status = nf90_put_att(ncid,VidTime, 'long_name', 'data-time')
          status = nf90_put_att(ncid,VidTime, 'units', 'days since 2012-01-01 00:00:00')
!
          status = nf90_def_var(ncid, 'date',     nf90_INT, (/ NtDimID /), VidDate)
          status = nf90_put_att(ncid,VidDate, 'long_name', 'current date (YYYYMMDD)')

          status = nf90_def_var(ncid, 'datesec',     nf90_INT, (/ NtDimID /), VidDsec)
          status = nf90_put_att(ncid,VidDsec, 'long_name', 'current seconds of current date')
          status = nf90_put_att(ncid,VidDsec, 'units', 'seconds')


!
          status = nf90_def_var(ncid, 'lat',     nf90_float, (/ NyDimID /),  VidLat)
          status = nf90_put_att(ncid,Vidlat, 'long_name', 'data-latitude')
          status = nf90_put_att(ncid,Vidlat, 'units', 'degrees_north')
!
          status = nf90_def_var(ncid, 'lon',     nf90_float, (/ NxDimID /), Vidlon)
          status = nf90_put_att(ncid,Vidlon, 'long_name', 'data-longitude')
          status = nf90_put_att(ncid,Vidlon, 'units', 'degrees_east')
!
          status = nf90_def_var(ncid, 'lev',     nf90_float, (/ NzDimID /), VidLev)
          status = nf90_put_att(ncid,Vidlev, 'units', 'level')
          status = nf90_put_att(ncid,Vidlev, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
          status = nf90_put_att(ncid,Vidlev, 'long_name', 'hybrid level at midpoint (1000*(A+B))')
!          status = nf90_put_att(ncid,Vidlev, 'units', 'mb')
          status = nf90_put_att(ncid,Vidlev, 'positive', 'down') !   ncdf_attput, id, ilev_id, 'positive', 'down'
          status = nf90_def_var(ncid, 'hyam',     nf90_float, (/ NzDimID /), VidHyam)
          status = nf90_put_att(ncid,VidHyam, 'long_name', 'hybrid level at midpoint A')
          status = nf90_put_att(ncid,VidHyam, 'units', ' dimensionless ')

          status = nf90_def_var(ncid, 'hybm',     nf90_float, (/ NzDimID /), VidHybm)
          status = nf90_put_att(ncid,VidHybm, 'long_name', 'hybrid level at midpoint B')
          status = nf90_put_att(ncid,VidHybm, 'units', ' dimensionless ')
!
          status = nf90_def_var(ncid, 'ilev',     nf90_float, (/ NzIDimID /), VidiLev)
          status = nf90_put_att(ncid,VidIlev, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
          status = nf90_put_att(ncid,VidIlev, 'long_name', 'hybrid level at interface (1000*(A+B))')
          status = nf90_put_att(ncid,VidIlev, 'units', 'lev')
          status = nf90_put_att(ncid,VidIlev, 'positive', 'down')
          status = nf90_def_var(ncid, 'hyai',     nf90_float, (/ NzIDimID /), VidHyaI)
          status = nf90_put_att(ncid,VidHyam, 'long_name', 'hybrid level at interface Ai')
          status = nf90_put_att(ncid,VidHyam, 'units', ' dimensionless ')

          status = nf90_def_var(ncid, 'hybi',     nf90_float, (/ NzIDimID /), VidHybI)
          status = nf90_put_att(ncid,VidHybI, 'long_name', 'hybrid level at interface Bi')
          status = nf90_put_att(ncid,VidHybI, 'units', ' dimensionless ')
!
          status = nf90_def_var(ncid, 'P0',     nf90_float, (/ scalDimID /), Vid_sp0)
          status = nf90_put_att(ncid,Vid_sp0, 'long_name', ' Reference surface pressure')
          status = nf90_put_att(ncid,Vid_sp0, 'units', ' Pa ')

!          status = nf90_def_var(ncid, 'day_obsed',   nf90_double, (/ NtDimID /),  VIDdate)
!============================================
!=========================================== 
!23456
      status = nf90_def_var(ncid, 'HS', nf90_float, (/ NxDimId, NyDimId, NtDimID /), VidHs)
      status = nf90_put_att(ncid,VidT, 'long_name', ' Surface elevation ')
      status = nf90_put_att(ncid,VidT, 'units', 'm')
!
      status = nf90_def_var(ncid, 'PS', nf90_float, (/ NxDimId, NyDimId, NtDimID /), VidHs)
      status = nf90_put_att(ncid,VidT, 'long_name', ' Surface pressure ')
      status = nf90_put_att(ncid,VidT, 'units', ' Pa ')
!
      status = nf90_def_var(ncid, 'U', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidU)
      status = nf90_put_att(ncid,VidU, 'long_name', ' Zonal wind ')
      status = nf90_put_att(ncid,VidU, 'units', 'm/s')
!
      status = nf90_def_var(ncid, 'V', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidV)
      status = nf90_put_att(ncid,VidV, 'long_name', ' Meridional wind ')
      status = nf90_put_att(ncid,VidV, 'units', 'm/s' )
!
      status = nf90_def_var(ncid, 'T', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidT)
      status = nf90_put_att(ncid,VidT, 'long_name', ' Temperature')
      status = nf90_put_att(ncid,VidT, 'units', 'K') 
!
      status = nf90_def_var(ncid, 'P3D', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidT)
      status = nf90_put_att(ncid,VidT, 'long_name', ' Pressure-3D')
      status = nf90_put_att(ncid,VidT, 'units', 'Pa')
!
      status = nf90_def_var(ncid, 'DP', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidT)
      status = nf90_put_att(ncid,VidT, 'long_name', ' P-thickness')
      status = nf90_put_att(ncid,VidT, 'units', 'Pa')
!
      status = nf90_def_var(ncid, 'Q', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidQ)
      status = nf90_put_att(ncid,VidQ, 'long_name', ' H2O, mmr')
      status = nf90_put_att(ncid,VidQ, 'units', ' kg/kg ') 

      status = nf90_def_var(ncid, 'O3', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidQo3)
      status = nf90_put_att(ncid,VidQo3, 'long_name', ' O3, mmr')
      status = nf90_put_att(ncid,VidQo3, 'units', ' kg/kg ') 

      status = nf90_def_var(ncid, 'O2', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidQo2)
      status = nf90_put_att(ncid,VidQo2, 'long_name', ' O2, mmr')
      status = nf90_put_att(ncid,VidQo2, 'units', ' kg/kg ') 
!
      status = nf90_def_var(ncid, 'O', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidQo)
      status = nf90_put_att(ncid,VidQo, 'long_name', ' O, mmr')
      status = nf90_put_att(ncid,VidQo, 'units', ' kg/kg ') 
!
      status = nf90_def_var(ncid, 'ZGKM', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidZgkm)
      status = nf90_put_att(ncid,VidZgkm, 'long_name', ' Geo-Height')
      status = nf90_put_att(ncid,VidZgkm, 'units', ' km ') 
!
! Tides 3-vars (U,V,T) * 6-coef = 18-VARS
!

! End Names & Attributes
!
      status = nf90_enddef(ncid)

!========================================
! Put Ymd and Model Coordiantes
!========================================

      time = 24.0           ! daily mean
      Datesec = ihr * 3600     
!
! Grids:     see glats_physics.f
!
      do i=1, lonr
         lonwam(i) = 0.+ (i-1)*360./float(lonr)
      enddo

      do i=1, latr
!         latwam(i) = 90.- (i-1)*180./float(latr-1)
         latwam(i)= lat_gaus_t62(i)     !colrad_r(i)*45./atan(1.)      ! Gaussian Lats
      enddo
!
! need to define hyai/hybi on edges  with zeroes at the top
! VG-read/write in the dynamical core
!
      do k=1, levs+1
        hyai(k)=ak5(k)*1.e-2         !*1000/(Ps = 1.e5) = 1.e-2
        hybi(k)=bk5(k)
        pmbi(k) =1.e5*(hyai(k)+hybi(k))
      enddo

      do k=1, levs
        pmbm(k)=0.5*(pmbi(k)+pmbi(k+1))
        hyam(k)=0.5*(hyai(k)+hyai(k+1))
        hybm(k)=0.5*(hybi(k)+hybi(k+1))
      enddo
!
! Grids completed
!
      status = nf90_put_var(ncid,  VidDate,   YMD)      !, start = start1,  count = count1)
      status = nf90_put_var(ncid,  VidDsec,   Datesec) !
      status = nf90_put_var(ncid,  VidTime,   Time)
 
      status = nf90_put_var(ncid, VidLon, LonWam)
      status = nf90_put_var(ncid, VidLat, LatWam)

      status = nf90_put_var(ncid, VidTime, Time)
      status = nf90_put_var(ncid, VidLev, Pmbm)
      status = nf90_put_var(ncid, VidHyam, hyam)
      status = nf90_put_var(ncid, VidHybm, hybm)

      status = nf90_put_var(ncid, VidiLev, Pmbi)
      status = nf90_put_var(ncid, VidHyai, hyai)
      status = nf90_put_var(ncid, VidHybi, hybi)
      Pref = 1.e5
      status = nf90_put_var(ncid, Vid_sp0, Pref)


      ENDIF     ! (IOPROC)
!=================================================
! Now Gain 2D/3D and write on ME=IOPROC
!=================================================

        call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp2d, 1, wPs, GLOBAL_LATS_R,LONSPERLAR)
!
       call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then    
        iernc=nf90_inq_varid( ncid, 'PS', vid )
        iernc= nf90_put_var( ncid, vid, tmp2d)
      endif
        call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp2d, 1, wPhis, GLOBAL_LATS_R,LONSPERLAR)
      if(me == ioproc) then    
        iernc=nf90_inq_varid( ncid, 'HS', vid )
        iernc= nf90_put_var( ncid, vid, tmp2d)
      endif
!U
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wU, GLOBAL_LATS_R,LONSPERLAR)

       call mpi_barrier(mpi_comm_all,ierr)

      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'U', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
       endif
!V
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wV, GLOBAL_LATS_R,LONSPERLAR)

       call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'V', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif
!T
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wT, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'T', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif
!       print *, 'VAY-NCHIST-T ', maxval(tmp3d), minval(tmp3d)
!
! ADD P3D, DP, ZPHIL + TIDES
!
!P3d
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wP3d, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'P3D', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif
!DP
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wDp, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)

      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'DP', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif

!Z3D in kM
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wZ3d, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)

      if(me == ioproc) then
        iernc=nf90_inq_varid( ncid, 'ZGKM', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d )
      endif
!=================================================================
! Composition:    Q, Qo3, Qo, Qo2
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wQ, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'Q', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif
!O3
       call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wO3, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'O3', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif     
!OP
       call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wOP,  GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'O', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d )
      endif
!OP
       call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wO2, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'O2', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d )
      ENDIF
!
! close file
!
        if(me == ioproc) then
           status = nf90_close(ncid)

!
! data TYPES: NF90_BYTE, NF90_CHAR, NF90_SHORT, NF90_INT, NF90_FLOAT, and NF90_DOUBLE
!
           print *, ' NO-collects for NC-file '
        endif    ! IOPROC

        deallocate ( tmp2d, tmp3d )


         RETURN

          
!==============================================================================
!         on ME=0 open NC_FILE gather "ALL" PEs
!            write-out NC-daily FILE 36+7+2 instances
!            vs "7x24=168" having Hourly output  ~4 times less !!!!
!==============================================================================
!
! =>       call unsplit2d_phys(ioproc,wrkga,buffo,global_lats_r)
!
!          CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!
!          WRK3D & WRK2D => write into NC_file
!
!==============================================================================
          end subroutine idea_write_wamhour
!         