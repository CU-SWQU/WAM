! input forcing parameter module
      module wam_ifp_mod

      use comio
      use wam_ifp_class
      use mpi_def, only: MPI_COMM_ALL, MPI_INFO_NULL
      use layout1, only: me

      implicit none

      type(param_t) :: params
      type(farr_t)  :: farr
! Legacy variables
      real    :: f107_fix, f107d_fix, kp_fix
      real    :: swpcf107_fix, swpcf107d_fix, swpckp_fix
!

      contains

      subroutine read_ifp

      character(len=19), parameter :: filename = "input_parameters.nc"
      class(COMIO_T), allocatable :: io
      integer, parameter      :: fmt =  COMIO_FMT_PNETCDF
      integer, pointer :: dims(:)

      call COMIO_Create(io, fmt, &
                        comm=MPI_COMM_ALL, &
                        info=MPI_INFO_NULL)

      call dealloc()

      farr % default_euvfac = .false.
      farr % default_jhfac  = .false.

      call io % open(filename, "r")
      call io % description("skip", params % skip)
      call io % description("ifp_interval", params % ifp_interval)

      call io % domain("f107", dims)
      call alloc(dims(1))
      call io % read("f107",    farr % f107)
      call io % read("f107d",   farr % f107d)
      call io % read("kp",      farr % kp)
      call io % read("kpa",     farr % kpa)
      call io % read("nhp",     farr % nhp)
      call io % read("nhpi",    farr % nhpi)
      call io % read("shp",     farr % shp)
      call io % read("shpi",    farr % shpi)
      call io % read("swden",   farr % swden)
      call io % read("swang",   farr % swang)
      call io % read("swvel",   farr % swvel)
      call io % read("swbz",    farr % swbz)
      call io % read("swbt",    farr % swbt)

      call io % read("euvfac", farr % euvfac)
      if (io % err % check()) farr % default_euvfac = .true.
      call io % read("jhfac",   farr % jhfac)
      if (io % err % check()) farr % default_jhfac  = .true.

      call io % close()

      end subroutine read_ifp

      subroutine alloc(dim)
        integer, intent(in) :: dim

        if (.not.allocated(farr%f107))   allocate(farr%f107 (dim))
        if (.not.allocated(farr%f107d))  allocate(farr%f107d(dim))
        if (.not.allocated(farr%kp))     allocate(farr%kp   (dim))
        if (.not.allocated(farr%kpa))    allocate(farr%kpa  (dim))
        if (.not.allocated(farr%nhp))    allocate(farr%nhp  (dim))
        if (.not.allocated(farr%nhpi))   allocate(farr%nhpi (dim))
        if (.not.allocated(farr%shp))    allocate(farr%shp  (dim))
        if (.not.allocated(farr%shpi))   allocate(farr%shpi (dim))
        if (.not.allocated(farr%swden))  allocate(farr%swden(dim))
        if (.not.allocated(farr%swvel))  allocate(farr%swvel(dim))
        if (.not.allocated(farr%swang))  allocate(farr%swang(dim))
        if (.not.allocated(farr%swbz))   allocate(farr%swbz (dim))
        if (.not.allocated(farr%swbt))   allocate(farr%swbt (dim))
        if (.not.allocated(farr%jhfac))  allocate(farr%jhfac(dim))
        if (.not.allocated(farr%euvfac)) allocate(farr%euvfac(dim))

      end subroutine alloc

      subroutine dealloc()

        if (allocated(farr%f107))   deallocate(farr%f107)
        if (allocated(farr%f107d))  deallocate(farr%f107d)
        if (allocated(farr%kp))     deallocate(farr%kp)
        if (allocated(farr%kpa))    deallocate(farr%kpa)
        if (allocated(farr%nhp))    deallocate(farr%nhp)
        if (allocated(farr%nhpi))   deallocate(farr%nhpi)
        if (allocated(farr%shp))    deallocate(farr%shp)
        if (allocated(farr%shpi))   deallocate(farr%shpi)
        if (allocated(farr%swden))  deallocate(farr%swden)
        if (allocated(farr%swvel))  deallocate(farr%swvel)
        if (allocated(farr%swang))  deallocate(farr%swang)
        if (allocated(farr%swbz))   deallocate(farr%swbz)
        if (allocated(farr%swbt))   deallocate(farr%swbt)
        if (allocated(farr%jhfac))  deallocate(farr%jhfac)
        if (allocated(farr%euvfac)) deallocate(farr%euvfac)

      end subroutine dealloc
! legacy code below, not sure this is still needed
!==========================================================
! Below two service subs to disable "read_wam_f107_kp_txt"
! during model tune-ups
!==========================================================
      subroutine fix_spweather_data
!=======================================================================
!VAY 2016: This is temporal "substitue" for "read_wam_f107_kp_txt"
!    with fixed Kp and F107 data to work with long-term WAM run
! TO DO "advance_solar" KP-F107 drivers with WAM calendar
!=======================================================================
      swpcf107_fix = 100.
      swpckp_fix   = 1.
      swpcf107d_fix = swpcf107_fix

      f107_fix = 100.
      kp_fix   = 1.
      f107d_fix = f107_fix
      end subroutine fix_spweather_data
!
      subroutine read_spweather_real_data
!=======================================================================
!VAY 2016: This is temporal "substitue" for "read_wam_f107_kp_txt"
!    with fixed Kp and F107 data to work with long-term WAM run
! TO do "advance_solar" KP-F107 drivers with WAM calendar
!=======================================================================
      f107_fix = 100.
      kp_fix   = 1.
      f107d_fix = f107_fix

      end subroutine read_spweather_real_data
!

      end module wam_ifp_mod
