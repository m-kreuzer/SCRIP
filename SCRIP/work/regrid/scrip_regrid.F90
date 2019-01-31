!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this program is a short driver that reads a variable read from
!     a netCDF file (CF format), applies regridding according to an
!     interpolation file and stores the result in a netCDF file
!     (CF format)
!
!       based on:
!     CVS: $Id: scrip_test.f,v 1.6 2000/04/19 21:45:09 pwjones Exp $
!
!       which comes with the following license:
!     ------------------------------------------------------------------
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!***********************************************************************

      program scrip_regrid

!-----------------------------------------------------------------------

      use SCRIP_KindsMod   ! defines common data types
      use SCRIP_CommMod    ! for SCRIP parallel environment
      use SCRIP_ErrorMod   ! SCRIP error checking and logging
      use SCRIP_IOUnitsMod ! manages I/O units
      use SCRIP_InitMod    ! initialize SCRIP
      use SCRIP_ConfigMod  ! SCRIP configuration module
      use SCRIP_NetcdfMod  ! netcdf error handler
      use netcdf
      use constants    ! defines common constants
      use grids        ! module containing grid info
      use remap_vars   ! module containing remapping info
      use remap_mod    ! module containing remapping routines
      use remap_read   ! routines for reading remap files

      implicit none


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: &
         errorCode

      character (SCRIP_charLength) ::  &  ! namelist variables
              input_file,       &! name of input file containing data
              field_name,       &! input_file variable to be regridded
              interp_file,      &! filename of regridding weights
              output_file        ! name of file with regridded output


      character (SCRIP_charLength) ::  &
              map_name      ! name for mapping from grid1 to grid2

      integer (SCRIP_i4) ::   & ! netCDF ids for files and arrays
              ncstat, &
              nc_inputfile_id, nc_fieldname_id, &
              nc_outfile_id,  &
              nc_srcgrdcntrlat_id, nc_srcgrdcntrlon_id, &
              nc_dstgrdcntrlat_id, nc_dstgrdcntrlon_id, &
              nc_srcgrdimask_id, nc_dstgrdimask_id, &
              nc_srcgrdarea_id, nc_dstgrdarea_id, &
              nc_srcgrdfrac_id, nc_dstgrdfrac_id, &
              nc_srcarray_id, & 
              nc_dstarray_id

      integer (SCRIP_i4), dimension(:), allocatable ::  &
              nc_grid1size_id, nc_grid2size_id


      real (SCRIP_r8) ::   &
              integral_grid1, &
              integral_grid2

!-----------------------------------------------------------------------

      character (SCRIP_charLength) ::  &
                dim_name    ! netCDF dimension name

      integer (SCRIP_i4) :: i,j,n,imin,imax,idiff,  &
          iunit                  ! unit number for input configuration file

      integer (SCRIP_i4), dimension(:), allocatable ::  &
          grid1_imask, grid2_imask, grid2_count

      real (SCRIP_r8), dimension(:), allocatable :: &
          grid1_array_in,       &
          grid1_array_norm,     &
          grid2_array_out,      &
          grid2_array_norm,     &
          grid2_temp

      ! temporary arrays for saving data to ncfiles (real & integer)
      real (SCRIP_i4), dimension(:,:), allocatable ::  &
          grid1itmp2d, grid2itmp2d

      real (SCRIP_r8), dimension(:,:), allocatable ::  &
          grid1tmp2d, grid2tmp2d

      ! used for uniform output format
      character (20) :: format_str 

      character (12) ::   &
          rtnName = 'scrip_regrid'

!-----------------------------------------------------------------------
!
!     initialize SCRIP
!
!-----------------------------------------------------------------------

      errorCode = SCRIP_Success
      call SCRIP_CommInitMessageEnvironment
      call SCRIP_Initialize(errorCode)
      if (SCRIP_errorCheck(errorCode, rtnName, &
                           'error initializing SCRIP')) then 
        call SCRIP_RegridExit(errorCode)
      endif


!-----------------------------------------------------------------------
!
!  read input namelist   
!
!-----------------------------------------------------------------------

      call SCRIP_ConfigOpen(iunit, errorCode, 'scrip_regrid_in')
      if (SCRIP_ErrorCheck(errorCode, rtnName, &
          'error opening config file')) call SCRIP_RegridExit(errorCode)

      call SCRIP_ConfigRead(iunit, 'remapInputs',                      &
          'input_file', input_file, 'unknown', errorCode,              &
           outStringBefore=' -> file containing data for regridding: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_RegridExit(errorCode)


      call SCRIP_ConfigRead(iunit, 'remapInputs',                      &
             'field_name', field_name, 'unknown', errorCode,    &
             outStringBefore= ' -> variable/field to be regridded: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_RegridExit(errorCode)

      call SCRIP_ConfigRead(iunit, 'remapInputs',                      &
                    'interp_file', interp_file, 'unknown', errorCode,  &
                    outStringBefore=                                   &
          ' -> interpolation file containing weights for regridding: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_RegridExit(errorCode) 
      
      call SCRIP_ConfigRead(iunit, 'remapInputs',                      &
                'output_file', output_file, 'unknown', errorCode,      &
                outStringBefore=' -> output file with regridded data: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_RegridExit(errorCode) 

      call SCRIP_ConfigClose(iunit, errorCode)
      if (SCRIP_ErrorCheck(errorCode, rtnName, &
          'error closing config file')) call SCRIP_RegridExit(errorCode)

      print *,''
!-----------------------------------------------------------------------
!
!     read remapping data
!
!-----------------------------------------------------------------------

      ! reading grid specification data from interpolation file and 
      !   saving 'title' metadata to map_name
      call read_remap(map_name, interp_file, errorCode)
      if (SCRIP_ErrorCheck(errorCode, rtnName,  &
          'error reading remap file')) call SCRIP_RegridExit(errorCode)

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

      allocate (grid1_array_in    (grid1_size),  &
                grid1_array_norm  (grid1_size),  &
                grid1_imask       (grid1_size),  &
                grid2_array_out   (grid2_size),  &
                grid2_array_norm  (grid2_size),  &
                grid2_temp        (grid2_size),  &
                grid2_imask       (grid2_size),  &
                grid2_count       (grid2_size))

      where (grid1_mask)
        grid1_imask = 1
      elsewhere
        grid1_imask = 0
      endwhere
      where (grid2_mask)
        grid2_imask = 1
      elsewhere
        grid2_imask = 0
      endwhere

!-----------------------------------------------------------------------
!
!     setup a NetCDF file for output
!
!-----------------------------------------------------------------------

      !***
      !*** create netCDF dataset 
      !***

      ncstat = nf90_create (output_file, NF90_CLOBBER, nc_outfile_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
          'error opening output file')) call SCRIP_RegridExit(errorCode)


      ncstat = nf90_put_att(nc_outfile_id, NF90_GLOBAL, 'title',     &
                            map_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,         &
                                 'error writing output file title')) &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define grid size dimensions
      !***

      allocate( nc_grid1size_id(grid1_rank),  &
                nc_grid2size_id(grid2_rank))


      format_str = '(a9,i1)'

      do n=1,grid1_rank
        write(dim_name, format_str) 'grid1_dim',n
        ncstat = nf90_def_dim(nc_outfile_id, dim_name, &
                              grid1_dims(n), nc_grid1size_id(n))
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
            'error defining grid1 dims')) &
            call SCRIP_RegridExit(errorCode)
      end do

      do n=1,grid2_rank
        write(dim_name, format_str) 'grid2_dim',n
        ncstat = nf90_def_dim(nc_outfile_id, dim_name, &
                              grid2_dims(n), nc_grid2size_id(n))
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
            'error defining grid2 dims')) &
            call SCRIP_RegridExit(errorCode)
      end do
      
      !***
      !*** define grid center latitude array
      !***

      ncstat = nf90_def_var(nc_outfile_id, 'src_grid_center_lat', &
                            NF90_DOUBLE, nc_grid1size_id, &
                            nc_srcgrdcntrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
          'error defining src grid center lat')) &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_srcgrdcntrlat_id, &  
                            'units', 'radians')  
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
          'error adding units src grid center lat')) &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_def_var(nc_outfile_id, 'dst_grid_center_lat', &
                            NF90_DOUBLE, nc_grid2size_id, &
                            nc_dstgrdcntrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
          'error defining dst grid center lat'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_dstgrdcntrlat_id,  &
                            'units', 'radians')  
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                    'error adding units dst grid center lon'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define grid center longitude array
      !***

      ncstat = nf90_def_var(nc_outfile_id, 'src_grid_center_lon',  &
                            NF90_DOUBLE, nc_grid1size_id,  &
                            nc_srcgrdcntrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                        'error defining src grid center lon'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_srcgrdcntrlon_id, &
                            'units', 'radians')
      if (scrip_netcdferrorcheck(ncstat, errorcode, rtnname, &
                    'error adding units src grid center lon'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_def_var(nc_outfile_id, 'dst_grid_center_lon',  &
                            NF90_DOUBLE, nc_grid2size_id,  &
                            nc_dstgrdcntrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                        'error defining dst grid center lon'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_dstgrdcntrlon_id, &
                            'units', 'radians')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                    'error adding units dst grid center lon'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define grid mask
      !***

      ncstat = nf90_def_var(nc_outfile_id, 'src_grid_imask', NF90_INT, &
                            nc_grid1size_id, nc_srcgrdimask_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                             'error defining src grid mask'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_srcgrdimask_id, &
                            'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                             'error adding units to src mask'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_def_var(nc_outfile_id, 'dst_grid_imask', NF90_INT, &
                            nc_grid2size_id, nc_dstgrdimask_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error defining dst mask'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_dstgrdimask_id, &
                            'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error adding units to dst mask'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define grid area arrays
      !***

      ncstat = nf90_def_var(nc_outfile_id, 'src_grid_area', &
                            NF90_DOUBLE, nc_grid1size_id, &
                            nc_srcgrdarea_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error defining src area'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_srcgrdarea_id, &  
                            'units', 'square radians')  
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
          'error adding units src grid area')) &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_def_var(nc_outfile_id, 'dst_grid_area',  &
                            NF90_DOUBLE, nc_grid2size_id,  &
                            nc_dstgrdarea_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error defining dst area'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_dstgrdarea_id, &  
                            'units', 'square radians')  
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
          'error adding units dst grid area')) &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define grid fraction arrays
      !***

      ncstat = nf90_def_var(nc_outfile_id, 'src_grid_frac',  &
                            NF90_DOUBLE, nc_grid1size_id,  &
                            nc_srcgrdfrac_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error defining src grid frac'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_def_var(nc_outfile_id, 'dst_grid_frac',  &
                            NF90_DOUBLE, nc_grid2size_id,  &
                            nc_dstgrdfrac_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error defining dst grid frac'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define source array
      !***

      ncstat = nf90_def_var(nc_outfile_id, 'src_array',  &
                            NF90_DOUBLE, nc_grid1size_id,  &
                            nc_srcarray_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error defining src field'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define destination array
      !***

      ncstat = nf90_def_var(nc_outfile_id, 'dst_array',  &
                            NF90_DOUBLE, nc_grid2size_id,  &
                            nc_dstarray_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                       'error defining 1st order dst field'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** end definition stage
      !***

      ncstat = nf90_enddef(nc_outfile_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                       'error ending netCDF definition phase'))  &
          call SCRIP_RegridExit(errorCode)

!-----------------------------------------------------------------------
!
!     write some grid info
!
!-----------------------------------------------------------------------

      !***
      !*** write grid center latitude array
      !***

      allocate(grid1tmp2d(grid1_dims(1),grid1_dims(2)), &
               grid2tmp2d(grid2_dims(1),grid2_dims(2)))
      n = 0
      do j=1,grid1_dims(2)
      do i=1,grid1_dims(1)
         n = n+1
         grid1tmp2d(i,j) = grid1_center_lat(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_srcgrdcntrlat_id, &
                            grid1tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                       'error writing src grid center lat'))  &
          call SCRIP_RegridExit(errorCode)

      n = 0
      do j=1,grid2_dims(2)
      do i=1,grid2_dims(1)
         n = n+1
         grid2tmp2d(i,j) = grid2_center_lat(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_dstgrdcntrlat_id, &
                            grid2tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                       'error writing dst grid center lat'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** write grid center longitude array
      !***

      n = 0 
      do j=1,grid1_dims(2)
      do i=1,grid1_dims(1)
         n = n+1
         grid1tmp2d(i,j) = grid1_center_lon(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_srcgrdcntrlon_id, &
                            grid1tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing src grid lon'))  &
          call SCRIP_RegridExit(errorCode)

      n = 0
      do j=1,grid2_dims(2)
      do i=1,grid2_dims(1)
         n = n+1
         grid2tmp2d(i,j) = grid2_center_lon(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_dstgrdcntrlon_id, &
                            grid2tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing dst grid lon'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** write grid mask
      !***
      
      allocate(grid1itmp2d(grid1_dims(1),grid1_dims(2)), &
               grid2itmp2d(grid2_dims(1),grid2_dims(2)))
      n = 0
      do j=1,grid1_dims(2)
      do i=1,grid1_dims(1)
         n = n+1
         grid1itmp2d(i,j) = grid1_imask(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_srcgrdimask_id, &
                            grid1itmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing src grid mask'))  &
          call SCRIP_RegridExit(errorCode)

      n = 0
      do j=1,grid2_dims(2)
      do i=1,grid2_dims(1)
         n = n+1
         grid2itmp2d(i,j) = grid2_imask(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_dstgrdimask_id, &
                            grid2itmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing dst grid mask'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define grid area arrays
      !***

      n = 0
      do j=1,grid1_dims(2)
      do i=1,grid1_dims(1)
         n = n+1
         grid1tmp2d(i,j) = grid1_area(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_srcgrdarea_id, &
                            grid1tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing src grid area'))  &
          call SCRIP_RegridExit(errorCode)

      n = 0
      do j=1,grid2_dims(2)
      do i=1,grid2_dims(1)
         n = n+1
         grid2tmp2d(i,j) = grid2_area(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_dstgrdarea_id, &
                            grid2tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing dst grid area'))  &
          call SCRIP_RegridExit(errorCode)

      !***
      !*** define grid fraction arrays
      !***

      n = 0
      do j=1,grid1_dims(2)
      do i=1,grid1_dims(1)
         n = n+1
         grid1tmp2d(i,j) = grid1_frac(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_srcgrdfrac_id,  &
                            grid1tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing src grid frac'))  &
          call SCRIP_RegridExit(errorCode)

      n = 0
      do j=1,grid2_dims(2)
      do i=1,grid2_dims(1)
         n = n+1
         grid2tmp2d(i,j) = grid2_frac(n)
      end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_dstgrdfrac_id, &
                            grid2tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing dst grid frac'))  &
          call SCRIP_RegridExit(errorCode)


!-----------------------------------------------------------------------
!
!    reading field_name variable from input_file 
!       -> storing in grid1_array_in
!
!-----------------------------------------------------------------------

      ncstat = nf90_open(input_file, NF90_NOWRITE, nc_inputfile_id) 
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,   &
                                 'error opening input_file'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_inq_varid(nc_inputfile_id, field_name,  &
                                nc_fieldname_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,     &
                                'error getting field_name id'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_get_var(nc_inputfile_id, nc_fieldname_id,  &
                            grid1_array_in, start = (/1, 1, 1, 1/),  &
                            count = (/grid1_dims(1),grid1_dims(2),1,1/))
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,   &
                                 'error reading field_name'))  &
          call SCRIP_RegridExit(errorCode)

      ncstat = nf90_close(nc_inputfile_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,   &
                                 'error closing input_file'))  &
          call SCRIP_RegridExit(errorCode)


!-----------------------------------------------------------------------
!
!     test a first-order map from grid1 to grid2
!
!-----------------------------------------------------------------------

      !!! FLUX NORMALISATION
      ! convert grid1 variable from unit 'flux' to unit 'flux per area'
      grid1_array_norm = grid1_array_in/grid1_area
     
      if (map_type /= map_type_bicubic) then
        ! do the actual regridding
        call remap(grid2_array_norm, wts_map1, grid2_add_map1, &
                   grid1_add_map1, grid1_array_norm)
      else
        call SCRIP_ErrorSet(errorCode, rtnName, &
                            'bicubic interpolation not supported')
      endif

      ! apply remapping normalisation
      if (map_type == map_type_conserv .or. &
          map_type == map_type_particle) then
        select case (norm_opt)
        case (norm_opt_none)
          grid2_temp = grid2_frac*grid2_area
          where (grid2_temp /= zero)
            grid2_array_norm = grid2_array_norm/grid2_temp
          else where
            grid2_array_norm = zero
          end where
        case (norm_opt_frcarea)
        case (norm_opt_dstarea)
          where (grid2_frac /= zero)
            grid2_array_norm = grid2_array_norm/grid2_frac
          else where
            grid2_array_norm = zero
          end where
        end select
      end if


      format_str = "(A20,F14.4,F14.4)"

      print *,'First order mapping from grid1 to grid2:'
      print *,'----------------------------------------'
      write(*,format_str) 'Grid1 min,max: ', minval(grid1_array_norm), &
                                             maxval(grid1_array_norm)
      write(*,format_str) 'Grid2 min,max: ', minval(grid2_array_norm), &
                                             maxval(grid2_array_norm)
      print *,''

      !***
      !*** Conservation Test
      !***
      integral_grid1  = sum(grid1_array_norm*grid1_area*grid1_frac)
      integral_grid2  = sum(grid2_array_norm*grid2_area*grid2_frac)

      !format_str = "(A20,F40.20)"
      format_str = "(A20,ES25.16)"

      print *,'Conservation:'
      print *,'-------------'
      write(*,format_str) 'Grid1 Integral: ', integral_grid1
      write(*,format_str) 'Grid2 Integral: ', integral_grid2
      print *,''
      write(*,format_str) 'absolute error: ', &
                                      abs(integral_grid2-integral_grid1)
      write(*,format_str) 'relative error: ', & 
          abs(integral_grid2 - integral_grid1) / abs(integral_grid1)
      print *,''

      !!! FLUX NORMALISATION
      ! convert grid2 variable from unit 'flux per area' to unit 'flux'
      !grid1_array = grid1_array*grid1_area
      grid2_array_out = grid2_array_norm*grid2_area*grid2_frac

!-----------------------------------------------------------------------
!
!     write results to NetCDF file
!
!-----------------------------------------------------------------------

      n = 0
      do j=1,grid1_dims(2)
        do i=1,grid1_dims(1)
          n = n+1
          grid1tmp2d(i,j) = grid1_array_in(n)
        end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_srcarray_id, &
                            grid1tmp2d)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                                 'error writing src field'))  &
          call SCRIP_RegridExit(errorCode)

      n = 0
      do j=1,grid2_dims(2)
        do i=1,grid2_dims(1)
          n = n+1
          grid2tmp2d(i,j) = grid2_array_out(n)
        end do
      end do
      ncstat = nf90_put_var(nc_outfile_id, nc_dstarray_id, &
                            grid2tmp2d  )
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
                            'error writing remapped field'))  &
          call SCRIP_RegridExit(errorCode)


!-----------------------------------------------------------------------
!
!     close netCDF file
!
!-----------------------------------------------------------------------

      ncstat = nf90_close(nc_outfile_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                                 'error closing output file'))  &
          call SCRIP_RegridExit(errorCode)

!-----------------------------------------------------------------------
!
!     calculate some statistics
!
!-----------------------------------------------------------------------

      grid2_count       = zero
      grid2_array_norm  = zero
      grid2_temp        = zero


      format_str = "(A, I10)"

      print *,'Regridding statistics:'
      print *,'----------------------'

      write(*,format_str) ' number of sparse matrix entries       ', &
                                                          num_links_map1
      do n=1,num_links_map1
        grid2_count(grid2_add_map1(n)) =   &
        grid2_count(grid2_add_map1(n)) + 1
        if (wts_map1(1,n) > one .or. wts_map1(1,n) < zero) then
          grid2_array_norm(grid2_add_map1(n)) =   &
          grid2_array_norm(grid2_add_map1(n)) + 1
          grid2_temp(grid2_add_map1(n)) = max(abs(wts_map1(1,n)),  &
          grid2_temp(grid2_add_map1(n)) )
        endif
      end do

      do n=1,grid2_size
        if (grid2_array_norm(n) > zero) then
          print *, n, nint(grid2_array_norm(n)), grid2_temp(n)
        endif
      end do

      imin = minval(grid2_count, mask=(grid2_count > 0))
      imax = maxval(grid2_count)
      idiff =  (imax - imin)/10 + 1
      write(*,format_str) ' total number of dest cells            ',  &
                                                          grid2_size
      write(*,format_str) ' number of cells participating in remap',  &
                                            count(grid2_count > zero)
      write(*,format_str) ' min no of entries/row                 ',imin
      write(*,format_str) ' max no of entries/row                 ',imax
      print *,''

      
      format_str = "(A,I6,A,I6,A,I6)"

      imax = imin + idiff
      do n=1,10
        write(*,format_str) ' num of rows with entries between ', &
              imin,' - ',imax-1, ' :', &
              count(grid2_count >= imin .and. grid2_count < imax)
        imin = imin + idiff
        imax = imax + idiff
      end do

!-----------------------------------------------------------------------

      end program scrip_regrid

!***********************************************************************

      subroutine SCRIP_RegridExit(errorCode)

      use SCRIP_KindsMod   ! for SCRIP data types
      use SCRIP_CommMod    ! for SCRIP parallel environment
      use SCRIP_ErrorMod   ! SCRIP error checking and logging

      integer (SCRIP_i4), intent(in) :: errorCode

      call SCRIP_ErrorPrint(errorCode, masterTask)
      call SCRIP_CommExitMessageEnvironment

      stop

      end subroutine SCRIP_RegridExit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
