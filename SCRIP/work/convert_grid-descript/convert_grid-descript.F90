!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this program is a short driver to read a variable from a netCDF
!     file with CF format and store it in a netCDF file with SCRIP format
!       -> needs interpolation weights file which stores grid
!          information in SCRIP format
!       -> grid of variable to read must correspond to destination(!)
!          grid of interpolation file
!       -> namelist to specify input/output/interp. file & field names:
!           convert_grid-descript_in
!       
!
!       based on:
!     CVS: $Id: scrip_test.f,v 1.6 2000/04/19 21:45:09 pwjones Exp $
!     
!       which comes with the following license:
!-----------------------------------------------------------------------
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

      program convert_grid_descript

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
              field_name_in,    &! variable name to be converted
              field_name_out,   &! name for converted variable in output_file
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
              nc_srcarray_id, nc_srcgradlat_id, nc_srcgradlon_id, &
              nc_dstarray1_id, nc_dstarray1a_id, nc_dstarray2_id, &
              nc_dsterror1_id, nc_dsterror1a_id, nc_dsterror2_id 

      !integer (SCRIP_i4), dimension(:), allocatable ::  &
      integer (SCRIP_i4)  ::  &
              nc_grid1size_id, nc_grid2size_id

!-----------------------------------------------------------------------

      character (SCRIP_charLength) ::  &
                dim_name    ! netCDF dimension name

      integer (SCRIP_i4) :: i,j,n,imin,imax,idiff,  &
          ip1,im1,jp1,jm1,nx,ny,  &  ! for computing bicub gradients
          in,is,ie,iw,ine,inw,ise,isw,  &
          iunit, &               ! unit number for input configuration file
          unit_attr_len          ! length of string which holds unit attribute
                                 ! (to be allocated)

      integer (SCRIP_i4), dimension(:), allocatable ::  &
          grid1_imask, grid2_imask, grid2_count

      real (SCRIP_r8) ::  &
          delew, delns,  &  ! variables for computing bicub gradients
          length            ! length scale for cosine hill test field

      real (SCRIP_r8), dimension(:), allocatable :: &
          grid1_array,  &
          grid1_tmp,  &
          grad1_lat,  &
          grad1_lon,  &
          grad1_latlon,  &
          grad1_lat_zero,  &
          grad1_lon_zero,  &
          grid2_array,  &
          grid2_err, &
          grid2_tmp

      !! temporary arrays for saving data to ncfiles (real & integer)
      !real (SCRIP_i4), dimension(:,:), allocatable ::  &
      !    grid1itmp2d, grid2itmp2d

      !real (SCRIP_r8), dimension(:,:), allocatable ::  &
      !    grid1tmp2d, grid2tmp2d

      ! used for uniform output format
      character (10) :: format_str 

      character(2000) :: title 
      character(len=:), allocatable :: unit_attr

      character (12) ::   &
          rtnName = 'convert_grid-descript'

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
        call SCRIP_ConvertExit(errorCode)
      endif

!!-----------------------------------------------------------------------
!!
!!     read namelist for file and mapping info
!!
!!-----------------------------------------------------------------------
!
!      call SCRIP_IOUnitsGet(iunit)
!      open(iunit, file='scrip_test_in', status='old', form='formatted')
!      read(iunit, nml=remap_inputs)
!      call SCRIP_IOUnitsRelease(iunit)
!      write(*,nml=remap_inputs)

!-----------------------------------------------------------------------
!
!  read input namelist   
!
!-----------------------------------------------------------------------

      call SCRIP_ConfigOpen(iunit, errorCode, 'convert_grid-descript_in')
      if (SCRIP_ErrorCheck(errorCode, rtnName, &
          'error opening config file')) call SCRIP_ConvertExit(errorCode)

      call SCRIP_ConfigRead(iunit, 'convertInputs',                    &
            'input_file', input_file, 'unknown', errorCode,            &
             outStringBefore='file containing data for regridding: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_ConvertExit(errorCode)


      call SCRIP_ConfigRead(iunit, 'convertInputs',                    &
                'field_name_in', field_name_in, 'unknown', errorCode,  &
                outStringBefore= 'variable/field to be regridded: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_ConvertExit(errorCode)

      call SCRIP_ConfigRead(iunit, 'convertInputs',                    &
            'field_name_out', field_name_out, 'unknown', errorCode,    &
            outStringBefore=                                           &
            'interpolation file containing weights for regridding: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_ConvertExit(errorCode)

      call SCRIP_ConfigRead(iunit, 'convertInputs',                    &
                    'interp_file', interp_file, 'unknown', errorCode,  &
                    outStringBefore=                                   &
            'interpolation file containing weights for regridding: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_ConvertExit(errorCode)

      call SCRIP_ConfigRead(iunit, 'convertInputs',                    &
                    'output_file', output_file, 'unknown', errorCode,  &
                    outStringBefore='output file with regridded data: ')
      if (SCRIP_ErrorCheck(errorCode, rtnName,   &
          'error reading input_file')) call SCRIP_ConvertExit(errorCode)

      call SCRIP_ConfigClose(iunit, errorCode)
      if (SCRIP_ErrorCheck(errorCode, rtnName, &
          'error closing config file')) call SCRIP_ConvertExit(errorCode)

!-----------------------------------------------------------------------
!
!     read remapping data
!
!-----------------------------------------------------------------------

      ! reading grid specification data from interpolation file and 
      !   saving 'title' metadata to map_name
      call read_remap(map_name, interp_file, errorCode)
      if (SCRIP_ErrorCheck(errorCode, rtnName,  &
          'error reading remap file')) call SCRIP_ConvertExit(errorCode)

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

      allocate (grid1_array    (grid1_size),  &
                grid1_tmp      (grid1_size),  &
                grad1_lat      (grid1_size),  &
                grad1_lon      (grid1_size),  &
                grad1_lat_zero (grid1_size),  &
                grad1_lon_zero (grid1_size),  &
                grid1_imask    (grid1_size),  &
                grid2_array    (grid2_size),  &
                grid2_err      (grid2_size),  &
                grid2_tmp      (grid2_size),  &
                grid2_imask    (grid2_size),  &
                grid2_count    (grid2_size))

      !where (grid1_mask)
      !  grid1_imask = 1
      !elsewhere
      !  grid1_imask = 0
      !endwhere
      where (grid2_mask)
        grid2_imask = 1
      elsewhere
        grid2_imask = 0
      endwhere

!-----------------------------------------------------------------------
!
!    reading field_name_in variable from input_file 
!       -> storing in grid2_array
!
!-----------------------------------------------------------------------

      ncstat = nf90_open(input_file, NF90_NOWRITE, nc_inputfile_id) 
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,        &
                                 'error opening input_file'))       &
          call SCRIP_ConvertExit(errorCode)

      ncstat = nf90_inq_varid(nc_inputfile_id, field_name_in,       &
                                nc_fieldname_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,        &
                                'error getting field_name_in id'))  &
          call SCRIP_ConvertExit(errorCode)

      ncstat = nf90_get_var(nc_inputfile_id, nc_fieldname_id,       &
                        grid2_array, start = (/1, 1, 1, 1/),        &
                        count = (/grid2_dims(1), grid2_dims(2), 1, 1/) )
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,        &
                                 'error reading field_name_in'))    &
          call SCRIP_ConvertExit(errorCode)

      ncstat = nf90_inquire_attribute(nc_inputfile_id,              &
                        nc_fieldname_id, "units", len=unit_attr_len)
                      
      allocate( character( len=unit_attr_len ) :: unit_attr) 

      ncstat = nf90_get_att(nc_inputfile_id, nc_fieldname_id,      &
                                "units", unit_attr)
      

      ncstat = nf90_close(nc_inputfile_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,        &
                                 'error closing input_file'))       &
          call SCRIP_ConvertExit(errorCode)

!-----------------------------------------------------------------------
!
!     convert grid area values to unit steradians
!
!-----------------------------------------------------------------------

      !!! MOM & LAND radius are defined in:
      !     mom5.0.2/src/shared/mosaic/constant.h
      !!! PISM radius calculated by:
      !     PISM_WGS84_radius.py

        ! to do: omit missing values in land mask (-999) from conversion


      if( INDEX(output_file, "MOM") ) then
          if( unit_attr .EQ. "m^2") then
              grid2_array = grid2_array / (6371 * 1e3)**2  
              unit_attr = "square radians"
          else
              call SCRIP_ErrorSet(errorCode, rtnName,    &
                  'error reading unit of area variable')
              call SCRIP_ConvertExit(errorCode)
          endif
      else if( INDEX(output_file, "LAND") ) then
          if( unit_attr .EQ. "m2") then
              grid2_array = grid2_array / (6371 * 1e3)**2  
              unit_attr = "square radians"
          else
              call SCRIP_ErrorSet(errorCode, rtnName,    &
                  'error reading unit of area variable')
              call SCRIP_ConvertExit(errorCode)
          endif
      else if( INDEX(output_file, "PISM") ) then
          if( unit_attr .EQ. "km2") then
              grid2_array = grid2_array / ( 6359.0360484866882 )**2  
              unit_attr = "square radians"
          else
              call SCRIP_ErrorSet(errorCode, rtnName,    &
                  'error reading unit of area variable')    
              call SCRIP_ConvertExit(errorCode)
          endif
      else 
          call SCRIP_ErrorSet(errorCode, rtnName,    &
              'error determing grid time (none of MOM/PISM/Land)') 
          call SCRIP_ConvertExit(errorCode)
      endif




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
          'error opening output file')) call SCRIP_ConvertExit(errorCode)

      WRITE(title,*) 'variable >>', trim(field_name_in),               &
                 '<< from file ', trim(input_file),                    &
                 ' (CF format) defined on destination-grid of & 
                    &weighting file ',                                 &
                 trim(interp_file), ' with name: >>', trim(map_name),  &
                 '<< saved in variable >>',            &
                 trim(field_name_out), '<< (SCRIP format)'

      ncstat = nf90_put_att(nc_outfile_id, NF90_GLOBAL, 'title',     &
                            !map_name)
                            title) 
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,         &
                                 'error writing output file title')) &
          call SCRIP_ConvertExit(errorCode)

      !***
      !*** define grid size dimensions
      !***

      format_str = '(a9,i1)'

        ncstat = nf90_def_dim(nc_outfile_id, 'grid_size', &
                              grid2_size, nc_grid2size_id)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
                'error defining grid_size dimension'))          &
            call SCRIP_ConvertExit(errorCode)


      !ncstat = nf90_def_var(nc_outfile_id, 'dst_grid_center_lat', &
      ncstat = nf90_def_var(nc_outfile_id, 'grid_center_lat', &
                            NF90_DOUBLE, nc_grid2size_id, &
                            nc_dstgrdcntrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
          'error defining grid center lat'))  &
          call SCRIP_ConvertExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_dstgrdcntrlat_id,  &
                            'units', 'radians')  
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                    'error adding units grid center lon'))  &
          call SCRIP_ConvertExit(errorCode)

      !***
      !*** define grid center longitude array
      !***

      !ncstat = nf90_def_var(nc_outfile_id, 'dst_grid_center_lon',  &
      ncstat = nf90_def_var(nc_outfile_id, 'grid_center_lon',  &
                            NF90_DOUBLE, nc_grid2size_id,  &
                            nc_dstgrdcntrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                        'error defining grid center lon'))  &
          call SCRIP_ConvertExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_dstgrdcntrlon_id, &
                            'units', 'radians')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                    'error adding units grid center lon'))  &
          call SCRIP_ConvertExit(errorCode)

      !***
      !*** define grid mask
      !***

      !ncstat = nf90_def_var(nc_outfile_id, 'dst_grid_imask', NF90_INT, &
      ncstat = nf90_def_var(nc_outfile_id, 'grid_imask', NF90_INT, &
                            nc_grid2size_id, nc_dstgrdimask_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error defining dst mask'))  &
          call SCRIP_ConvertExit(errorCode)

      ncstat = nf90_put_att(nc_outfile_id, nc_dstgrdimask_id, &
                            'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error adding units to dst mask'))  &
          call SCRIP_ConvertExit(errorCode)


      !***
      !*** define destination arrays
      !***

      !ncstat = nf90_def_var(nc_outfile_id, 'dst_array1',    &
      ncstat = nf90_def_var(nc_outfile_id, field_name_out,  &
                            NF90_DOUBLE, nc_grid2size_id,   &
                            nc_dstarray1_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                       'error defining output field'))  &
          call SCRIP_ConvertExit(errorCode)


      ncstat = nf90_put_att(nc_outfile_id, nc_dstarray1_id, &
                            'units', unit_attr)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error adding units to dst mask'))  &
          call SCRIP_ConvertExit(errorCode)


      !***
      !*** end definition stage
      !***

      ncstat = nf90_enddef(nc_outfile_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                       'error ending netCDF definition phase'))  &
          call SCRIP_ConvertExit(errorCode)

!-----------------------------------------------------------------------
!
!     write some grid info
!
!-----------------------------------------------------------------------

      !***
      !*** write grid center latitude array
      !***

      !!! write one dimensional center latitude array
      ncstat = nf90_put_var(nc_outfile_id, nc_dstgrdcntrlat_id, &
                            grid2_center_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                       'error writing dst grid center lat'))  &
          call SCRIP_ConvertExit(errorCode)

      !***
      !*** write grid center longitude array
      !***

      !!! write one dimensional center latitude array
      ncstat = nf90_put_var(nc_outfile_id, nc_dstgrdcntrlon_id, &
                            grid2_center_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing dst grid lon'))  &
          call SCRIP_ConvertExit(errorCode)

      !***
      !*** write grid mask
      !***
      
      !!! write one dimensional grid mask 
      ncstat = nf90_put_var(nc_outfile_id, nc_dstgrdimask_id, &
                            grid2_imask)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                            'error writing dst grid mask'))  &
          call SCRIP_ConvertExit(errorCode)




!-----------------------------------------------------------------------
!
!     write results to NetCDF file
!
!-----------------------------------------------------------------------

      !!! write one dimensional grid2 array 
      ncstat = nf90_put_var(nc_outfile_id, nc_dstarray1_id, &
                            grid2_array  )
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,  &
                            'error writing remapped field'))  &
          call SCRIP_ConvertExit(errorCode)


!-----------------------------------------------------------------------
!
!     close netCDF file
!
!-----------------------------------------------------------------------

      ncstat = nf90_close(nc_outfile_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, &
                                 'error closing output file'))  &
          call SCRIP_ConvertExit(errorCode)


!-----------------------------------------------------------------------

      end program convert_grid_descript

!***********************************************************************

      subroutine SCRIP_ConvertExit(errorCode)

      use SCRIP_KindsMod   ! for SCRIP data types
      use SCRIP_CommMod    ! for SCRIP parallel environment
      use SCRIP_ErrorMod   ! SCRIP error checking and logging

      integer (SCRIP_i4), intent(in) :: errorCode

      call SCRIP_ErrorPrint(errorCode, masterTask)
      call SCRIP_CommExitMessageEnvironment

      stop

      end subroutine SCRIP_ConvertExit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
