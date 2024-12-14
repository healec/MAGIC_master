! =========================================================
      subroutine out3(meqn,mbc,mx,my,&
     &     mz,xlower,ylower,zlower,dx,dy,dz,q,t,iframe)
! =========================================================

! Routine to save compressed chunked output from MAGIC
! Note current filter is based on GZIP compression
! Old indexing for output is preserved

     use mpi
     use hdf5
     implicit double precision (a-h,o-z)
    
     !------------------ HDF variables ------------------!
     integer(hid_t) :: plist_id      ! property list identifier
     integer(hid_t) :: dcpl          ! property list identifier
     integer(hid_t) :: file_id       ! file identifier
     integer(hid_t) :: dataset_id    ! dataset identifier
     integer(hid_t) :: dataspace_id  ! dataspace identifier
     double precision, dimension(25) :: attr_datacur
     INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/25/) ! Attribute dimension
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
     INTEGER(HID_T) :: attr_id ! Attribute identifier
     INTEGER(HID_T) :: aspace_id ! Attribute dataspace identifier
     INTEGER(HID_T) :: atype_id ! Attribute dataspace identifier
     INTEGER :: arank = 1 ! Attribute rank
     integer(hid_t) :: filespace, memspace, memd
     !---------------------------------------------------! 
	 
      parameter   (nDim = 3)

     !------------------ Filter variables ---------------!

     integer(hsize_t), dimension(4) :: cdims = (/1,1,1,1/) ! chunks data dimensions
     !---------------------------------------------------!
     
     !------------------ miscellaneous ------------------!
     character(len=3) :: c                                     ! dataset name for specific rank
     character(len=10) :: dataset_name
     integer :: rank = 4                                       ! data rank. q is 4D
     character(mpi_max_processor_name) hostname
     dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,&
	 &               1-mbc:mz+mbc), mtotal(nDim)
     integer(hsize_t), dimension(4) :: dimsf ! data dataset dimensions
     integer :: i,j,k,l,m,info,idd
	 real*8 :: ngrids_out
     character*20 fname
     common /mpicomm/ mpi_comm_3d, lx, ly, lz, mtotal, mstart
     common /mpi_proc_info/ np, id
     !---------------------------------------------------!

     ! initialize HDF5 fortran interface
    call h5open_f(ierr)
    ngrids_out = 1.d0
    
    ! define size of q for every core
    dimsf(1) = meqn
    dimsf(2) = mx
    dimsf(3) = my
    dimsf(4) = mz
    
!    allocate (data(dimsf(1),dimsf(2),dimsf(3),dimsf(4)))
      
     info = mpi_info_null
     
     fname = 'fort.qq' &
     & // char(ichar('0') + mod(iframe/1000,10)) &
     & // char(ichar('0') + mod(iframe/100,10)) &
     & // char(ichar('0') + mod(iframe/10,10)) &
     & // char(ichar('0') + mod(iframe,10)) &
     & // 'id' &
     & // char(ichar('0') + mod(id/1000,10)) &
     & // char(ichar('0') + mod(id/100,10)) &
     & // char(ichar('0') + mod(id/10,10)) &
     & // char(ichar('0') + mod(id,10)) &
     & // '.h5'
    
    
     ! create datatype for the attribute
 	 ! Copy existing datatype
 	 ! H5T_NATIVE_CHARACTER - copy this datatype
     ! atype_id - copy datatype to this variable
	call h5tcopy_f(h5t_native_double,atype_id,ierr)
    
    ! Current version sets chunks size as whole dataset of q
    cdims(1) = 1 !dimsf(1)
    cdims(2) = dimsf(2)
    cdims(3) = dimsf(3)
	cdims(4) = dimsf(4)
    
    	 ! create scalar dataspace for the attribute
    	 call h5screate_simple_f(arank,adims,aspace_id, ierr)
	
         ! create the hdf5 file
         ! filename: name of current file
         ! h5f_acc_trunc_f: rewrite if file exists
         ! file_id: identifier for file
         call h5fcreate_f(fname, h5f_acc_trunc_f, file_id, ierr)
     
         ! create the dataspace for the dataset
         ! rank: rank of datasets
         ! dimsf: size of every dimension
         ! dataspace_id: identifier of dataspace
         call h5screate_simple_f(rank, dimsf, dataspace_id, ierr)
         
         ! create properties variable for the data
         ! h5p_dataset_create_f: property to create a data dataset
         ! dcpl: variable for data dataset properties
         call h5pcreate_f(h5p_dataset_create_f, dcpl, ierr)
        
         ! attribute the chunk size
         ! dcpl: link this property variable with chunk size
         ! 4: dimension of q
         ! cdims: dimensions of chunk in every direction
         call h5pset_chunk_f(dcpl, 4, cdims, ierr)
         
         ! attribute the compression type (GZIP compression)
         ! dcpl: link this property variable with filter
         ! 6: compression rank
         call h5pset_deflate_f(dcpl, 6, ierr)

         ! attribute time of allocation of space for data in datasets
         ! h5d_alloc_time_early_f - allocate all space when the dataset is created
         call h5pset_alloc_time_f(dcpl, h5d_alloc_time_early_f, ierr)
        
! create name for every dataset

         dataset_name = "Pid"

         ! create dataset for this processor (based on id)
         ! file_id: name of file where to create dataset
         ! dataset_name: name of the dataset
         ! h5t_native_integer: type of data in dataset
         ! dataspace_id: dataspace used for dataset
         ! dataset_id: identifier for dataset
         ! dcpl_id: dataset creation property list
         call h5dcreate_f(file_id, dataset_name, h5t_native_double, &
                             dataspace_id, dataset_id, ierr, dcpl_id=dcpl)

        ! Attributes list
	    attr_datacur(1) = ngrids_out
        attr_datacur(2) = nDim
        attr_datacur(3) = t
        attr_datacur(4) = meqn
        attr_datacur(5) = 1.
        attr_datacur(6) = mtotal(1) ! Full mx (all processors)
        attr_datacur(7) = mtotal(2) ! Full my
        attr_datacur(8) = mtotal(3) ! Full mz
        attr_datacur(9) = 0.
        attr_datacur(10) = xlower ! lower boundary x
        attr_datacur(11) = ylower ! lower boundary y
        attr_datacur(12) = zlower ! lower boundary z 
        attr_datacur(13) = 0
        attr_datacur(14) = xlower+mtotal(1)*dx ! upper boundary x
        attr_datacur(15) = ylower+mtotal(2)*dy ! upper boundary y
        attr_datacur(16) = zlower+mtotal(3)*dz ! upper boundary z
        attr_datacur(17) = 0.
        attr_datacur(18) = dx
        attr_datacur(19) = dy
        attr_datacur(20) = dz
        attr_datacur(21) = 0.
        ! Next variable are not included for HDF4 version but needed for slicing in HDF5
        attr_datacur(22) = lx
        attr_datacur(23) = ly
        attr_datacur(24) = lz
        attr_datacur(25) = iframe ! Control variable
        
         call h5acreate_f(dataset_id,"Parameters",atype_id,aspace_id,attr_id, ierr)
         
         data_dims(1) = 25
         call h5awrite_f(attr_id, atype_id, attr_datacur, data_dims, ierr)
     	 
         ! close attribute
         call h5aclose_f(attr_id, ierr)    
         
         ! close access to the dataspace for attribute
         call h5sclose_f(aspace_id, ierr)    
         
         call h5tclose_f(atype_id, ierr)  


      ! setup file access property variable with parallel i/o access
      ! plist_id: property variable
      ! comm - mpi communicator for mpi/io
      ! info - info regarding file access patterns and file system specifications
     call h5pcreate_f(h5p_file_access_f, plist_id, ierr)
     call h5pset_fapl_mpio_f(plist_id, mpi_comm_3d, info, ierr)
     
     ! close the property list
     call h5pclose_f(plist_id, ierr) 
     

     call h5dwrite_f(dataset_id, h5t_native_double, q(1:10,&
     & 1:mx, 1:my, 1:mz),& 
     & dimsf, ierr)
	 
     call h5dclose_f(dataset_id,ierr)
     
     
     call h5fclose_f(file_id, ierr)
	
     ! close fortran interface
     call h5close_f(ierr)
     
     return
     end
