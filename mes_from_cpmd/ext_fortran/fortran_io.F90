      subroutine fortran_read_unformatted(n_x,n_y,n_z,fn_inp,arr_out)
      implicit none

      integer, parameter :: dp = selected_real_kind(15, 307)

      integer, intent(in) :: n_x, n_y, n_z
      real(kind=dp),  intent(inout)  :: arr_out(n_x, n_y, n_z)
      character(len=255), intent(in) :: fn_inp

      integer ix, iy, iz, iostat
      iostat=0
      open(1,file=trim(fn_inp),status='old',form='unformatted',
     &     iostat=iostat)
      if(iostat.ne.0) then
        write(*,'(A,A)') 'Failed to open file ', trim(fn_inp)
        call exit(1)
      endif

      do iz=1,n_z
        do iy=1,n_y
          read(1,iostat=iostat) (arr_out(ix,iy,iz), ix=1,n_x)
        enddo
      enddo
      close(1)

      end subroutine fortran_read_unformatted

      subroutine fortran_write_unformatted(n_x,n_y,n_z,fn_out,arr_in)
      implicit none

      integer, parameter :: dp = selected_real_kind(15, 307)

      integer, intent(in) :: n_x, n_y, n_z
      real(kind=dp),  intent(inout)  :: arr_in(n_x, n_y, n_z)
      character(len=255), intent(in) :: fn_out

      integer ix, iy, iz, iostat
      iostat=0
      open(1,file=trim(fn_out),form='unformatted',iostat=iostat)
      if(iostat.ne.0) then
        write(*,'(A,A)') 'Failed to open file ', trim(fn_out)
        call exit(1)
      endif

      do iz=1,n_z
        do iy=1,n_y
          write(1,iostat=iostat) (arr_in(ix,iy,iz), ix=1,n_x)
        enddo
      enddo
      close(1)

      end subroutine fortran_write_unformatted
      
      subroutine fortran_read_pert(n_x,n_y,n_z,fn_inp,arr_out)
      implicit none

      integer, parameter :: dp = selected_real_kind(15, 307)

      integer, intent(in) :: n_x, n_y, n_z
      real(kind=dp),  intent(inout)  :: arr_out(n_x, n_y, n_z)
      character(len=255), intent(in) :: fn_inp

      integer ix, iy, iz, iostat
      iostat=0
      open(1,file=trim(fn_inp),status='old',form='unformatted',
     &     iostat=iostat)
      if(iostat.ne.0) then
        write(*,'(A,A)') 'Failed to open file ', trim(fn_inp)
        call exit(1)
      endif

      do ix=1,n_z
        do iy=1,n_y
          read(1,iostat=iostat) (arr_out(ix,iy,iz), iz=1,n_x)
        enddo
      enddo
      close(1)

      end subroutine fortran_read_pert

      subroutine fortran_write_pert(n_x,n_y,n_z,fn_out,arr_in)
      implicit none

      integer, parameter :: dp = selected_real_kind(15, 307)

      integer, intent(in) :: n_x, n_y, n_z
      real(kind=dp),  intent(inout)  :: arr_in(n_x, n_y, n_z)
      character(len=255), intent(in) :: fn_out

      integer ix, iy, iz, iostat
      iostat=0
      open(1,file=trim(fn_out),form='unformatted',iostat=iostat)
      if(iostat.ne.0) then
        write(*,'(A,A)') 'Failed to open file ', trim(fn_out)
        call exit(1)
      endif

      do ix=1,n_z
        do iy=1,n_y
          write(1,iostat=iostat) (arr_in(ix,iy,iz), iz=1,n_x)
        enddo
      enddo
      close(1)

      end subroutine fortran_write_pert
