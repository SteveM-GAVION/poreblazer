!--------------------------------------------------------------------------------------
! Distribution of the grid over MPI ranks (fork only).
!
! The grid is split into slabs of whole z planes, one per rank: rank r holds planes
! z0..z1. The atoms are not split; every rank reads the whole structure, so no step
! needs the atoms of another rank. What crosses ranks is small: one plane of cluster
! labels between neighbouring slabs (the halo of the percolation analysis), the labels
! of the clusters that touch a slab boundary, and sums, minima and maxima.
!
! Built without USE_MPI (the default Makefile target) every routine here is a stub for
! a single rank, and the program is the serial one. Built with -DUSE_MPI (make mpi), it
! runs under mpirun or srun; on one rank it does exactly what the serial build does.
!--------------------------------------------------------------------------------------

module parallel
#ifdef USE_MPI
    use mpi
#endif
    implicit none
    private

    integer, public :: par_rank = 0, par_size = 1

    public :: par_init, par_finalize, par_abort, par_slab, par_owner
    public :: par_sum_i8, par_sum_i4_array, par_sum_i8_array, par_max_d, par_min_d, par_or_i1_array
    public :: par_exscan_i8, par_allgather_i8, par_shift_down_i8
    public :: par_send_d, par_recv_d, par_send_r4, par_recv_r4, par_send_i8, par_recv_i8
    public :: par_allgather_d, par_barrier, par_max_d_array, par_bcast_line

contains

    subroutine par_init()
#ifdef USE_MPI
        integer :: ierr, provided
        call mpi_init_thread(MPI_THREAD_FUNNELED, provided, ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, par_rank, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, par_size, ierr)
#endif
    end subroutine par_init

    subroutine par_finalize()
#ifdef USE_MPI
        integer :: ierr
        call mpi_finalize(ierr)
#endif
    end subroutine par_finalize

    ! Stop every rank, with a message from the one that found the problem
    subroutine par_abort(message)
        character(len=*), intent(in) :: message
#ifdef USE_MPI
        integer :: ierr
#endif
        write(0,*) trim(message)
#ifdef USE_MPI
        call mpi_abort(MPI_COMM_WORLD, 1, ierr)
#endif
        stop 1
    end subroutine par_abort

    subroutine par_barrier()
#ifdef USE_MPI
        integer :: ierr
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_barrier

    ! The planes z0..z1 (of nz) that rank r holds: as even a split as whole planes allow
    subroutine par_slab(nz, r, z0, z1)
        integer, intent(in)  :: nz, r
        integer, intent(out) :: z0, z1
        z0 = int((int(r, 8)*nz)/par_size) + 1
        z1 = int((int(r + 1, 8)*nz)/par_size)
    end subroutine par_slab

    ! The rank that holds plane z (of nz)
    integer function par_owner(nz, z)
        integer, intent(in) :: nz, z
        integer :: r, z0, z1
        r = int((int(z - 1, 8)*par_size)/nz)
        r = max(0, min(par_size - 1, r))
        do
            call par_slab(nz, r, z0, z1)
            if(z < z0) then
                r = r - 1
            else if(z > z1) then
                r = r + 1
            else
                exit
            end if
        end do
        par_owner = r
    end function par_owner

    integer*8 function par_sum_i8(x)
        integer*8, intent(in) :: x
#ifdef USE_MPI
        integer :: ierr
        call mpi_allreduce(x, par_sum_i8, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        par_sum_i8 = x
#endif
    end function par_sum_i8

    subroutine par_sum_i4_array(a, n)
        integer, intent(inout) :: a(*)
        integer, intent(in)    :: n
#ifdef USE_MPI
        integer :: ierr
        if(par_size > 1 .and. n > 0) &
            call mpi_allreduce(MPI_IN_PLACE, a, n, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_sum_i4_array

    subroutine par_sum_i8_array(a)
        integer*8, intent(inout) :: a(:)
#ifdef USE_MPI
        integer :: ierr
        if(par_size > 1 .and. size(a) > 0) &
            call mpi_allreduce(MPI_IN_PLACE, a, size(a), MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_sum_i8_array

    ! Element-wise OR of 0/1 bytes over the ranks
    subroutine par_or_i1_array(a, n)
        integer*1, intent(inout) :: a(*)
        integer, intent(in)      :: n
#ifdef USE_MPI
        integer :: ierr
        if(par_size > 1 .and. n > 0) &
            call mpi_allreduce(MPI_IN_PLACE, a, n, MPI_BYTE, MPI_BOR, MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_or_i1_array

    real*8 function par_max_d(x)
        real*8, intent(in) :: x
#ifdef USE_MPI
        integer :: ierr
        call mpi_allreduce(x, par_max_d, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
        par_max_d = x
#endif
    end function par_max_d

    real*8 function par_min_d(x)
        real*8, intent(in) :: x
#ifdef USE_MPI
        integer :: ierr
        call mpi_allreduce(x, par_min_d, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
        par_min_d = x
#endif
    end function par_min_d

    ! The sum of x over the ranks below this one (0 on rank 0)
    integer*8 function par_exscan_i8(x)
        integer*8, intent(in) :: x
#ifdef USE_MPI
        integer*8 :: all(0:par_size-1)
        integer :: ierr
        call mpi_allgather(x, 1, MPI_INTEGER8, all, 1, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
        par_exscan_i8 = 0
        if(par_rank > 0) par_exscan_i8 = sum(all(0:par_rank-1))
#else
        par_exscan_i8 = 0
#endif
    end function par_exscan_i8

    ! Every rank's a, concatenated in rank order, into out (allocated here)
    subroutine par_allgather_i8(a, out)
        integer*8, intent(in)               :: a(:)
        integer*8, allocatable, intent(out) :: out(:)
#ifdef USE_MPI
        integer :: counts(0:par_size-1), displs(0:par_size-1), ierr, r
        call mpi_allgather(size(a), 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        displs(0) = 0
        do r = 1, par_size - 1
            displs(r) = displs(r-1) + counts(r-1)
        end do
        allocate(out(max(1, sum(counts))))
        call mpi_allgatherv(a, size(a), MPI_INTEGER8, out, counts, displs, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
        if(sum(counts) == 0) then
            deallocate(out)
            allocate(out(0))
        end if
#else
        allocate(out(size(a)))
        out = a
#endif
    end subroutine par_allgather_i8

    subroutine par_allgather_d(a, out)
        real*8, intent(in)               :: a(:)
        real*8, allocatable, intent(out) :: out(:)
#ifdef USE_MPI
        integer :: counts(0:par_size-1), displs(0:par_size-1), ierr, r
        call mpi_allgather(size(a), 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        displs(0) = 0
        do r = 1, par_size - 1
            displs(r) = displs(r-1) + counts(r-1)
        end do
        allocate(out(max(1, sum(counts))))
        call mpi_allgatherv(a, size(a), MPI_DOUBLE_PRECISION, out, counts, displs, MPI_DOUBLE_PRECISION, &
                            MPI_COMM_WORLD, ierr)
        if(sum(counts) == 0) then
            deallocate(out)
            allocate(out(0))
        end if
#else
        allocate(out(size(a)))
        out = a
#endif
    end subroutine par_allgather_d

    ! The halo: send a to the rank below (periodic) and receive the rank above's into b.
    ! On one rank b is a.
    subroutine par_shift_down_i8(a, b)
        integer*8, intent(in)  :: a(:)
        integer*8, intent(out) :: b(:)
#ifdef USE_MPI
        integer :: ierr, below, above
        below = modulo(par_rank - 1, par_size)
        above = modulo(par_rank + 1, par_size)
        call mpi_sendrecv(a, size(a), MPI_INTEGER8, below, 11, b, size(b), MPI_INTEGER8, above, 11, &
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#else
        b = a
#endif
    end subroutine par_shift_down_i8


    ! Point-to-point transfers of n values, for pipelines and for gathering output on rank 0
    subroutine par_send_d(x, n, dest, tag)
        real*8, intent(in)  :: x(*)
        integer, intent(in) :: n, dest, tag
#ifdef USE_MPI
        integer :: ierr
        call mpi_send(x, n, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_send_d

    subroutine par_recv_d(x, n, source, tag)
        real*8, intent(inout) :: x(*)
        integer, intent(in)   :: n, source, tag
#ifdef USE_MPI
        integer :: ierr
        call mpi_recv(x, n, MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
    end subroutine par_recv_d

    subroutine par_send_r4(x, n, dest, tag)
        real, intent(in)    :: x(*)
        integer, intent(in) :: n, dest, tag
#ifdef USE_MPI
        integer :: ierr
        call mpi_send(x, n, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_send_r4

    subroutine par_recv_r4(x, n, source, tag)
        real, intent(inout) :: x(*)
        integer, intent(in) :: n, source, tag
#ifdef USE_MPI
        integer :: ierr
        call mpi_recv(x, n, MPI_REAL, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
    end subroutine par_recv_r4

    subroutine par_send_i8(x, n, dest, tag)
        integer*8, intent(in) :: x(*)
        integer, intent(in)   :: n, dest, tag
#ifdef USE_MPI
        integer :: ierr
        call mpi_send(x, n, MPI_INTEGER8, dest, tag, MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_send_i8

    subroutine par_recv_i8(x, n, source, tag)
        integer*8, intent(inout) :: x(*)
        integer, intent(in)      :: n, source, tag
#ifdef USE_MPI
        integer :: ierr
        call mpi_recv(x, n, MPI_INTEGER8, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
    end subroutine par_recv_i8

    ! Rank 0's line, on every rank (only rank 0 can read stdin under mpirun)
    subroutine par_bcast_line(line)
        character(len=*), intent(inout) :: line
#ifdef USE_MPI
        integer :: ierr
        if(par_size > 1) call mpi_bcast(line, len(line), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_bcast_line

    ! Element-wise maximum over the ranks
    subroutine par_max_d_array(a, n)
        real*8, intent(inout) :: a(*)
        integer, intent(in)   :: n
#ifdef USE_MPI
        integer :: ierr
        if(par_size > 1 .and. n > 0) &
            call mpi_allreduce(MPI_IN_PLACE, a, n, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
    end subroutine par_max_d_array

end module parallel
