! The parallel exact percolation analysis (exact_percolation) against the serial one
! (clusteranalysis_exact, then span or span_simple) on random periodic lattices, on any
! number of threads and MPI ranks. Every rank makes the same lattices, analyses its own
! slab of each, and compares the slab with the same planes of the serial answer.
!
!   gfortran -O2 -fopenmp -cpp src/parallel.f90 src/percolation.f90 tests/test_percolation_mpi.f90 -o t
!   mpif90 -O2 -fopenmp -cpp -DUSE_MPI src/parallel.f90 src/percolation.f90 tests/test_percolation_mpi.f90 -o t
!   OMP_NUM_THREADS=3 mpirun -n 4 ./t
program test_percolation_mpi
    use parallel
    use percolation
    implicit none
    integer :: failures, trial, LX, LY, LZ, z0, z1, nc, nn_ref, nn, span_ref, span_par, i, nl
    integer*1, allocatable :: grid(:,:,:), ref(:,:,:), slab(:,:,:)
    integer, allocatable :: cluster(:,:,:), cl(:), sites_ref(:), sites(:), summary(:,:)
    real, allocatable :: r(:,:,:)
    real :: p
    integer :: seed(64)
    logical :: ok

    call par_init()
    failures = 0
    percolation_labelling = 1
    allocate(cl(5000000), summary(20, 2))
    seed = 4242
    call random_seed(put=seed(1:size_of_seed()))
    do trial = 1, 300
        LX = 3 + mod(trial * 7, 19)
        LY = 3 + mod(trial * 11, 17)
        LZ = max(par_size, 2 + mod(trial * 5, 23))
        p = 0.15 + 0.7 * mod(trial * 13, 100) / 100.0
        allocate(grid(LX, LY, LZ), ref(LX, LY, LZ), r(LX, LY, LZ), cluster(LX, LY, LZ))
        allocate(sites_ref(LX*LY*LZ))
        call random_number(r)
        grid = 0
        where (r < p) grid = 1
        call par_slab(LZ, par_rank, z0, z1)
        nl = LX*LY*(z1 - z0 + 1)
        allocate(slab(LX, LY, z1 - z0 + 1), sites(nl))

        ! span: every spanning cluster's sites, and their list
        ref = grid
        cluster = 0; cl = 0; nc = 0
        call clusteranalysis_exact(ref, cluster, cl, nc)
        nn_ref = -1
        call reference_span(ref, cluster, cl, nc, sites_ref, nn_ref)
        slab = grid(:, :, z0:z1)
        nn = -1
        call percolation_set_slab(LZ, z0)
        call exact_percolation(slab, summary, span_par, .True., sites, nn)
        ok = all(slab == ref(:, :, z0:z1))
        if (nn_ref >= 0) then
            ! the reference's sites in this slab, renumbered from the slab's first plane
            ok = ok .and. nn == count(sites_ref(1:nn_ref) > LX*LY*(z0-1) .and. sites_ref(1:nn_ref) <= LX*LY*z1)
            if (ok) ok = all(sites(1:nn) == pack(sites_ref(1:nn_ref), sites_ref(1:nn_ref) > LX*LY*(z0-1) .and. &
                                             sites_ref(1:nn_ref) <= LX*LY*z1) - LX*LY*(z0-1))
        else
            ok = ok .and. nn == -1
        end if
        if (.not. ok) then
            failures = failures + 1
            print '(a, i4, a, i3, a, 3i4, f6.2)', ' FAIL span: trial', trial, ' rank', par_rank, ' L', LX, LY, LZ, p
        end if

        ! span_simple: the first spanning cluster's directions
        ref = grid
        cluster = 0; cl = 0; nc = 0
        call clusteranalysis_exact(ref, cluster, cl, nc)
        span_ref = reference_simple(cluster, cl, nc)
        slab = grid(:, :, z0:z1)
        call exact_percolation(slab, summary, span_par, .False.)
        if (span_par /= span_ref) then
            failures = failures + 1
            print '(a, i4, a, i3, a, 2i3)', ' FAIL span_simple: trial', trial, ' rank', par_rank, ' spanning', &
                span_par, span_ref
        end if
        deallocate(grid, ref, r, cluster, sites_ref, slab, sites)
    end do
    failures = int(par_sum_i8(int(failures, 8)))
    if (par_rank == 0) then
        if (failures == 0) then
            print '(a, i3, a)', ' PASS parallel exact percolation matches the serial one on 300 lattices (', &
                par_size, ' ranks)'
        else
            print '(a, i6, a)', ' FAIL ', failures, ' mismatches'
        end if
    end if
    call par_finalize()
    if (failures /= 0) stop 1

contains

    integer function size_of_seed()
        call random_seed(size=size_of_seed)
    end function size_of_seed

    ! span, independently: clusters of at least min(LX, LY, LZ) sites that occupy every
    ! plane across x, y or z; the lattice keeps only their sites (and is left as it is
    ! when no cluster is that big)
    subroutine reference_span(g, cluster, cl, nc, sites, nn)
        integer*1, intent(inout) :: g(:,:,:)
        integer, intent(in) :: cluster(:,:,:), cl(:), nc
        integer, intent(out) :: sites(:)
        integer, intent(inout) :: nn
        logical :: keep(nc), any_cand
        integer :: n, i, j, k, s
        any_cand = .false.
        do n = 1, nc
            keep(n) = .false.
            if (cl(n) < min(size(g,1), size(g,2), size(g,3))) cycle
            any_cand = .true.
            keep(n) = directions(cluster, n) > 0
        end do
        if (.not. any_cand) return
        nn = 0
        s = 0
        do k = 1, size(g,3)
            do j = 1, size(g,2)
                do i = 1, size(g,1)
                    s = s + 1
                    g(i,j,k) = 0
                    if (cluster(i,j,k) > 0) then
                        if (keep(cluster(i,j,k))) then
                            g(i,j,k) = 1
                            nn = nn + 1
                            sites(nn) = s
                        end if
                    end if
                end do
            end do
        end do
    end subroutine reference_span

    integer function reference_simple(cluster, cl, nc)
        integer, intent(in) :: cluster(:,:,:), cl(:), nc
        integer :: n
        reference_simple = 0
        do n = 1, nc
            if (cl(n) < min(size(cluster,1), size(cluster,2), size(cluster,3))) cycle
            reference_simple = directions(cluster, n)
            if (reference_simple > 0) return
        end do
        reference_simple = 0
    end function reference_simple

    integer function directions(cluster, n)
        integer, intent(in) :: cluster(:,:,:), n
        integer :: i
        logical :: x, y, z
        x = .true.; y = .true.; z = .true.
        do i = 1, size(cluster, 1)
            if (.not. any(cluster(i,:,:) == n)) x = .false.
        end do
        do i = 1, size(cluster, 2)
            if (.not. any(cluster(:,i,:) == n)) y = .false.
        end do
        do i = 1, size(cluster, 3)
            if (.not. any(cluster(:,:,i) == n)) z = .false.
        end do
        directions = count((/ x, y, z /))
    end function directions

end program test_percolation_mpi
