!----------------------------------------------------------
! This module performs percolation analisys of lattice: returns number and sizes of lattice clusters
! and wheather the system is percolated
!----------------------------------------------------------

Module percolation

    Use parallel, only: par_size, par_rank, par_abort, par_shift_down_i8, par_allgather_i8, &
                        par_sum_i8, par_sum_i8_array, par_or_i1_array
    !$ Use omp_lib
    Implicit None
    Save

    Private
    Public :: percolation_calc, percolation_calc_simple, percolation_set_slab
    Public :: clusteranalysis, clusteranalysis_exact     ! for tests/test_percolation.f90
    Public :: exact_percolation                          ! for tests/test_percolation_mpi.f90

    ! Cluster labelling: 0 is Poreblazer 3.0.5's (the default, so results match upstream),
    ! 1 is exact union-find. Poreblazer's labelling can split one connected cluster into
    ! several; see FORK.md.
    Integer, Public :: percolation_labelling = 0

    ! The lattices given to the percolation analysis are this rank's slab of the grid:
    ! planes perc_z0 .. perc_z0 + size(lattice, 3) - 1 of perc_nz. One rank holds the
    ! whole grid (perc_nz = 0 means the lattice given).
    Integer :: perc_nz = 0, perc_z0 = 1


Contains

    Subroutine percolation_set_slab(nz, z0)
        Integer, Intent(In) :: nz, z0
        perc_nz = nz
        perc_z0 = z0
    End Subroutine percolation_set_slab

!-----------------------------------------------------------------------------------------------
!!subroutine that initiates percolation analysis
!-----------------------------------------------------------------------------------------------

    Subroutine percolation_calc(lattice_in, n_sites, nn_sites, cl_summary)

        Integer*1, Dimension(:,:,:), Intent(InOut)                                     :: lattice_in
        Integer, Dimension(:), Intent(InOut), Optional                               :: n_sites
        Integer, Intent(InOut), Optional                                             :: nn_sites
        Integer, Dimension(:,:), Intent(InOut)                                       :: cl_summary
        Integer, Dimension(:,:,:), allocatable                                       :: cluster
        Integer, Dimension(:), allocatable                                           :: cl
        Integer, Dimension(:), allocatable                                           :: trcl
        Integer                                                                      :: nc, spanning

        If(percolation_labelling == 1) Then
            Call exact_percolation(lattice_in, cl_summary, spanning, .True., n_sites, nn_sites)
            Return
        End If
        Call require_whole_grid(lattice_in)
        allocate(cluster(size(lattice_in, 1), size(lattice_in, 2), size(lattice_in, 3)))
        allocate(cl(5000000),trcl(5000000))
        nc = 0
        cluster = 0
        cl = 0
        trcl = 0

        Call clusteranalysis(lattice_in,cluster,cl,trcl,nc)
        Call span(lattice_in, n_sites, nn_sites, cluster, cl, nc, cl_summary)
        deallocate(cluster)
        deallocate(cl)
        deallocate(trcl)
    End Subroutine percolation_calc

!-----------------------------------------------------------------------------------------------
!!subroutine that initiates percolation analysis
!-----------------------------------------------------------------------------------------------

    Subroutine percolation_calc_simple(lattice_in, cl_summary, spanning)

        Integer*1, Dimension(:,:,:), Intent(InOut)                                       :: lattice_in
        Integer, Dimension(:,:), Intent(InOut)                                         :: cl_summary
        Integer, Dimension(:,:,:), allocatable                                         :: cluster
        Integer, Dimension(:), allocatable                                             :: cl,trcl
        Integer, Intent(InOut)                                                         :: spanning
        Integer                                                                        :: nc

        If(percolation_labelling == 1) Then
            Call exact_percolation(lattice_in, cl_summary, spanning, .False.)
            Return
        End If
        Call require_whole_grid(lattice_in)
        allocate(cluster(size(lattice_in, 1), size(lattice_in, 2), size(lattice_in, 3)))
        allocate(cl(100000),trcl(100000))
        nc = 0
        cluster = 0
        cl = 0
        trcl = 0
        spanning = 0

        Call clusteranalysis(lattice_in,cluster,cl,trcl,nc)

        Call span_simple(lattice_in, cluster, cl, nc, cl_summary, spanning)
        deallocate(cluster, cl, trcl)
    End Subroutine percolation_calc_simple



!---------------------------------------------------------------------
! Subroutine that takes a 3D grid of sites (0 or 1) and performs cluster
! analysis
!---------------------------------------------------------------------

    Subroutine clusteranalysis(ngrid,cluster,cl,trcl,nc)

        Integer*1, Dimension(:,:,:), Intent(InOut)     :: ngrid
        Integer, Dimension(:,:,:), Intent(InOut)  :: cluster
        Integer, Dimension(:), Intent(InOut)      :: cl,trcl
        Integer, Intent(InOut)                       :: nc
        Integer, Dimension(26)                       :: local
        Integer                                      :: i, j, k, l, LX, LY, LZ, scenario, ncold, trcli, icount, i1
        Logical :: loop
        Real :: xr, yr, zr
        Integer, Dimension(:), allocatable                                             :: newlabel


        ! newlabel(t) is the final label of the clusters whose root label is t (0 until the
        ! first site of that cluster is met). Upstream searched a list of the labels so far
        ! for every site, which is quadratic in the number of clusters and overflowed after
        ! 100000 clusters; the labels are the same.
        allocate(newlabel(size(trcl)))
        newlabel = 0

        LX = size(ngrid,1)
        LY = size(ngrid,2)
        LZ = size(ngrid,3)

        do k=1, LZ
            do j=1, LY
                do i=1, LX
                   If(ngrid(i,j,k) == 1) Then
                       If(cluster(i,j,k) /= 0) cycle

                       Call reveal_local(i,j,k,cluster,ngrid,local, scenario)

                       If(scenario == 1) Then
                           Call scenario_1(i,j,k, cluster, ngrid, cl, trcl, nc)
                       Else
                           Call scenario_2(i,j,k, cluster, ngrid, cl, trcl, nc)
                       End If

                    End If

                end do
            end do
        end do
     
        icount = 0 
        nc = 0

        do k=1, LZ
            do j=1, LY
                do i=1, LX

                        If(ngrid(i,j,k) == 0) cycle
                        icount =  icount + 1
                        trcli = trcl(cluster(i, j, k))
                        if(newlabel(trcli) == 0) then
                            nc = nc + 1
                            newlabel(trcli) = nc
                        end if
                        cluster(i, j, k) = newlabel(trcli)
                        cl(newlabel(trcli)) = cl(newlabel(trcli)) + 1
                end do
            end do
        end do
    
     do i=1, nc
     end do

    ! Several additional loops may be required to reconnect the neighbouring clusters

!        do l=1,100
!            ncold=nc
!            do k=1, LZ
!                do j=1, LY
!                    do i=1, LX
!                     If(ngrid(i,j,k) == 0) cycle
!                     write(20,*) i, j, k, cluster(i,j,k), cl(cluster(i,j,k)), nc
!
!                        If(ngrid(i,j,k) == 0) cycle
!                        Call update(i,j,k,cluster,ngrid,cl,nc)
!
!                    end do
!                end do
!            end do
!
!            If(nc==ncold) exit
!         end do
!     print*, "hello2"
     deallocate(newlabel)
     End Subroutine clusteranalysis

!---------------------------------------------------------------------
! Exact cluster labelling: union-find over the occupied sites with the same periodic
! six-neighbour connectivity, clusters numbered in order of their first site in scan
! order (the order clusteranalysis numbers them in). Two sites share a label exactly
! when a path of occupied neighbours joins them.
!---------------------------------------------------------------------

    Subroutine clusteranalysis_exact(ngrid, cluster, cl, nc)
        Integer*1, Dimension(:,:,:), Intent(In)   :: ngrid
        Integer, Dimension(:,:,:), Intent(InOut)  :: cluster
        Integer, Dimension(:), Intent(InOut)      :: cl
        Integer, Intent(InOut)                    :: nc
        Integer, Dimension(:), allocatable        :: parent
        Integer                                   :: i, j, k, LX, LY, LZ, s, r

        LX = size(ngrid,1)
        LY = size(ngrid,2)
        LZ = size(ngrid,3)
        allocate(parent(LX*LY*LZ))
        do s=1, LX*LY*LZ
            parent(s) = s
        end do
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    if(ngrid(i,j,k) /= 1) cycle
                    s = i + LX*(j-1) + LX*LY*(k-1)
                    if(ngrid(modulo(i,LX)+1,j,k) == 1) call unite(s, modulo(i,LX)+1 + LX*(j-1) + LX*LY*(k-1))
                    if(ngrid(i,modulo(j,LY)+1,k) == 1) call unite(s, i + LX*modulo(j,LY) + LX*LY*(k-1))
                    if(ngrid(i,j,modulo(k,LZ)+1) == 1) call unite(s, i + LX*(j-1) + LX*LY*modulo(k,LZ))
                end do
            end do
        end do
        ! Number the clusters in order of their first site; a numbered root holds -number
        cluster = 0
        nc = 0
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    if(ngrid(i,j,k) /= 1) cycle
                    r = find(i + LX*(j-1) + LX*LY*(k-1))
                    if(parent(r) > 0) then
                        nc = nc + 1
                        if(nc > size(cl)) stop "clusteranalysis_exact: more clusters than the cluster size array holds"
                        parent(r) = -nc
                    end if
                    cluster(i,j,k) = -parent(r)
                    cl(-parent(r)) = cl(-parent(r)) + 1
                end do
            end do
        end do
        deallocate(parent)

    contains

        ! The root of x, halving the path on the way (roots point to themselves, or hold
        ! -number once numbered)
        integer function find(x)
            integer, intent(in) :: x
            integer :: y
            y = x
            do while(parent(y) > 0 .and. parent(y) /= y)
                if(parent(parent(y)) > 0) parent(y) = parent(parent(y))
                y = parent(y)
            end do
            find = y
        end function find

        subroutine unite(a, b)
            integer, intent(in) :: a, b
            integer :: ra, rb
            ra = find(a)
            rb = find(b)
            if(ra /= rb) parent(max(ra, rb)) = min(ra, rb)
        end subroutine unite

    End Subroutine clusteranalysis_exact

!---------------------------------------------------------------------
! Subroutine which reveals status of the neighbouring sites
!---------------------------------------------------------------------

    Subroutine reveal_local(i0,j0,k0,cluster,ngrid,local, scenario)
        Integer, Intent(In)                      :: i0,j0,k0
        Integer, Dimension(:,:,:), Intent(In)     :: cluster
        Integer*1, Dimension(:,:,:), Intent(In)     ::ngrid
        Integer, Dimension(:), Intent(InOut)        :: local
        Integer, Intent(InOut)                     :: scenario
        Integer                                  :: i1, j1, k1, i, j, k, ic

    ! Default scenario 1: neighbouring sites are not occupied
    ! or occupied but not assigned to any clusters

     scenario = 1
     ic = 1

    ! Six neighbouring sites

     k1=k0-1; j1=j0; i1=i0    
     If(k1<1) k1=size(cluster,3)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then
           local(ic) = cluster(i1,j1,k1)

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters

           scenario = 2
           Return
           End If
        End If

     k1=k0+1; j1=j0; i1=i0
     If(k1>size(cluster,3)) k1= 1
        If(ngrid(i1,j1,k1)==1) Then                 
           If(cluster(i1,j1,k1)/=0) Then             
           local(ic) = cluster(i1,j1,k1)                 

           ! Scenario 2 where some of the sites are occupied and                 
           ! assigned to some clusters                 

           scenario = 2                 
           Return                 
           End If             
        End If

     k1=k0; j1=j0-1; i1=i0
     If(j1<1) j1=size(cluster,2)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then                 
           local(ic) = cluster(i1,j1,k1)                 

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters

           scenario = 2
           Return
           End If
        End If
      
     k1=k0; j1=j0+1; i1=i0
     If(j1>size(cluster,2)) j1= 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then                 
           local(ic) = cluster(i1,j1,k1)                 

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters

           scenario = 2
           Return
           End If
        End If

     k1=k0; j1=j0; i1=i0-1
     If(i1<1) i1=size(cluster,1)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then
           local(ic) = cluster(i1,j1,k1)

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters  

           scenario = 2
           Return
           End If
        End If

     k1=k0; j1=j0; i1=i0+1
     If(i1>size(cluster,1)) i1 = 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then
           local(ic) = cluster(i1,j1,k1)

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters

           scenario = 2
           Return
           End If
       End If

     Return

    End Subroutine reveal_local

!---------------------------------------------------------------------
! Subroutine which considers scenario 1 where neighbouring sites are not occupied
! or occupied but not assigned to any clusters
!---------------------------------------------------------------------

    Subroutine scenario_1(i0,j0,k0, cluster,ngrid, cl, trcl, nc)
        Integer, Intent(In)                         :: i0,j0,k0
        Integer, Dimension(:,:,:), Intent(InOut) :: cluster
        Integer*1, Dimension(:,:,:), Intent(InOut) :: ngrid
        Integer, Dimension(:), Intent(InOut)      :: cl, trcl
        Integer, Intent(InOut)                      :: nc

        nc = nc + 1
        trcl(nc) = nc
        cluster(i0,j0,k0) = nc
    End Subroutine scenario_1

!---------------------------------------------------------------------
! Subroutine which considers scenario 2 where some of the sites areoccupied and
! assigned to some clusters
!---------------------------------------------------------------------

    Subroutine scenario_2(i0,j0,k0, cluster,ngrid, cl, trcl, nc)
        Integer, Intent(In)                         :: i0,j0,k0
        Integer, Dimension(:,:,:), Intent(InOut) :: cluster
        Integer*1, Dimension(:,:,:), Intent(InOut) :: ngrid
        Integer, Dimension(:), Intent(InOut)      :: cl, trcl
        Integer, Intent(InOut)                      :: nc
        Integer                                     :: i,j,k,i1,j1,k1, lowest, current, trlowest

        lowest = huge(0); trlowest = huge(0)
        ! Six neighbouring sites

        k1=k0-1; j1=j0; i1=i0            
        If(k1<1) k1=size(cluster,3)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))

           End If
        End If

        k1=k0+1; j1=j0; i1=i0            
        If(k1>size(cluster,3)) k1= 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)               
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))

           End If
        End If

        k1=k0; j1=j0-1; i1=i0
        If(j1<1) j1=size(cluster,2)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)               
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))           

           End If
        End If

        k1=k0; j1=j0+1; i1=i0
        If(j1>size(cluster,2)) j1= 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))

           End If
        End If

        k1=k0; j1=j0; i1=i0-1
        If(i1<1) i1=size(cluster,1)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)               
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))           

           End If
        End If

        k1=k0; j1=j0; i1=i0+1
        If(i1>size(cluster,1)) i1 = 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))

           End If
        End If

        If(lowest< huge(0)) then
        cluster(i0, j0, k0) = lowest
        trcl(lowest) = trlowest
        else
        nc = nc + 1
        trcl(nc) = nc
        cluster(i0,j0,k0) = nc
        return
        end if

        k1=k0-1; j1=j0; i1=i0
        If(k1<1) k1=size(cluster,3)
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End If         

        k1=k0+1; j1=j0; i1=i0 
        If(k1>size(cluster,3)) k1 = 1
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest 
        End If

        k1=k0; j1=j0-1; i1=i0 
        If(j1<1) j1=size(cluster,2)
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End If

        k1=k0; j1=j0+1; i1=i0        
        If(j1>size(cluster,2)) j1 = 1
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End IF


        k1=k0; j1=j0; i1=i0-1           
        If(i1<1) i1=size(cluster,1)
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End If  

        k1=k0; j1=j0; i1=i0+1           
        If(i1>size(cluster,1)) i1 = 1
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End If

        Return


    End Subroutine scenario_2

!---------------------------------------------------------------------
! Subroutine to update the cluster structure
!---------------------------------------------------------------------

    subroutine permutate(lowest, current, cluster, nc, cl)
        Integer, Intent(In)                         :: lowest, current
        Integer, Dimension(:,:,:), Intent(InOut) :: cluster
        Integer, Dimension(:), Intent(InOut)      :: cl
        Integer, Intent(InOut)                      :: nc
        Integer                                     :: i,j,k, LX, LY, LZ

        LX = size(cluster,1)
        LY = size(cluster,2)
        LZ = size(cluster,3)

        do k=1, LZ
            do j=1, LY
                do i=1, LX

                    If(cluster(i,j,k)==current) Then
                        cl(lowest) = cl(lowest) + 1
                        cluster(i,j,k) = lowest
                    End If

                end do
            end do
        end do

        ! Eliminate current cluster by copying last cluster (nc) into it

        cl(current) = cl(nc)

        do k=1, LZ
            do j=1, LY
                do i=1, LX

                    If(cluster(i,j,k)==nc) Then
                        cluster(i,j,k) = current
                    End If

                end do
            end do
        end do
        cl(nc) = 0
        nc = nc - 1
    End Subroutine permutate

!----------------------------------------------------------------------------
!Subroutine to determine if a percolating cluster exists in x y or z-direction
!----------------------------------------------------------------------------

!There could be more than one spanning cluster in each direction
!Useful output: if there exists a spanning cluster

!Look at cluster information and find out if any one of the clusters has an available site in every x y or z-axis position.

    subroutine span(lattice_in, n_sites, nn_sites, cluster, cl, nc, cl_summary)

        Integer*1, Dimension(:,:,:), Intent(InOut) :: lattice_in
        Integer, Dimension(:), Intent(InOut), Optional :: n_sites     ! the sites of the spanning clusters, if wanted
        Integer, Intent(InOut), Optional          :: nn_sites
        Integer, Dimension(:,:,:), Intent(In)     :: cluster
        Integer, Dimension(:), Intent(In)         :: cl
        Integer, Intent(In)                       :: nc
        Integer, Dimension(:,:), Intent(Out)      :: cl_summary
        Integer                                   :: x_span, y_span, z_span,potentialspan !x_span=0 if there is no spanning cluster
        Integer                                   :: n,i,j,k, LX, LY, LZ, ic, icount, spanning, i1, nn
        Integer, Dimension(:), allocatable        :: cand
        Logical(kind=1), Dimension(:,:), allocatable :: xo, yo, zo
        Integer, Dimension(:), allocatable :: x_array, y_array, z_array
        Logical :: attempt

        LX = size(cluster,1)
        LY = size(cluster,2)
        LZ = size(cluster,3)
        allocate(x_array(LX), y_array(LY), z_array(LZ))
        attempt = .False.
        ic = 0
        cl_summary = 0

        ! Which planes each candidate cluster occupies, found in one pass over the grid
        ! (upstream scanned the grid up to three times per candidate); the test per
        ! candidate, in label order, is unchanged
        Call plane_occupancy(cluster, cl, nc, cand, xo, yo, zo)
        do n=1, nc
            if(cand(n) == 0) cycle
            attempt = .True.
            potentialspan = n
            spanning = 0
            if(all(xo(:, cand(n)))) spanning = spanning + 1
            if(all(yo(:, cand(n)))) spanning = spanning + 1
            if(all(zo(:, cand(n)))) spanning = spanning + 1
            if(spanning>0) then
                ic = ic + 1
                cl_summary(1, 1)      = ic
                cl_summary(ic+1, 1)   = potentialspan
                cl_summary(ic+1, 2)   = spanning
            end if
        end do
        deallocate(cand, xo, yo, zo)

        if(attempt.eqv..False.) then
            write(*,*) " The system is NOT percolated in ANY direction "
            spanning = 0

            deallocate(x_array, y_array, z_array)
            return
        end if

        lattice_in = 0
        nn = 0
        if(present(n_sites)) n_sites = 0
        icount = 0
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    icount = icount + 1
                    do i1=1, ic
                        if(cluster(i,j,k) ==  cl_summary(i1+1, 1)) then
                            lattice_in(i,j,k) = 1
                            nn = nn + 1
                            if(present(n_sites)) n_sites(nn) = icount
                        end if
                    end do
                end do
            end do
        end do
        if(present(nn_sites)) nn_sites = nn
        deallocate(x_array, y_array, z_array)
        
    end subroutine span

!----------------------------------------------------------------------------
!Subroutine to determine if a percolating cluster exists in x y or z-direction
!----------------------------------------------------------------------------

!There could be more than one spanning cluster in each direction
!Useful output: if there exists a spanning cluster

!Look at cluster information and find out if any one of the clusters has an available site in every x y or z-axis position.

    subroutine span_simple(lattice_in, cluster, cl, nc, cl_summary, spanning)

        Integer*1, Dimension(:,:,:), Intent(InOut)  :: lattice_in
        Integer, Dimension(:,:,:), Intent(InOut)  :: cluster(:,:,:)
        Integer, Dimension(:), Intent(In)         :: cl
        Integer, Intent(In)                       :: nc
        Integer, Dimension(:,:), Intent(Out)      :: cl_summary
        Integer, Intent(InOut)                    :: spanning
        Integer                                   :: x_span, y_span, z_span,potentialspan !x_span=0 if there is no spanning cluster
        Integer                                   :: n,i,j,k, LX, LY, LZ, ic
        Integer, Dimension(:), allocatable        :: cand
        Logical(kind=1), Dimension(:,:), allocatable :: xo, yo, zo
        Integer, Dimension(:), allocatable  :: x_array, y_array, z_array
        Logical :: attempt

        LX = size(cluster,1)
        LY = size(cluster,2)
        LZ = size(cluster,3)
        allocate(x_array(LX), y_array(LY), z_array(LZ))
        attempt = .False.
        ic = 0
        cl_summary = 0

        ! Which planes each candidate cluster occupies, found in one pass over the grid
        ! (upstream scanned the grid up to three times per candidate); the test per
        ! candidate, in label order, is unchanged
        Call plane_occupancy(cluster, cl, nc, cand, xo, yo, zo)
        do n=1, nc
            if(cand(n) == 0) cycle
            attempt = .True.
            potentialspan = n
            spanning = 0
            if(all(xo(:, cand(n)))) spanning = spanning + 1
            if(all(yo(:, cand(n)))) spanning = spanning + 1
            if(all(zo(:, cand(n)))) spanning = spanning + 1
            if(spanning>0) then
                ic = ic + 1
                cl_summary(1, 1)      = ic
                cl_summary(ic+1, 1)   = potentialspan
                cl_summary(ic+1, 2)   = spanning
                deallocate(x_array, y_array, z_array, cand, xo, yo, zo)
                return
            end if
        end do
        deallocate(cand, xo, yo, zo)

        if(attempt.eqv..False.) then
        !    write(*,*) " The system is NOT percolated in ANY direction "
            spanning = 0
        end if

        deallocate(x_array, y_array, z_array)

    end subroutine span_simple

!----------------------------------------------------------------------------
! For each cluster big enough to span (at least LX, LY or LZ sites), cand(n) numbers it
! and xo(i, cand(n)), yo(j, ...) and zo(k, ...) say whether it has a site in plane
! i, j or k. One pass over the grid.
!----------------------------------------------------------------------------

    subroutine plane_occupancy(cluster, cl, nc, cand, xo, yo, zo)
        Integer, Dimension(:,:,:), Intent(In)                  :: cluster
        Integer, Dimension(:), Intent(In)                      :: cl
        Integer, Intent(In)                                    :: nc
        Integer, Dimension(:), allocatable, Intent(Out)        :: cand
        Logical(kind=1), Dimension(:,:), allocatable, Intent(Out) :: xo, yo, zo
        Integer                                                :: n, m, i, j, k, c, LX, LY, LZ

        LX = size(cluster,1)
        LY = size(cluster,2)
        LZ = size(cluster,3)
        allocate(cand(max(1, nc)))
        cand = 0
        m = 0
        do n=1, nc
            if(cl(n)>=LX.or.cl(n)>=LY.or.cl(n)>=LZ) then
                m = m + 1
                cand(n) = m
            end if
        end do
        allocate(xo(LX, max(1, m)), yo(LY, max(1, m)), zo(LZ, max(1, m)))
        xo = .False.
        yo = .False.
        zo = .False.
        if(m == 0) return
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    c = cluster(i,j,k)
                    if(c < 1 .or. c > nc) cycle
                    c = cand(c)
                    if(c == 0) cycle
                    xo(i, c) = .True.
                    yo(j, c) = .True.
                    zo(k, c) = .True.
                end do
            end do
        end do
    end subroutine plane_occupancy


!---------------------------------------------------------------------
! Subroutine which searches for two clusters sharing a site and
! reconnects them into one large cluster
!---------------------------------------------------------------------

    subroutine update(i0,j0,k0, cluster,ngrid, cl, nc)
        Integer, Intent(In)                         :: i0,j0,k0
        Integer, Dimension(:,:,:), Intent(InOut) :: cluster
        Integer*1, Dimension(:,:,:), Intent(InOut) :: ngrid
        Integer, Dimension(:), Intent(InOut)      :: cl
        Integer, Intent(InOut)                      :: nc
        Integer                                     :: i,j,k,i1,j1,k1, lowest

        lowest = huge(0)
        !First we find the lowest cluster observed among the neighbours
        do k=k0-1, k0+1, 1
            do j=j0-1, j0+1, 1
                do i=i0-1, i0+1, 1

                    if (k/=k0.and.(i/=i0.or.j/=j0)) cycle
                    if (j/=j0.and.(k/=k0.or.i/=i0)) cycle
                    if (i/=i0.and.(j/=j0.or.k/=k0)) cycle

                    i1 = i
                    j1 = j
                    k1 = k

                    If(i1<1) i1=size(cluster,1)
                    If(i1>size(cluster,1)) i1= 1
                    If(j1<1) j1=size(cluster,2)
                    If(j1>size(cluster,2)) j1= 1
                    If(k1<1) k1=size(cluster,3)
                    If(k1>size(cluster,3)) k1= 1

                    If(ngrid(i1,j1,k1)==1) Then
                        If(cluster(i1,j1,k1)/=0) Then
                            If(lowest>cluster(i1,j1,k1)) Then
                                lowest = cluster(i1,j1,k1)
                            End If
                        End If
                    End If

                end do
            end do
        end do

        If(lowest==0) Then
            print*, 'something is wrong, it is not scenario 2'
            stop
        End If

        do k=k0-1, k0+1, 1
            do j=j0-1, j0+1, 1
                do i=i0-1, i0+1, 1

                    if (k/=k0.and.(i/=i0.or.j/=j0)) cycle
                    if (j/=j0.and.(k/=k0.or.i/=i0)) cycle
                    if (i/=i0.and.(j/=j0.or.k/=k0)) cycle

                    i1 = i
                    j1 = j
                    k1 = k

                    If(i1<1) i1=size(cluster,1)
                    If(i1>size(cluster,1)) i1= 1
                    If(j1<1) j1=size(cluster,2)
                    If(j1>size(cluster,2)) j1= 1
                    If(k1<1) k1=size(cluster,3)
                    If(k1>size(cluster,3)) k1= 1

                    If((i1/=i0.or.j1/=j0.or.k1/=k0).and.cluster(i1,j1,k1)==0) cycle

                    If(ngrid(i1,j1,k1)==1) Then
                        If(cluster(i1,j1,k1)==0) Then
                            cl(lowest) = cl(lowest) + 1
                            cluster(i1,j1,k1) = lowest
                        ElseIf(cluster(i1,j1,k1)==lowest)Then
                            cycle
                        Else
                            Call permutate(lowest, cluster(i1,j1,k1), cluster, nc, cl)
                        End If
                    End If

                end do
            end do
        end do
    End Subroutine update

!---------------------------------------------------------------------
! Poreblazer's own labelling needs the whole grid on one rank
!---------------------------------------------------------------------

    Subroutine require_whole_grid(lattice_in)
        Integer*1, Dimension(:,:,:), Intent(In) :: lattice_in
        If(perc_nz /= 0 .and. perc_nz /= size(lattice_in, 3)) &
            Call par_abort("Poreblazer's percolation labelling needs one MPI rank; use exact labelling "// &
                           "(1 after the visualisation option in defaults.dat) to split the grid")
    End Subroutine require_whole_grid

!---------------------------------------------------------------------
! Exact percolation analysis of this rank's slab of the lattice, in parallel.
!
! The result is that of clusteranalysis_exact followed by span (all_spanning) or
! span_simple on the whole grid, on any number of threads and ranks:
!
!  1. Each thread labels its own block of planes with union-find, in which the root of
!     every set is its smallest site: so a cluster's root is its first site in scan
!     order, the site clusteranalysis_exact numbers clusters by.
!  2. The blocks are joined across the planes between them, serially, and every site
!     is pointed straight at its root, in parallel.
!  3. Across ranks: each rank sends the roots of its bottom plane to the rank below (the
!     halo; periodic), which pairs them with its top plane. The pairs, a few per
!     cluster that crosses a slab boundary, are shared by every rank, which joins them
!     with union-find in the same way: a cluster's global root is its first site in the
!     whole grid. Clusters inside one slab need nothing from other ranks.
!  4. A cluster is a candidate if it has at least min(LX, LY, LZ) sites, and spans a
!     direction if it occupies every plane across it: as span does. Clusters crossing
!     ranks add up their sizes and their occupied planes over the ranks.
!  5. Clusters are taken in order of their first site, which is clusteranalysis_exact's
!     label order: span_simple's answer is the first spanning cluster's.
!
! Memory: 4 bytes per site of the slab (parent), instead of the 8 of clusteranalysis_exact
! and span (cluster labels and parent), and no table sized by the number of clusters.
!---------------------------------------------------------------------

    Subroutine exact_percolation(lattice_in, cl_summary, spanning, all_spanning, n_sites, nn_sites)
        Integer*1, Dimension(:,:,:), Intent(InOut)     :: lattice_in
        Integer, Dimension(:,:), Intent(Out)           :: cl_summary
        Integer, Intent(Out)                           :: spanning
        Logical, Intent(In)                            :: all_spanning
        Integer, Dimension(:), Intent(InOut), Optional :: n_sites
        Integer, Intent(InOut), Optional               :: nn_sites

        Integer                   :: LX, LY, LZ, NZ, z0, i, j, k, s, r, c, m, n, nl, ns, ic, minl, nedge, nsp
        Integer*8                 :: plane, base, nloc, gid, ncand
        Logical                   :: wrapz
        Integer, allocatable      :: parent(:), sidx(:), gpar(:)
        Integer*8, allocatable    :: sendp(:), recvp(:), edges(:), alledges(:), nodes(:), gsize(:), spans(:), allspans(:)
        Integer*1, allocatable    :: locc(:,:), socc(:,:)
        Integer, allocatable      :: ldims(:), sdims(:)
        Integer*8, allocatable    :: ea(:), eb(:)

        LX = size(lattice_in, 1)
        LY = size(lattice_in, 2)
        LZ = size(lattice_in, 3)
        NZ = LZ
        z0 = 1
        If(perc_nz /= 0) Then
            NZ = perc_nz
            z0 = perc_z0
        End If
        wrapz = (NZ == LZ)
        plane = int(LX, 8)*LY
        nloc = plane*LZ
        If(nloc > huge(0)) Call par_abort("exact_percolation: more than 2**31 grid cubes on one rank; use more ranks")
        base = int(z0 - 1, 8)*plane
        minl = min(LX, LY, NZ)
        cl_summary = 0
        spanning = 0

        allocate(parent(nloc))
        Call label_local(lattice_in, parent, wrapz)

        ! 3. Clusters that cross slab boundaries: pairs of roots (gids, global site numbers)
        nedge = 0
        allocate(nodes(0), gpar(0))
        If(.not. wrapz) Then
            allocate(sendp(plane), recvp(plane), ea(plane), eb(plane))
            do j=1, LY
                do i=1, LX
                    s = i + LX*(j-1)
                    sendp(s) = 0
                    If(lattice_in(i,j,1) == 1) sendp(s) = base + parent(s)
                end do
            end do
            Call par_shift_down_i8(sendp, recvp)
            do j=1, LY
                do i=1, LX
                    s = i + LX*(j-1)
                    If(recvp(s) == 0 .or. lattice_in(i,j,LZ) /= 1) cycle
                    nedge = nedge + 1
                    ea(nedge) = base + parent(s + int(plane)*(LZ-1))
                    eb(nedge) = recvp(s)
                end do
            end do
            deallocate(sendp, recvp)
            Call sort_pairs(ea, eb, nedge)
            ! unique pairs, packed for sharing
            allocate(edges(2*nedge))
            m = 0
            do n=1, nedge
                If(n > 1) Then
                    If(ea(n) == ea(n-1) .and. eb(n) == eb(n-1)) cycle
                End If
                edges(2*m+1) = ea(n)
                edges(2*m+2) = eb(n)
                m = m + 1
            end do
            deallocate(ea, eb)
            Call par_allgather_i8(edges(1:2*m), alledges)
            deallocate(edges)
            ! the roots that cross, sorted, joined by the pairs (the smallest gid is the root)
            deallocate(nodes, gpar)
            allocate(nodes(size(alledges)))
            nodes = alledges
            n = size(nodes)
            Call sort_values_i8(nodes, n)
            m = 0
            do i=1, n
                If(i > 1) Then
                    If(nodes(i) == nodes(m)) cycle
                End If
                m = m + 1
                nodes(m) = nodes(i)
            end do
            Call shrink_i8(nodes, m)
            allocate(gpar(m))
            do i=1, m
                gpar(i) = i
            end do
            do n=1, size(alledges)/2
                Call unite(gpar, find_i8(nodes, alledges(2*n-1)), find_i8(nodes, alledges(2*n)))
            end do
            do i=1, m
                gpar(i) = findroot(gpar, i)
            end do
            deallocate(alledges)
        End If

        ! 4. Sizes: each root's entry becomes -(its local size)
        !$omp parallel do default(shared) private(s)
        do s=1, int(nloc)
            If(parent(s) == s) parent(s) = -1
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(s, r)
        do s=1, int(nloc)
            r = parent(s)
            If(r <= 0) cycle
            !$omp atomic
            parent(r) = parent(r) - 1
        end do
        !$omp end parallel do

        ! Sizes of the clusters that cross ranks, by global root, summed over the ranks
        allocate(gsize(size(nodes)))
        gsize = 0
        If(size(nodes) > 0) Then
            do s=1, int(nloc)
                If(parent(s) >= 0) cycle
                c = find_i8(nodes, base + s)
                If(c > 0) gsize(gpar(c)) = gsize(gpar(c)) - parent(s)
            end do
            Call par_sum_i8_array(gsize)
        End If

        ! Candidates: shared ones (clusters crossing ranks, by global root) and local ones.
        ! A root's entry becomes -1 (not a candidate), -(1 + l) for local candidate l, or
        ! -(1 + nl + p) for shared candidate p.
        allocate(sidx(size(nodes)))
        sidx = 0
        ns = 0
        do i=1, size(nodes)
            If(gpar(i) == i .and. gsize(i) >= minl) Then
                ns = ns + 1
                sidx(i) = ns
            End If
        end do
        nl = 0
        do s=1, int(nloc)
            If(parent(s) >= 0) cycle
            If(find_i8(nodes, base + s) > 0) cycle
            If(-parent(s) >= minl) nl = nl + 1
        end do
        ncand = par_sum_i8(int(nl, 8)) + ns
        m = 0
        do s=1, int(nloc)
            If(parent(s) >= 0) cycle
            c = find_i8(nodes, base + s)
            If(c > 0) Then
                If(sidx(gpar(c)) > 0) Then
                    parent(s) = -(1 + nl + sidx(gpar(c)))
                Else
                    parent(s) = -1
                End If
            Else If(-parent(s) >= minl) Then
                m = m + 1
                parent(s) = -(1 + m)
            Else
                parent(s) = -1
            End If
        end do

        If(ncand == 0) Then
            ! No cluster could span (span leaves the lattice as it is)
            If(all_spanning) write(*,*) " The system is NOT percolated in ANY direction "
            deallocate(parent, nodes, gpar, gsize, sidx)
            Return
        End If

        ! Planes occupied by each candidate: x, then y, then z (local planes for local
        ! candidates, every plane for shared ones)
        allocate(locc(LX+LY+LZ, max(1, nl)), socc(LX+LY+NZ, max(1, ns)))
        locc = 0
        socc = 0
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    If(lattice_in(i,j,k) /= 1) cycle
                    s = i + LX*(j-1) + int(plane)*(k-1)
                    r = s
                    If(parent(s) > 0) r = parent(s)
                    c = -parent(r) - 1
                    If(c <= 0) cycle
                    If(c <= nl) Then
                        locc(i, c) = 1
                        locc(LX+j, c) = 1
                        locc(LX+LY+k, c) = 1
                    Else
                        socc(i, c-nl) = 1
                        socc(LX+j, c-nl) = 1
                        socc(LX+LY+z0-1+k, c-nl) = 1
                    End If
                end do
            end do
        end do
        If(ns > 0) Call par_or_i1_array(socc, size(socc, 1)*ns)

        ! Spanning: the number of directions in which every plane is occupied
        allocate(ldims(max(1, nl)), sdims(max(1, ns)))
        ldims = 0
        sdims = 0
        do c=1, nl
            If(all(locc(1:LX, c) == 1)) ldims(c) = ldims(c) + 1
            If(all(locc(LX+1:LX+LY, c) == 1)) ldims(c) = ldims(c) + 1
            If(wrapz .and. all(locc(LX+LY+1:LX+LY+LZ, c) == 1)) ldims(c) = ldims(c) + 1
        end do
        do c=1, ns
            If(all(socc(1:LX, c) == 1)) sdims(c) = sdims(c) + 1
            If(all(socc(LX+1:LX+LY, c) == 1)) sdims(c) = sdims(c) + 1
            If(all(socc(LX+LY+1:LX+LY+NZ, c) == 1)) sdims(c) = sdims(c) + 1
        end do
        deallocate(locc, socc)

        ! The spanning clusters, by global root and number of directions: local ones from
        ! this rank, shared ones from the rank that holds their root
        allocate(spans(2*(nl + ns)))
        nsp = 0
        m = 0
        do s=1, int(nloc)
            If(parent(s) >= 0) cycle
            c = -parent(s) - 1
            If(c >= 1 .and. c <= nl) Then
                If(ldims(c) > 0) Then
                    spans(2*nsp+1) = base + s
                    spans(2*nsp+2) = ldims(c)
                    nsp = nsp + 1
                End If
            End If
        end do
        do i=1, size(nodes)
            If(sidx(i) == 0) cycle
            If(sdims(sidx(i)) == 0) cycle
            If(nodes(i) <= base .or. nodes(i) > base + nloc) cycle
            spans(2*nsp+1) = nodes(i)
            spans(2*nsp+2) = sdims(sidx(i))
            nsp = nsp + 1
        end do
        Call par_allgather_i8(spans(1:2*nsp), allspans)
        deallocate(spans)
        nsp = size(allspans)/2
        allocate(ea(max(1, nsp)), eb(max(1, nsp)))
        do n=1, nsp
            ea(n) = allspans(2*n-1)
            eb(n) = allspans(2*n)
        end do
        deallocate(allspans)
        Call sort_pairs(ea, eb, nsp)

        ic = 0
        do n=1, nsp
            If(.not. all_spanning .and. n > 1) exit
            ic = ic + 1
            If(ic + 1 <= size(cl_summary, 1)) Then
                cl_summary(1, 1)    = ic
                cl_summary(ic+1, 1) = ic
                cl_summary(ic+1, 2) = int(eb(n))
            End If
            spanning = int(eb(n))
        end do
        deallocate(ea, eb)

        If(all_spanning) Then
            ! Keep only the sites of spanning clusters; the root's entry becomes -2 if its
            ! cluster spans
            do s=1, int(nloc)
                If(parent(s) >= 0) cycle
                c = -parent(s) - 1
                m = 0
                If(c >= 1 .and. c <= nl) Then
                    m = ldims(c)
                Else If(c > nl) Then
                    m = sdims(c - nl)
                End If
                parent(s) = -1
                If(m > 0) parent(s) = -2
            end do
            n = 0
            If(present(n_sites)) n_sites = 0
            do k=1, LZ
                do j=1, LY
                    do i=1, LX
                        s = i + LX*(j-1) + int(plane)*(k-1)
                        r = s
                        If(lattice_in(i,j,k) == 1) Then
                            If(parent(s) > 0) r = parent(s)
                            If(parent(r) == -2) Then
                                n = n + 1
                                If(present(n_sites)) n_sites(n) = s
                                cycle
                            End If
                        End If
                        lattice_in(i,j,k) = 0
                    end do
                end do
            end do
            If(present(nn_sites)) nn_sites = n
        End If
        deallocate(parent, nodes, gpar, gsize, sidx, ldims, sdims)

    End Subroutine exact_percolation

!---------------------------------------------------------------------
! Union-find labelling of one slab, in parallel: parent(s) is the smallest site of s's
! cluster (0 for an empty site), with six-neighbour connectivity, periodic in x and y,
! and in z if wrapz. Sites are numbered i + LX*(j-1) + LX*LY*(k-1).
!---------------------------------------------------------------------

    Subroutine label_local(ngrid, parent, wrapz)
        Integer*1, Dimension(:,:,:), Intent(In) :: ngrid
        Integer, Dimension(:), Intent(Out)      :: parent
        Logical, Intent(In)                     :: wrapz
        Integer                                 :: LX, LY, LZ, nt, t, kb, ke, i, j, k, s, r, q, plane
        Integer, allocatable                    :: kend(:)

        LX = size(ngrid, 1)
        LY = size(ngrid, 2)
        LZ = size(ngrid, 3)
        plane = LX*LY
        nt = 1
        !$ nt = omp_get_max_threads()
        nt = max(1, min(nt, LZ))
        allocate(kend(0:nt-1))
        do t=0, nt-1
            kend(t) = ((t+1)*LZ)/nt
        end do

        ! 1. Each thread labels its own planes kb..ke
        !$omp parallel do default(shared) private(t, kb, ke, i, j, k, s) schedule(static, 1) num_threads(nt)
        do t=0, nt-1
            kb = (t*LZ)/nt + 1
            ke = kend(t)
            do k=kb, ke
                do j=1, LY
                    do i=1, LX
                        s = i + LX*(j-1) + plane*(k-1)
                        parent(s) = 0
                        If(ngrid(i,j,k) == 1) parent(s) = s
                    end do
                end do
            end do
            do k=kb, ke
                do j=1, LY
                    do i=1, LX
                        If(ngrid(i,j,k) /= 1) cycle
                        s = i + LX*(j-1) + plane*(k-1)
                        If(ngrid(modulo(i,LX)+1,j,k) == 1) Call unite(parent, s, modulo(i,LX)+1 + LX*(j-1) + plane*(k-1))
                        If(ngrid(i,modulo(j,LY)+1,k) == 1) Call unite(parent, s, i + LX*modulo(j,LY) + plane*(k-1))
                        If(k < ke) Then
                            If(ngrid(i,j,k+1) == 1) Call unite(parent, s, s + plane)
                        End If
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! 2. Join the blocks across the planes between them (and z's periodic boundary)
        do t=0, nt-1
            k = kend(t)
            If(k == LZ) Then
                If(.not. wrapz) cycle
                q = 1
            Else
                q = k + 1
            End If
            do j=1, LY
                do i=1, LX
                    If(ngrid(i,j,k) == 1 .and. ngrid(i,j,q) == 1) &
                        Call unite(parent, i + LX*(j-1) + plane*(k-1), i + LX*(j-1) + plane*(q-1))
                end do
            end do
        end do

        ! Point every site at its root. A parent is always a smaller site, so a thread
        ! going up its own planes finds its own sites' parents already done; a parent in
        ! another thread's planes is only ever replaced by one of its ancestors, so the
        ! chain still ends at the root.
        !$omp parallel do default(shared) private(t, s, r, q) schedule(static, 1) num_threads(nt)
        do t=0, nt-1
            do s=((t*LZ)/nt)*plane + 1, kend(t)*plane
                !$omp atomic read
                r = parent(s)
                If(r == 0 .or. r == s) cycle
                do
                    !$omp atomic read
                    q = parent(r)
                    If(q == r) exit
                    r = q
                end do
                !$omp atomic write
                parent(s) = r
            end do
        end do
        !$omp end parallel do
        deallocate(kend)
    End Subroutine label_local

    ! The root of x (the smallest member of its set), halving the path
    Integer Function findroot(parent, x)
        Integer, Dimension(:), Intent(InOut) :: parent
        Integer, Intent(In)                  :: x
        Integer                              :: y
        y = x
        do while(parent(y) /= y)
            parent(y) = parent(parent(y))
            y = parent(y)
        end do
        findroot = y
    End Function findroot

    Subroutine unite(parent, a, b)
        Integer, Dimension(:), Intent(InOut) :: parent
        Integer, Intent(In)                  :: a, b
        Integer                              :: ra, rb
        ra = findroot(parent, a)
        rb = findroot(parent, b)
        If(ra /= rb) parent(max(ra, rb)) = min(ra, rb)
    End Subroutine unite

    ! The position of x in the sorted list v, or 0
    Integer Function find_i8(v, x)
        Integer*8, Dimension(:), Intent(In) :: v
        Integer*8, Intent(In)               :: x
        Integer                             :: lo, hi, mid
        find_i8 = 0
        lo = 1
        hi = size(v)
        do while(lo <= hi)
            mid = lo + (hi - lo)/2
            If(v(mid) == x) Then
                find_i8 = mid
                return
            Else If(v(mid) < x) Then
                lo = mid + 1
            Else
                hi = mid - 1
            End If
        end do
    End Function find_i8

    Subroutine shrink_i8(v, n)
        Integer*8, allocatable, Intent(InOut) :: v(:)
        Integer, Intent(In)                   :: n
        Integer*8, allocatable                :: w(:)
        allocate(w(n))
        w = v(1:n)
        Call move_alloc(w, v)
    End Subroutine shrink_i8

    ! Heap sorts (n log n, in place): of v(1:n), and of the pairs (a, b)(1:n) by a then b
    Subroutine sort_values_i8(v, n)
        Integer*8, Dimension(:), Intent(InOut) :: v
        Integer, Intent(In)                    :: n
        Integer*8, allocatable                 :: dummy(:)
        allocate(dummy(max(1, n)))
        dummy = 0
        Call sort_pairs(v, dummy, n)
        deallocate(dummy)
    End Subroutine sort_values_i8

    Subroutine sort_pairs(a, b, n)
        Integer*8, Dimension(:), Intent(InOut) :: a, b
        Integer, Intent(In)                    :: n
        Integer                                :: i, last
        Integer*8                              :: t
        do i=n/2, 1, -1
            Call sift(i, n)
        end do
        do last=n, 2, -1
            t = a(1); a(1) = a(last); a(last) = t
            t = b(1); b(1) = b(last); b(last) = t
            Call sift(1, last-1)
        end do
    contains
        Logical Function less(i, j)
            Integer, Intent(In) :: i, j
            less = a(i) < a(j) .or. (a(i) == a(j) .and. b(i) < b(j))
        End Function less
        Subroutine sift(start, m)
            Integer, Intent(In) :: start, m
            Integer             :: root, child
            Integer*8           :: t
            root = start
            do while(2*root <= m)
                child = 2*root
                If(child < m) Then
                    If(less(child, child+1)) child = child + 1
                End If
                If(.not. less(root, child)) return
                t = a(root); a(root) = a(child); a(child) = t
                t = b(root); b(root) = b(child); b(child) = t
                root = child
            end do
        End Subroutine sift
    End Subroutine sort_pairs

End Module percolation
