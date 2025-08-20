!> Routines for calving from tabular ice shelves
module ice_shelf_tabular_calving

use constants_mod, only: pi, omega, HLF
use mpp_mod, only : mpp_send, mpp_recv, mpp_max, mpp_npes, mpp_pe, mpp_root_pe, NULL_PE
use mpp_mod, only: COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
use mpp_mod, only: COMM_TAG_5, COMM_TAG_6, COMM_TAG_7, COMM_TAG_8, mpp_sync_self
use mpp_domains_mod, only: mpp_update_domains
use ice_bergs_framework, only : ij_component_of_id, icebergs_gridded, tabular_calving_state
use ice_bergs_framework, only : icebergs, iceberg, bond, delete_iceberg_from_list
use ice_bergs_framework, only : add_new_berg_to_list
use ice_bergs_framework, only : spread_variable_across_cells, sum_up_spread_fields
use ice_bergs_framework, only : hexagon_into_quadrants_using_triangles, Rearth
use ice_bergs_framework, only : square_into_quadrants_using_triangles
use ice_bergs_framework, only : initialize_iceberg_bonds, count_bonds
use ice_bergs_framework, only : find_cell, pos_within_cell, generate_id, rho_seawater
use ice_bergs_framework, only : debug, footloose, connect_all_bonds, delete_all_bonds
use ice_bergs_framework, only : update_halo_calved_tabular_icebergs, assign_n_bonds,transfer_mts_bergs
use fms_mod, only : error_mesg, FATAL, WARNING, stderr

implicit none ; private

real, parameter :: pi_180=pi/180.  !< Converts degrees to radians

public initialize_tabular_calving, process_tabular_calving, ice_shelf_calving_end

contains

!> Initializes ice shelf tabular calving data, parameters, and diagnostics
subroutine initialize_tabular_calving(TC, grd)
  type(tabular_calving_state), pointer    :: TC !< A pointer to the tabular calving structure
  type(icebergs_gridded), pointer :: grd

  if (associated(TC)) then
    call error_mesg('KID, initialize_tabular_calving: ' , &
                          'called with an associated tabular_calving pointer.',FATAL)
    return
  endif

  allocate(TC)
  allocate(TC%c_id(grd%isd:grd%ied,grd%jsd:grd%jed)             , source=0   )
  allocate(TC%calve_mask(grd%isd:grd%ied,grd%jsd:grd%jed)       , source=0.0 )
  allocate(TC%h_shelf(grd%isd:grd%ied,grd%jsd:grd%jed)          , source=0.0 )
  allocate(TC%frac_shelf(grd%isd:grd%ied,grd%jsd:grd%jed)       , source=0.0 )
  allocate(TC%frac_cberg_calved(grd%isd:grd%ied,grd%jsd:grd%jed), source=0.0 )
  allocate(TC%frac_cberg(grd%isd:grd%ied,grd%jsd:grd%jed)       , source=0.0 )
  if (grd%id_pf_area>0) &
    allocate(TC%saved_pf_area(grd%isd:grd%ied,grd%jsd:grd%jed,10), source=0.0 )
end subroutine initialize_tabular_calving

!> Deallocates all memory associated with this module
subroutine ice_shelf_calving_end(TC)
  type(tabular_calving_state), pointer :: TC !< A pointer to the ice shelf calving structure

  if (.not.associated(TC)) return

  if (associated(TC%berg_list))         deallocate(TC%berg_list)
  if (associated(TC%c_id))              deallocate(TC%c_id)
  if (associated(TC%calve_mask))        deallocate(TC%calve_mask)
  if (associated(TC%h_shelf))           deallocate(TC%h_shelf)
  if (associated(TC%frac_shelf))        deallocate(TC%frac_shelf)
  if (associated(TC%frac_cberg_calved)) deallocate(TC%frac_cberg_calved)
  if (associated(TC%frac_cberg))        deallocate(TC%frac_cberg)
  if (associated(TC%saved_pf_area))     deallocate(TC%saved_pf_area)
  if (associated(TC))                   deallocate(TC)
end subroutine ice_shelf_calving_end

!>  routine to initialize iKID icebergs from a tabular calving mask, called from icebergs_run
subroutine process_tabular_calving(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(tabular_calving_state), pointer :: TC !< A pointer to the tabular calving structure

  TC => bergs%TC
  if (associated(TC%berg_list)) deallocate(TC%berg_list)

  !We want to initialize (calve) bonded-particle tabular bergs that cover the "berg cells" defined where
  !mask>0. There may be multiple tabular bergs to calve in the domain, which may overlap multiple cells and
  !PEs. The approach is to create a field of unique IDs that are each associated with a different tabular
  !berg. This ID is the same on all PEs that the corresponding berg may overlap. Then, the max/min x and y
  !coordinates of each berg are passed between the PEs that the berg overlaps. Next, a rectangular array of
  !bonded particles in initialized that spans these max/min coordinates, and which is defined identically
  !and redundantly for each PE that the berg overlaps. On each PE, excess particles which do not overlap any
  !berg cell associated with the current berg are then trimmed off. Initialization of the fully-bonded berg
  !is completed in the iceberg module, where the bonded particles of the berg are connected across PE
  !boundaries to produce the full iceberg.
  ! if (mpp_pe().eq.mpp_root_pe()) print *,'start processing bergs'
  !1) Generate a unique label (TC%c_id) for each berg on the computational domain of a PE
  call initialize_tabular_calving_labels_1PE(bergs%grd, TC)
  !2) fill halos with the unique berg labels. Find connected berg cells, and update them
  !   with the lowest berg label (TC%c_id) of all of the connected berg cells. Repeat until no changes.
  ! if (mpp_pe().eq.mpp_root_pe()) print *,'updating over pes'
  call update_tabular_calving_labels_over_pes(bergs%grd, TC) !%calve_mask, TC%c_id)
  !3) Make a list of each of the i bergs on the local PE domain (TC%berg_list(i,1)), their
  !   min (TC%berg_list(i,2)) max longitude (TC%berg_list(i,3)), and their
  !   min (TC%berg_list(i,4)) max latitude (TC%berg_list(i,5))
  ! if (mpp_pe().eq.mpp_root_pe()) print *,'tabular_berg_info'
  call tabular_berg_info(bergs%grd, TC)
  !4) Initialize iKID icebergs over these bounds, and remove excess particles.
  call ice_shelf_to_bonded_bergs(bergs, TC)
  ! if (mpp_pe().eq.mpp_root_pe()) print *,'done processing bergs'
  !TODO: Note that the code here assumes that any two neighboring cells, each with mask>0, must be part of
  !the same berg. If we want to calve adjacent bergs, we can calve one on the first timestep and the other
  !on the next timestep after the first is fully initialized. Or if using damage, calve both as one berg,
  !interp the damage to the new berg, and break bonds where there is a rift.

  !Alternatively, if we want to allow multiple, adjacent bergs, to calve on a single time step, then we
  !could assign different (positive, non-zero) mask values to differentiate the bergs. After the single berg
  !is calved, we can break bonds where the mask is different (we would not break over halo cells unless the
  !mask indicated, though the mask values may differ on other PEs).  In other words, we would initialize the
  !multiple and adjacent bergs as a single berg, which is which is subsequently broken up into the multiple
  !adjecent bergs. This approach guarantees that the adjacent bergs are initialized without inter-berg
  !particle overlap/separation.

  !Alternatively, we could just have a "background" set of bonded particles that are kept constant over
  !time, and initialize bergs by copying over a subset of these bonded particles as needed to represent the
  !new bergs. However, this approach is not as versatile for controlling particle size. Also, it is not
  !simple to guarantee that the particles extents of the (meridionally-aligned) first and last column would
  !align (zonally).

  !comment out halo-filling of calve_mask for now. If there are multiple, adjacent calving events, we may
  !define them with different calve_mask>0 on each PE, so that the value of calve_mask may differ between
  !PEs for the same berg. However, we make sure each PE calculates its own calve_mask in its own halo cells,
  !and then make a copy of the resulting calve_mask. Then mpp_update_domains the original calve_mask and
  !compare to the copy to backcalculate a consistent value for calve_mask for each berg that can be shared
  !between each PE later.  call mpp_update_domains(TC%calve_mask, G%domain)

  !TODO: Icebergs may overlap with each other or the ice shelf, so make sure mass scaling is done
  !appropriately on new iceberg particles to avoid issues with pressure on the ocean...Also should strongly
  !force bergs away from ice shelves. Perhaps need to ensure that the total pressure on the ocean does not
  !exceed that of a combination of the iceberg(s) and ice shelf mass should it fully-cover the cell (with
  !adjustments for mass-weighting and the percentage of the total unadjusted mass that the pressure of each
  !component exerts on the cell...). Mass will not be conserved at that instant, but ultimately is over
  !time.
end subroutine process_tabular_calving


!> Initializes labels for the grid cells that comprise a tabular iceberg that is about to calve from an ice
!! shelf (i.e. a group of all neighboring grid cells where calving mask > 0). Considers the current PE only.
subroutine initialize_tabular_calving_labels_1PE(grd, TC)
  type(icebergs_gridded), pointer :: grd
  type(tabular_calving_state), pointer :: TC !< A pointer to the tabular calving structure
  integer :: i, j

  TC%c_id(:,:)=0
  do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
    if ((TC%calve_mask(i,j) > 0) .and. (TC%c_id(i,j) == 0)) then
      TC%c_id(i,j) = ij_component_of_id(grd,i,j)
      ! print *,'c_id(i,j) 1',TC%c_id(i,j)
      call label_tabular_bergs(grd, i, j, TC)
      ! print *,'c_id(i,j) 2',TC%c_id(i,j)
    endif
  enddo; enddo

  if (maxval(TC%c_id)>0) then
    print *,'pe,max_cid,ymin,ymax', mpp_pe(),maxval(TC%c_id),&
      minval(grd%lat(grd%isc:grd%iec,grd%jsc:grd%jec)),maxval(grd%lat(grd%isc:grd%iec,grd%jsc:grd%jec))
  endif
end subroutine initialize_tabular_calving_labels_1PE

!> Assigns the same label to all grid cells that comprise a tabular iceberg that is about to calve
!! from an ice shelf (i.e. a group of all neighboring grid cells where calving mask > 0).
recursive subroutine label_tabular_bergs(grd, ic, jc, TC) !, mask, c_id)
  type(icebergs_gridded), pointer :: grd
  type(tabular_calving_state), pointer :: TC !< A pointer to the tabular calving structure
  integer,               intent(in)    :: ic  !< The i-index of the input cell
  integer,               intent(in)    :: jc  !< The j-index of the input cell

  integer :: i, j, is, ie, js, je

  is=max(ic-1,grd%isd); ie=min(ic+1,grd%ied)
  js=max(jc-1,grd%jsd); je=min(jc+1,grd%jed)

  do j=js,je; do i=is,ie

    if ((TC%calve_mask(i,j) > 0) .and. (TC%c_id(i,j) /= TC%c_id(ic,jc))) then

      TC%c_id(i,j) = max(TC%c_id(ic,jc),TC%c_id(i,j))
      ! print *,'c_id(i,j) 3',TC%c_id(i,j)
      call label_tabular_bergs(grd, i, j, TC)
    endif
  enddo; enddo
end subroutine label_tabular_bergs

!> Adjusts labels of tabular icebergs on the grid so that they are consistent between all PEs that
!! the tabular bergs may overlap
subroutine update_tabular_calving_labels_over_pes(grd, TC) !mask, c_id)
  type(icebergs_gridded), pointer :: grd
  type(tabular_calving_state), pointer :: TC !< A pointer to the tabular calving structure
  integer :: i, j, k, i2, j2
  integer :: change,max_c_id

  change=1
  do while (change>0)
    change=0
    call mpp_update_domains(TC%c_id, grd%domain)

    do k=1,2
      if (k==1) then
        i=grd%isc; i2=i-1
      else
        i=grd%iec; i2=i+1
      endif

      do j=grd%jsd,grd%jed
        if (TC%c_id(i,j) /= 0) then
          if (TC%c_id(i2,j) > TC%c_id(i,j)) then
            TC%c_id(i,j) = TC%c_id(i2,j)
            ! print *,'TC%c_id(i,j) 4',TC%c_id(i,j)
            change=change+1
            call label_tabular_bergs(grd, i, j, TC)
          endif
        endif
      enddo

      if (k==1) then
        j=grd%jsc; j2=j-1
      else
        j=grd%jec; j2=j+1
      endif

      do i=grd%isd,grd%ied
        if (TC%c_id(i,j) /= 0) then
          if (TC%c_id(i,j2) > TC%c_id(i,j)) then
            ! print *,'TC%c_id(i,j) 5a',TC%c_id(i,j)
            TC%c_id(i,j) = TC%c_id(i,j2)
            ! print *,'TC%c_id(i,j) 5b',TC%c_id(i,j)
            change=change+1
            call label_tabular_bergs(grd, i, j, TC)
          endif
        endif
      enddo
    enddo

    if (change>0) print *,'pe,change,max_c_id',mpp_pe(),change,maxval(TC%c_id)

    call mpp_max(change)
    max_c_id=maxval(TC%c_id)
    call mpp_max(max_c_id)
    if (change>0 .and. mpp_pe()==mpp_root_pe()) print *,'change,max_c_id',change,max_c_id
  enddo

end subroutine update_tabular_calving_labels_over_pes

!> Finds the unique bergs of a computational domain field, and their global extent
!! (min and max latitude and longitude)
subroutine tabular_berg_info(grd, TC)
  type(icebergs_gridded), pointer :: grd
  type(tabular_calving_state), pointer, intent(inout) :: TC !< A pointer to the tabular calving structure
  integer :: bcount ! number of unique bergs
  ! local variables
  real, dimension(:), allocatable :: tmp
  !integer, pointer :: c_id(:,:)
  integer :: n, i
  real :: min_val, max_val
  integer :: pe_N,pe_S,pe_E,pe_W

  !c_id=>TC%c_id
  n=(grd%iec-grd%isc) * (grd%jec-grd%jsc)
  allocate(tmp(n))
  tmp(:) = 0.
  min_val = minval(TC%c_id(grd%isc:grd%iec,grd%jsc:grd%jec))
  max_val = maxval(TC%c_id(grd%isc:grd%iec,grd%jsc:grd%jec))
  bcount = 0
  ! print *,'pe,min_val,max_val',mpp_pe(),int(min_val),int(max_val)

  if (max_val>0) then
    do while (min_val<max_val)
      bcount = bcount+1
      min_val = minval(TC%c_id(grd%isc:grd%iec,grd%jsc:grd%jec), mask=TC%c_id(grd%isc:grd%iec,grd%jsc:grd%jec)>min_val)
      tmp(bcount) = real(min_val)
    enddo
  endif

  if (associated(TC%berg_list)) deallocate(TC%berg_list)

  if (bcount>0) then
    allocate(TC%berg_list(bcount,5))

    do i = 1,bcount

      !the unique berg
      TC%berg_list(i,1) = real(tmp(i))

      !minlon
      TC%berg_list(i,2) = minval(grd%lon(grd%isc-1:grd%iec-1,grd%jsc-1:grd%jec-1), &
        mask=TC%c_id(grd%isc:grd%iec,grd%jsc:grd%jec)==tmp(i))
      !maxlon
      TC%berg_list(i,3) = maxval(grd%lon(grd%isc:grd%iec,grd%jsc:grd%jec), &
        mask=TC%c_id(grd%isc:grd%iec,grd%jsc:grd%jec)==tmp(i))
      !minlat
      TC%berg_list(i,4) = minval(grd%lat(grd%isc-1:grd%iec-1,grd%jsc-1:grd%jec-1), &
        mask=TC%c_id(grd%isc:grd%iec,grd%jsc:grd%jec)==tmp(i))
      !maxlat
      TC%berg_list(i,5) = maxval(grd%lat(grd%isc:grd%iec,grd%jsc:grd%jec), &
        mask=TC%c_id(grd%isc:grd%iec,grd%jsc:grd%jec)==tmp(i))
    enddo
  else
    TC%berg_list=>NULL()
  endif

  TC%berg_pe_count=bcount

  !update the local coordinate bounds for each berg in TC%berg_list with the
  !global coordinate bounds for each berg across all PEs
  call update_berg_lists_on_all_pes(grd, TC)

  deallocate(tmp)
end subroutine tabular_berg_info

!> Updates TC%berg_list
recursive subroutine update_berg_lists_on_all_pes(grd, TC)
  type(icebergs_gridded), pointer :: grd
  type(tabular_calving_state), pointer, intent(inout) :: TC !< A pointer to the tabular calving structure
  integer :: i1(4), i2(4), j1(4), j2(4), nbergs(4)
  integer, allocatable :: tracker(:,:)
  real, allocatable :: buffer(:)
  integer :: n,k
  integer :: changes, localchanges
  integer :: data_rcvd_from_e, data_rcvd_from_w, data_rcvd_from_n, data_rcvd_from_s
  !integer, pointer :: c_id(:,:)
  !real, pointer :: pebl(:,:)


  !Halo range for each cardinal direction

  i1(1)=grd%iec+1; i2(1)=grd%ied;   j1(1)=grd%jsd;   j2(1)=grd%jed   !E
  i1(2)=grd%isd;   i2(2)=grd%isc-1; j1(2)=grd%jsd;   j2(2)=grd%jed   !W
  i1(3)=grd%isd;   i2(3)=grd%ied;   j1(3)=grd%jec+1; j2(3)=grd%jed   !N
  i1(4)=grd%isd;   i2(4)=grd%ied;   j1(4)=grd%jsd;   j2(4)=grd%jsc-1 !S

  !c_id=>TC%c_id
  !pebl=>TC%berg_list

  !each row of tracker is a berg on the PE, and each column corresponds to a cardinal direction.
  !where tracker == 1, the berg in that row is to be sent to the PE in the corresponding direction.
  allocate(tracker(max(TC%berg_pe_count,1),4)); tracker(:,:)=0

  if (associated(TC%berg_list)) then
    do n = 1,TC%berg_pe_count
      do k=1,4
        if (any(TC%c_id(i1(k):i2(k),j1(k):j2(k))==int(TC%berg_list(n,1)))) tracker(n,k)=1
        if (any(TC%c_id(i1(k):i2(k),j1(k):j2(k))==int(TC%berg_list(n,1)))) tracker(n,k)=1
      enddo
    enddo
  endif

  !number of bergs to send to each direction
  nbergs(1:4)=sum(tracker,1)

  ! if (mpp_pe().eq.mpp_root_pe()) then
  !   print *,'tracker(:,1)',tracker(:,1)
  !   print *,'tracker(:,2)',tracker(:,2)
  !   print *,'tracker(:,3)',tracker(:,3)
  !   print *,'tracker(:,4)',tracker(:,4)
  ! endif

  ! print *,'pe,berg_count,nbergs',mpp_pe(),TC%berg_pe_count,nbergs(1:4)
  call mpp_sync_self()

  changes = 1
  localchanges = 1
  do while (changes/=0)
    ! if (mpp_pe().eq.mpp_root_pe()) print *,'changes',changes
    changes=0

    !send bergs east/west
    if (grd%pe_E.ne.NULL_PE) then
      if (localchanges>0) then
        call mpp_send(nbergs(1)*5, plen=1, to_pe=grd%pe_E, tag=COMM_TAG_1)
        ! if (mpp_pe().eq.mpp_root_pe()) print *,'done send E'
        if (nbergs(1).gt.0) then
          allocate(buffer(nbergs(1)*5))
          call pack_tabular_buffer(TC%berg_pe_count, nbergs(1)*5, 1, TC%berg_list, tracker, buffer)
          call mpp_send(buffer, nbergs(1)*5, grd%pe_E, tag=COMM_TAG_2)
          deallocate(buffer)
        endif
      else
        call mpp_send(0, plen=1, to_pe=grd%pe_E, tag=COMM_TAG_1)
      endif
    endif
    if (grd%pe_W.ne.NULL_PE) then
      if (localchanges>0) then
        call mpp_send(nbergs(2)*5, plen=1, to_pe=grd%pe_W, tag=COMM_TAG_3)
        ! if (mpp_pe().eq.mpp_root_pe()) print *,'done send W'
        if (nbergs(2).gt.0) then
          allocate(buffer(nbergs(2)*5))
          call pack_tabular_buffer(TC%berg_pe_count, nbergs(2)*5, 2, TC%berg_list, tracker, buffer)
          call mpp_send(buffer, nbergs(2)*5, grd%pe_W, tag=COMM_TAG_4)
          deallocate(buffer)
        endif
      else
        call mpp_send(0, plen=1, to_pe=grd%pe_W, tag=COMM_TAG_3)
      endif
    endif

    !receive bergs from west/east
    if (grd%pe_W.ne.NULL_PE) then
      data_rcvd_from_w=0
      call mpp_recv(data_rcvd_from_w, glen=1, from_pe=grd%pe_W, tag=COMM_TAG_1)
      ! if (mpp_pe().eq.mpp_root_pe()) print *,'done receive W'
      if (data_rcvd_from_w.gt.0) then
        allocate(buffer(data_rcvd_from_w))
        call mpp_recv(buffer, data_rcvd_from_w, grd%pe_W, tag=COMM_TAG_2)
        call unpack_tabular_buffer_and_update_bounds(TC%berg_pe_count, data_rcvd_from_w, &
                                                     TC%berg_list, buffer, changes)
        deallocate(buffer)
      endif
    endif
    if (grd%pe_E.ne.NULL_PE) then
      data_rcvd_from_e=0
      call mpp_recv(data_rcvd_from_e, glen=1, from_pe=grd%pe_E, tag=COMM_TAG_3)
      ! if (mpp_pe().eq.mpp_root_pe()) print *,'done receive E'
      if (data_rcvd_from_e.gt.0) then
        allocate(buffer(data_rcvd_from_e))
        call mpp_recv(buffer, data_rcvd_from_e, grd%pe_E, tag=COMM_TAG_4)
        call unpack_tabular_buffer_and_update_bounds(TC%berg_pe_count, data_rcvd_from_e, TC%berg_list, buffer, changes)
        deallocate(buffer)
      endif
    endif

    call mpp_sync_self()
    localchanges=max(localchanges,changes)
    ! if (mpp_pe().eq.mpp_root_pe()) print *,'localchanges 0.5',localchanges

    !send bergs north/south
    if (grd%pe_N.ne.NULL_PE) then
      if (localchanges>0) then
        call mpp_send(nbergs(3)*5, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_5)
        ! if (mpp_pe().eq.mpp_root_pe()) print *,'done send N'
        if (nbergs(3).gt.0) then
          allocate(buffer(nbergs(3)*5))
          call pack_tabular_buffer(TC%berg_pe_count, nbergs(3)*5, 3, TC%berg_list, tracker, buffer)
          call mpp_send(buffer, nbergs(3)*5, grd%pe_N, tag=COMM_TAG_6)
          deallocate(buffer)
        endif
      else
        call mpp_send(0, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_5)
      endif
    endif
    if (grd%pe_S.ne.NULL_PE) then
      if (localchanges>0) then
        call mpp_send(nbergs(4)*5, plen=1, to_pe=grd%pe_S, tag=COMM_TAG_7)
        ! if (mpp_pe().eq.mpp_root_pe()) print *,'done send S'
        if (nbergs(4).gt.0) then
          allocate(buffer(nbergs(4)*5))
          call pack_tabular_buffer(TC%berg_pe_count, nbergs(4)*5, 4, TC%berg_list, tracker, buffer)
          call mpp_send(buffer, nbergs(4)*5, grd%pe_S, tag=COMM_TAG_8)
          deallocate(buffer)
        endif
      else
        call mpp_send(0, plen=1, to_pe=grd%pe_S, tag=COMM_TAG_7)
      endif
    endif

    !receive bergs north/south
    if (grd%pe_S.ne.NULL_PE) then
      data_rcvd_from_s=0
      call mpp_recv(data_rcvd_from_s, glen=1, from_pe=grd%pe_S, tag=COMM_TAG_5)
      ! if (mpp_pe().eq.mpp_root_pe()) print *,'done receive S'
      if (data_rcvd_from_s.gt.0) then
        allocate(buffer(data_rcvd_from_s))
        call mpp_recv(buffer, data_rcvd_from_s, grd%pe_S, tag=COMM_TAG_6)
        call unpack_tabular_buffer_and_update_bounds(TC%berg_pe_count, data_rcvd_from_s, TC%berg_list, buffer, changes)
        deallocate(buffer)
      endif
    endif
    if (grd%pe_N.ne.NULL_PE) then
      data_rcvd_from_n=0
      call mpp_recv(data_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_7)
      ! if (mpp_pe().eq.mpp_root_pe()) print *,'done receive N'
      if (data_rcvd_from_n.gt.0) then
        allocate(buffer(data_rcvd_from_n))
        call mpp_recv(buffer, data_rcvd_from_n, grd%pe_N, tag=COMM_TAG_8)
        call unpack_tabular_buffer_and_update_bounds(TC%berg_pe_count, data_rcvd_from_n, TC%berg_list, buffer, changes)
        deallocate(buffer)
      endif
    endif

    call mpp_sync_self()
    localchanges=changes
    ! if (mpp_pe().eq.mpp_root_pe()) print *,'localchanges 1',localchanges
    call mpp_max(changes)
    ! if (mpp_pe().eq.mpp_root_pe()) print *,'changes 1',changes
  enddo

  if (allocated(tracker)) deallocate(tracker)
end subroutine update_berg_lists_on_all_pes

!> pack the buffer with the info for the bergs in the c_id array being sent to another PE
subroutine pack_tabular_buffer(lbergs, bberg_data, dir, pebl, tracker, buffer)
  integer :: lbergs !< number of local bergs on the current PE
  integer :: bberg_data !< count of berg data being sent in the buffer to another PE
  integer :: dir !< direction the bergs are being sent
  real, pointer :: pebl(:,:) !< array of bergs on the current PE and their lat/lon bounds
  integer :: tracker(lbergs,4) !< tracks the direction to send bergs from the current PE
  real :: buffer(bberg_data) !< the buffer being packed with bergs to send to another PE
  integer :: k, i

  buffer(:)=0.
  i=1
  do k=1,lbergs
    if (tracker(k,dir)==1) then
      buffer(i:i+4)=pebl(k,1:5)
      i=i+5
    endif
  enddo

end subroutine pack_tabular_buffer

!> Unpack the buffer with the info for the bergs in the c_id array being sent from another PE
!! if unpacked berg has more extreme bounds than the same berg in the current PE, update the
!! bounds of the current PE berg to match
subroutine unpack_tabular_buffer_and_update_bounds(lbergs, bberg_data, pebl, buffer, changes)
  integer :: lbergs !< number of local bergs on the current PE
  integer :: bberg_data !< count of berg data being received in the buffer from another PE
  real, pointer :: pebl(:,:) !< array of bergs on the current PE and their lat/lon bounds
  real :: buffer(bberg_data) !< the buffer of bergs being received on the current PE
  integer :: changes !< tracks number of bound changes. Calling loop ends when changes = 0.
  real :: buffer2(bberg_data/5,5)
  integer :: m, n

  !reshape the buffer
  do n=1,bberg_data/5
    buffer2(n,1:5) = buffer((n-1)*5+1:(n-1)*5+5)
  enddo

  do m=1,lbergs
    do n=1,bberg_data/5
      if (pebl(m,1)==buffer2(n,1)) then
        !same berg ID detected in current PE berg list and berg buffer from other PE
        !extend current berg bounds if needed, to reflect more extensive bounds from other PE
        if (pebl(m,2)>buffer2(n,2)) then; pebl(m,2)=buffer2(n,2); changes=changes+1; endif
        if (pebl(m,3)<buffer2(n,3)) then; pebl(m,3)=buffer2(n,3); changes=changes+1; endif
        if (pebl(m,4)>buffer2(n,4)) then; pebl(m,4)=buffer2(n,4); changes=changes+1; endif
        if (pebl(m,5)<buffer2(n,5)) then; pebl(m,5)=buffer2(n,5); changes=changes+1; endif
      endif
    enddo
  enddo

end subroutine unpack_tabular_buffer_and_update_bounds

!> Initializes new iKID bonded icebergs over the lat/lon bounds specified in the tabular iceberg list, and trims
!! excess particles so that the shape of each new tabular iceberg is consistent with
!! its corresponding gridded iceberg mask. This approach guarantees consistent iceberg initialization across PEs
!! Pressure on the ocean is slowly transitioned between ice shelf and iceberg over time, so also returns
!! the fraction to reduce shelf pressure
subroutine ice_shelf_to_bonded_bergs(bergs, TC)
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(tabular_calving_state), pointer :: TC !< A pointer to the tabular calving structure
  integer :: bcount ! number of tabular bergs to initialize on this PE
  type(icebergs_gridded), pointer :: grd
  real :: diameter, dlat, dlon, lon, lat
  real :: minlon, maxlon, minlat, maxlat
  real :: minlon0, maxlon0, minlat0, maxlat0
  real :: minx, miny, cos_lat_ref, dlonscale
  integer :: stderrunit, i, k, nbonds
  logical :: check_bond_quality

  !The number of tabular icebergs to initialize on the current PE
  bcount = TC%berg_pe_count

  !Note: account for the calving mask to be between 0 and 1
  !Initialize bergs over all cells with mask>0. Eliminate a berg if its groundfrac is greater than some threshold
  !(or maybe the grid mask should deal with this), or if the majority of the berg does not overlap a mask>0 cell?
  !Then, calculate each cell's fraction of coverage by particles. If this fraction of coverage is greater
  !than the calving fraction, try eliminating the least-bonded particle in the cell to see if that helps.
  !If this elimination does not improve the match to the calving fraction (considering both the current cell
  !and any other cell that the eliminated particle overlaps), then undo the elimination.
  !Also, overlap between neighboring bergs (newly-or-previously calved) should be considered. Perhaps they can
  !be allowed to interact, but fracture is suppressed until they have fully separated... Or maybe defining an
  !iceberg particle configuration for each individual ice shelf would help in some way? Seems complicated.
  !Or maybe another, adjacent, calving event cannot happen (no pressure) until the first berg drifts away
  !(i.e. no overlap), which also captures backpressure of the existing berg on the emerging berg...

  !TODO: should you consider not calving and grounded ice?

  ! Get the stderr unit number
  stderrunit = stderr()

  grd=>bergs%grd

  if (bcount>0) then
    !Currently all bonded-elements that calve from ice shelves will have the same radius...
    !TODO: Calculate an "adaptive" element size where larger bergs are initialized with larger (and fewer) elements,
    !      thereby reducing the computational expense of iKID. This would not work with
    !      bergs%snap_tabular_calving_to_bonded_grid
    if (bergs%constant_radius_IS_berg>0) then
      diameter = 2*bergs%constant_radius_IS_berg
    else
      call error_mesg('KID, ice_shelf_to_bonded_bergs',&
        'tabular calving from ice shelves currently requires constant_radius_IS_berg>0.!',FATAL)
    endif

    !1) For each tabular berg (TC%berg_list(i,1)), initialize particles between the
    !   min and max longitude and latitude of the berg, (TC%berg_list(i,2:5))

    !calculate the particle spacing
    if (grd%grid_is_latlon) then
      if (bergs%hexagonal_icebergs) then
        dlat=sqrt(3.)*0.5*diameter*(180./pi)/Rearth
      else
        dlat=diameter*(180./pi)/Rearth
      endif
      dlonscale= diameter*(180./pi)/Rearth !used to calculate dlon according to local latitude
    else
      if (bergs%hexagonal_icebergs) then
        dlat=sqrt(3.)*0.5*diameter
      else
        dlat=diameter
      endif
      dlon=diameter
    endif

    !grid bounds
    maxlon0=maxval(grd%lon(          grd%iec, grd%jsc-1:grd%jec))
    maxlat0=maxval(grd%lat(grd%isc-1:grd%iec,           grd%jec))
    minlon0=minval(grd%lon(grd%isc-1        , grd%jsc-1:grd%jec))
    minlat0=minval(grd%lat(grd%isc-1:grd%iec, grd%jsc-1        ))

    !initialize the particles for each new tabular iceberg
    do i = 1,bcount

      minlon=TC%berg_list(i,2); maxlon=min(TC%berg_list(i,3),maxlon0)
      minlat=TC%berg_list(i,4); maxlat=min(TC%berg_list(i,5),maxlat0)

      ! print *,'lon bounds',minlon,maxlon
      ! print *,'lat bounds',minlat,maxlat

      !snap the lat/lon bounds to a constant cartesian grid with spacing equal to the
      !constant diameter of bonded bergs from ice shelves
      !This could be helpful to prevent overlapping when multiple, adjacent bergs are in the process of calving.
      !This will not always work where reentrant, unless grid_is_latlon=.false. and grid spacing is divisible by
      !iceberg diameter
      if (bergs%snap_tabular_calving_to_bonded_grid) then
        if (grd%grid_is_latlon) then
          miny=pi_180*Rearth*minlat !convert minlat to cartesian
          miny=miny-modulo(miny,diameter) !snap miny to the closest, smaller multiple of iceberg diameter
          minlat = miny/(pi_180*Rearth) !update minlat to miny, converted back to latitude

          !calculate minx, using the smallest minx calculated using either minlat or maxlat as reference latitude
          cos_lat_ref=min(cos(maxlat*pi_180),cos(minlat*pi_180))
          minx=pi_180*Rearth*cos_lat_ref*minlon
          minx=minx-modulo(minx,diameter) !snap minx to the closest, smaller multiple of iceberg diameter
          minlon=minx/(pi_180*Rearth*cos_lat_ref) !update minlon to minx, converted back to longitude
        else
          !Snap the cartesian "minlat" and "minlon" to the closest, smaller multiple of iceberg diameter
          minlat=minlat-modulo(minlat,diameter) ; minlon=minlon-modulo(minlon,diameter)
        endif
      endif

      !NOTE: Assuming we do not want to mess with particle shape (i.e. we still assume perfectly circular particles):
      !-To eliminate any potential for particle overlapping between berg conglomerates, which could occur
      ! over reentrant bounds for grid_is_latlon=.true., you could ensure that
      ! particle radius decreases as latitude increases, such that there are always a constant number of bergs in
      ! the x-direction (longitude), where the bergs on either side of the x-reentrant bound align perfectly...
      ! (i.e. a particle grid similar to a spherical grid).
      ! In this case, you could request a target particle diameter at 70 deg lat, e.g. dx_t=3000 m, which is what was
      ! used in the A68a simulation).
      ! Then, defining sigma(lat) = pi_180*Rearth*cos(lat*pi_180), the actual dx at 70 deg lat to prevent overlapping at
      ! reentrant bounds is dx(70) = sigma(70)*360 / n = 3000.5 , where n = floor(360*sigma/dx_t) = # of bergs
      ! at each latitude (for this example, n = 4555).
      ! The berg diameter at each latitude is dx(lat) = sigma(lat) * 360 / n ,
      ! which can then be used to backtrack the new particle "grid".
      ! However, the issue is that for this example, dx(70)=3000.5 m, then dx(80) = 1523.4 m and dx(60) = 4386.5 m.
      ! As particle size decreases, more sub-time steps are needed in the MTS scheme, so this may be unusable.
      ! Especially because the largest shelves (FR and Ross) will have the largest bergs, but too small of particles.
      !-One way around it may be to specify a lower range particle size, e.g. 3000 m, so that if you are below this,
      ! you double the particle size. But this is rather abrupt and awkward.
      !-Another approach is to use a slightly different, but reasonably similar radius, at each latitude to prevent
      ! overlapping, and then have a different bonding pattern. Or any other particle pattern (or lack thereof) with
      ! a different bonding pattern.
      ! While such an approach is not unusual for DEM/bonded-particle models, we have yet to test such an approach
      ! with the model.
      ! So for now, avoid it.
      !-You could also work with a hexagonal grid that tiles perfecly on a sphere. But this does not allow for square
      ! packing.
      ! Only square packing has been tuned. It also may not result in an ideal particle resolution.
      !-Overall, it's probably just easies to work with a constant particle size and constant packing
      ! (square preferred), and simply account for any potential overlapping of particles, which really should not
      ! happen that often anyway.

      lat = minlat
      k=0
      !Start calving bergs on the current PE, keeping in mind that bergs that fall on the eastern and northern
      !boundary of the PE should actually be included in the PEs to the east or north, respectively, instead of
      !the current PE
      do while (lat<=maxlat .and. lat<maxlat0)
        k=k+1

        if (lat>=minlat0) then

          if (grd%grid_is_latlon) dlon = dlonscale/cos(lat*(pi/180.))

          if (bergs%hexagonal_icebergs) then
            lon=minlon-mod(k,2)*0.5*dlon
          else
            lon=minlon
          endif

          do while (lon<=maxlon .and. lon<maxlon0)
            if (lon>=minlon0) then
              !Calve particles that overlap the mask.  Save the overlapping area of
              !each particle with neighboring cells. Their thickness and scaling
              !will be determined below in new_tabular_bergs_thickness_and_pressure.
              call begin_calving_tabular_iceberg_from_shelf(bergs, grd, lon, lat, &
                                                            TC%calve_mask, TC%frac_shelf, TC%h_shelf, 0.5*diameter)
            endif
            lon=lon+dlon
          enddo
        endif
        lat=lat+dlat
      enddo
    enddo
  endif !if bcount>0

  call mpp_update_domains(grd%area, grd%domain)

  !2) Initialize bonds and halo bergs. Eliminate bergs with zero thickness and which are 2 cells away from
  !the ice front. The particles within 2 cells of the front are kept for now, even if they currently have zero
  !thickness, because the pressure on the ocean is slowly transitioned from ice shelf to berg over some period of time;
  !over this period, ice shelf thickness is repeatedly interpolated to the particles, and some of the particles
  !near the front that do not currently have any thickness may eventually receive some thickness over this period
  !as the ice front advects.
  if (bergs%iceberg_bonds_on) then
    !call initialize_iceberg_bonds(bergs, tabular_calving_only=.true.)

    !only includes particles that just initialized
    call update_halo_calved_tabular_icebergs(bergs)

    !bond just the new tabular iceberg particles (which have static_berg=-2),
    !if they are within 1.25*diameter of each other
    !only includes particles that just initialized
    call initialize_iceberg_bonds(bergs, tabular_calving_only=.true.)
    ! if (bergs%mts) then
    ! call transfer_mts_bergs(bergs)
    ! else
    !call update_halo_icebergs(bergs)

    !includes all particles that just initialized
    call connect_all_bonds(bergs, match_bond_pairs=.true.,tabular_calving_only=.true.)
    ! endif
    !includes all particles that just initialized
    if (debug) then
      nbonds=0
      check_bond_quality=.True.
      call count_bonds(bergs, nbonds,check_bond_quality,tabular_calving_only=.true.)
    endif

    !only includes particles that just initialized
    call assign_n_bonds(bergs,tabular_calving_only=.true.)
  endif

  !This can be done with the rest of the bergs?
  !call transfer_mts_bergs(bergs)

  !3)
  !If there are partially-filled grid cells at the ice front (ice shelf CS%hmask==2), then we could end up with multiple
  !rows of partially-filled particles, which we do not want:
  !e.g. In the 1D example below, particles `cc` and `dd` may both end up as partially-filled because they both
  !     overlap the partially-filled grid cell with hmask==2. Instead, this partial fill should be fully given to
  !     'cc' first, and if `cc` becomes filled, any remaining fill can be distributed to `dd`.
  !     This way, we end up with only one edge row of partially-filled particles.
  !Grid cell domains with hmask: |   1   |   1   |   2   |   0   |
  !Labeled particles domains   : |  aa  |  bb  |  cc  |  dd  |
  !i.e. If thickness is smeared over several rows of partially-filled particles near the front, then
  !consolidate it to the particles closest to the fully-filled particles so that partially-filled particles
  !only exist on the outer edge of bonded-particle iceberg conglomerates. Note that this must be recalculated
  !each time step as the ice front may advect over time. Also, interpolate to the grid the fraction that bergs
  !contribute pressure to the ocean surface vs ice shelf pressure (to account for the transition between these
  !two pressures over a specified time scale).
  !Eliminate bergs that are at the end of this transition time and which have zero thickness.

  !For calculating ice shelf pressure on the ocean, the fraction of ice shelf in the cell becomes
  !modified as frac_shelf-frac_cberg, where frac_cberg will be 0 at the start of the transition time between
  !ice shelf and iceberg, and 1 at the end.

  !TODO: ice shelf calving mask could also be fractional over cells that have both calve and no-calve MPs. Could you
  !      treat these cells similarly to the ice front?

  !Calving mask must stay constant until calving is over, which can be detected when a cell obtains a value of
  !frac_cberg_calved>0. Then, eliminate the calving mask there. Probably easiest to allow multiple calving
  !masks on the same PE, but not if they are touching. Simply allow the first calving event to finish before
  !starting the second...

  call new_tabular_bergs_thickness_and_pressure(bergs)

  ! if (bergs%iceberg_bonds_on) then
  !   ! call initialize_iceberg_bonds(bergs, tabular_calving_only=.true.)
  !   ! call update_halo_calved_tabular_icebergs(bergs)

  !   if (bergs%mts) then
  !     call transfer_mts_bergs(bergs)
  !   else
  !     call update_halo_icebergs(bergs)
  !     call connect_all_bonds(bergs, match_bond_pairs=.true.)
  !   endif

  !   nbonds=0
  !   check_bond_quality=.True.
  !   call count_bonds(bergs, nbonds,check_bond_quality)
  !   call assign_n_bonds(bergs)
  ! endif
end subroutine ice_shelf_to_bonded_bergs

!> Calculates thickness new bonded-particle calved from an ice shelf. Also calculates pressure scaling for particles
!! and ice shelf as the calving part of the ice shelf transitions to bonded bergs over time.
subroutine new_tabular_bergs_thickness_and_pressure(bergs)
  type(icebergs), pointer :: bergs !< Container for all types and memory
  ! Local variables
  type(tabular_calving_state), pointer :: TC !< A pointer to the tabular calving structure
  type(icebergs_gridded), pointer :: grd
  type(iceberg), pointer :: berg, other_berg
  type(bond) , pointer :: current_bond
  real, dimension(:,:), pointer :: h_shelf ! The ice shelf thickness field (m)
  real, dimension(:,:), pointer :: frac_shelf ! The fraction of a grid cell covered by the ice shelf [nondim]
  real, dimension(:,:), pointer :: frac_cberg_calved ! Cell fraction of fully-calved bonded bergs from the ice sheet [nondim]
  real, dimension(:,:), pointer :: frac_cberg ! Cell fraction of partially-calved bonded bergs from the ice sheet [nondim]
  integer :: grdi, grdj, c1
  integer :: count, max_count, bcount, bcount_r, bcount_d
  real :: resid_area
  real, allocatable :: pf_area(:,:,:)
  real :: yUxL, yUxC, yUxR
  real :: yCxL, yCxC, yCxR
  real :: yDxL, yDxC, yDxR
  real :: T_scale, min_T_scale

  !A newly-calving iKID conglomerate may have both edge particles (i.e. with empty bond pairs) and
  !interior particles (i.e. with no empty bond pairs) that both overlap partially-filled ice-shelf
  !cells (where 0<frac_shelf<1). Consequently, these particles are only partially filled
  !(defined as having berg%mass_scaling<1) during the interpolation of gridded ice shelf fields to the
  !particles. However, only edge particles should be partially-filled, so here, some mass is transferred
  !from the edge particles to the interior particles to fill the interior particles completely.
  grd=>bergs%grd
  grd%frac_cberg_calved(:,:,:)=0.
  grd%frac_cberg(:,:,:)=0.

  TC=>bergs%TC
  h_shelf=>TC%h_shelf
  frac_shelf=>TC%frac_shelf
  frac_cberg_calved=>TC%frac_cberg_calved
  frac_cberg=>TC%frac_cberg

  !now that all bergs are initialized and bonded, all static_berg statuses should be positive
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg))
      berg%static_berg=abs(berg%static_berg)
      berg%sss=0
      berg=>berg%next
    enddo
  enddo; enddo

  !1) Determine how many bonds away from an initially "full" particle (completely overlaps filled ice shelf cells)
  !   each initially partially-full particle is. Partial fill particles closer to a full particle will preferentially
  !   be filled with ice shelf mass first, so that they can actually end up converting to full particles and only
  !   the outermost edge of a conglomerate will contain partially-full particles.
  count=1 !number of bonds a partially-full particle is away from a full particle
  max_count=0
  T_scale=0

  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied ! only process conglomerates overlapping the comp domain
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      !Start from a full particle (static_berg==2)
      !Processed bergs have their IDs set negative
      if (berg%static_berg==2 .and. berg%id>0) then
        bergs%new_tabular_list%first=>null()
        !returns a list of partially-full particles that are directly connected to full particles (i.e. count==1)
        berg%id=-berg%id
        call make_list_of_bonded_to_full(berg,bergs%new_tabular_list%first)
        if (associated(bergs%new_tabular_list%first)) then
          !For partially-full particles that are not directly bonded to a full particle
          !(i.e. count>1), determine how many bonds they are away from a full particle
          count=2
          call assign_bonds_from_full(bergs%new_tabular_list%first,count)
          max_count=max(count-1,max_count)
        endif
      endif
      berg=>berg%next
    enddo
  enddo; enddo
  !Also process conglomerates without any initially "full" particles. In this case, any partial-fill
  !particle that overlaps a full ice shelf cell will be given count==1 (as if it is one bond away from an
  !initially filled particle.
  do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg)) ! loop over all bergs
      if ((berg%static_berg==2.5) .and. berg%id>0) then
        bergs%new_tabular_list%first=>null()
        !this berg is partially-filled and overlaps a full ice shelf cell,
        !but does not eventually connect to a full particle
        !We can treat it as if it is bonded to a full particle
        berg%id=-berg%id
        count=2 !to make sure that max_count will be >=1
        call make_list_of_bonded_to_full2(berg,bergs%new_tabular_list%first)
        if (associated(bergs%new_tabular_list%first)) then
          call assign_bonds_from_full(bergs%new_tabular_list%first,count)
          max_count=max(count-1,max_count)
        endif
      endif
      berg=>berg%next
    enddo
  enddo; enddo

  call mpp_max(max_count)

  if (max_count==0) then
    !reset all berg ids
    !(if max_count>0, this is done elsewhere below)
    do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg))
        berg%id=abs(berg%id)
        berg=>berg%next
      enddo
    enddo; enddo
  endif

  if (max_count>0) then

    !2) use the gridded field "pf_area" field to calculate the total area of
    !   (initially) partially-full particles in each cell that are a
    !   certain "count" of bounds away from a full particle

   ! max_count=max_count+1
    allocate(pf_area(grd%isd:grd%ied,grd%jsd:grd%jed,max_count), source=0.0)
    do count=1,max_count
      grd%pf_area(:,:,:)=0.
      do grdj = grd%jsd,grd%jed ; do grdi = grd%isd,grd%ied
        berg=>bergs%list(grdi,grdj)%first
        do while (associated(berg))

          berg%id=abs(berg%id)

          if (grdj >= grd%jsc-1 .and. grdj <= grd%jec+1 .and. grdi >= grd%isc-1 .and. grdi <= grd%iec+1 .and. &
              berg%static_berg<=3 .and. berg%static_berg>2 .and. int(berg%sss)==count .and. berg%halo_berg<=1) then

            !add berg area in cell to pf_area
            yUxL = berg%sst
            yUxC = berg%uo
            yUxR = berg%vo
            yCxL = berg%ui
            yCxC = berg%vi
            yCxR = berg%ua
            yDxL = berg%va
            yDxC = berg%ssh_x
            yDxR = berg%ssh_y

            call spread_variable_across_cells(grd, grd%pf_area, berg%length * berg%width, grdi, grdj, &
                                              yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, 1.0)
          endif
          berg=>berg%next
        enddo
      enddo; enddo
      call sum_up_spread_fields(bergs, pf_area(grd%isc:grd%iec,grd%jsc:grd%jec,count), 'pf_area', ignore_mask_in=.true.)
      !save for future diagnostics? Size of TC%saved_pf_area is (grd%isd:grd%ied,grd%jsd:grd%jed,10)
      if (count<=10 .and. grd%id_pf_area>0) &
        TC%saved_pf_area(grd%isc:grd%iec,grd%jsc:grd%jec,count)=pf_area(grd%isc:grd%iec,grd%jsc:grd%jec,count)
    enddo

    !3) Convert pf_area from representing the area of partially-full particles with various "counts"
    !   overlapping the cell, to the scaling (between 0 and 1) used for the cell when interpolating
    !   grid thickness to the particles and determining mass scaling of the particles.
    !   There is a separate scaling for each "count" category.
    do grdj = grd%jsc-1,grd%jec+1 ; do grdi = grd%isc-1,grd%iec+1
      !If you wanted to use partially-masked cells, you would need to make sure the mask is retained until the particles
      !are released, and then multiply resid_area by the mask for each cell. But simpler for now to only initialize particles
      !that overlap a fully-masked cells.
      resid_area = frac_shelf(grdi,grdj) * grd%area_um(grdi,grdj)
      do count=1,max_count
        if (pf_area(grdi,grdj,count)>0) then
          if (pf_area(grdi,grdj,count)<resid_area) then
            !There is more cell area than area within the cell from overlapping particles with the current "count"
            !These particles will have a scaling factor of 1 from this cell
            resid_area=resid_area-pf_area(grdi,grdj,count)
            pf_area(grdi,grdj,count)=1
          else
            !The area within the cell from overlapping particles with the current "count" >= the remaining cell area.
            !These particles will have a scaling factor <= 1 from this cell
            pf_area(grdi,grdj,count)=resid_area/pf_area(grdi,grdj,count)
            if (count<max_count) pf_area(grdi,grdj,(count+1):max_count)=0
            resid_area=0
            exit
          endif
        endif
      enddo
    enddo; enddo
  endif !if max_count>0

  bcount=0
  bcount_r=0
  bcount_d=0
  min_T_scale=1
  T_scale=0
  !4) Using the scalings saved on pf_area, calculate thickness and mass scaling on the partially-full particles.
  !!  Interpolate berg thickness and time-scaling (percent berg pressure vs ice-shelf pressure) to grid.
  do grdj = grd%jsc-1,grd%jec+1 ; do grdi = grd%isc-1,grd%iec+1
    berg=>bergs%list(grdi,grdj)%first
    do while (associated(berg))
      if (berg%static_berg<2 .or. berg%halo_berg>1) then
        !non-calving berg or, if MTS, a berg that lies outside the current PE's grid, but is part of a conglomerate that
        !overlaps the current PE
        berg=>berg%next
      else
        bcount=bcount+1

        yUxL = berg%sst
        yUxC = berg%uo
        yUxR = berg%vo
        yCxL = berg%ui
        yCxC = berg%vi
        yCxR = berg%ua
        yDxL = berg%va
        yDxC = berg%ssh_x
        yDxR = berg%ssh_y

        if (grd%parity_x(grdi,grdj)<0.) then
          c1=-1
        else
          c1=1
        endif

        if (berg%static_berg<=3) then !Fully-filled bergs, or bergs that overlap a full cell:

          if (berg%static_berg==2) then

            berg%mass_scaling = 1
            berg%thickness =   yCxC*h_shelf(grdi   ,grdj   ) + &
                            (((yUxL*h_shelf(grdi-c1,grdj+c1) + yDxR*h_shelf(grdi+c1,grdj-c1))  + &
                              (yUxR*h_shelf(grdi+c1,grdj+c1) + yDxL*h_shelf(grdi-c1,grdj-c1))) + &
                             ((yUxC*h_shelf(grdi   ,grdj+c1) + yDxC*h_shelf(grdi   ,grdj-c1))  + &
                              (yCxL*h_shelf(grdi-c1,grdj   ) + yCxR*h_shelf(grdi+c1,grdj   ))))

          elseif (berg%static_berg>2) then !Initially partially-full bergs:
            count=int(berg%sss)
            berg%mass_scaling =   yCxC*pf_area(grdi   ,grdj   ,count) + &
                               (((yUxL*pf_area(grdi-c1,grdj+c1,count) + yDxR*pf_area(grdi+c1,grdj-c1,count))  + &
                                 (yUxR*pf_area(grdi+c1,grdj+c1,count) + yDxL*pf_area(grdi-c1,grdj-c1,count))) + &
                                ((yUxC*pf_area(grdi   ,grdj+c1,count) + yDxC*pf_area(grdi   ,grdj-c1,count))  + &
                                 (yCxL*pf_area(grdi-c1,grdj   ,count) + yCxR*pf_area(grdi+c1,grdj   ,count))))

            if (berg%mass_scaling>0) then
              berg%thickness = (  yCxC*pf_area(grdi   ,grdj   ,count)*h_shelf(grdi   ,grdj   )   + &
                               (((yUxL*pf_area(grdi-c1,grdj+c1,count)*h_shelf(grdi-c1,grdj+c1)   + &
                                  yDxR*pf_area(grdi+c1,grdj-c1,count)*h_shelf(grdi+c1,grdj-c1))  + &
                                 (yUxR*pf_area(grdi+c1,grdj+c1,count)*h_shelf(grdi+c1,grdj+c1)   + &
                                  yDxL*pf_area(grdi-c1,grdj-c1,count)*h_shelf(grdi-c1,grdj-c1))) + &
                                ((yUxC*pf_area(grdi   ,grdj+c1,count)*h_shelf(grdi   ,grdj+c1)   + &
                                  yDxC*pf_area(grdi   ,grdj-c1,count)*h_shelf(grdi   ,grdj-c1))  + &
                                 (yCxL*pf_area(grdi-c1,grdj   ,count)*h_shelf(grdi-c1,grdj   )   + &
                                  yCxR*pf_area(grdi+c1,grdj   ,count)*h_shelf(grdi+c1,grdj   ))))) / berg%mass_scaling
            else
              berg%thickness=0
            endif
          endif

          !--All bergs--

          berg%mass = berg%width * berg%length * berg%thickness * bergs%rho_bergs

          !The time-based pressure scaling factor. Over a timescale (hours) of bergs%shelf_to_tabular_hrs
          !(using the icebergs module "yearday" time convention), transition smoothly between:
          !  berg_scaling=0 for 0%   berg pressure on ocean and 100% ice shelf pressure
          !  berg_scaling=1 for 100% berg pressure on ocean and 0%   ice shelf pressure
          T_scale = min(((bergs%current_year*366.+bergs%current_yearday)-&
                         (berg%start_year*366.+berg%start_day))*24./bergs%shelf_to_tabular_hours, 1.0)

          !Interpolate the T_scale to the grid to modify the ice shelf pressure felt on
          !the ocean. This interpolation accounts for the possibility of multiple tabular bergs with different
          !T_scale that contribute to the same cell

          !(1) Process "calving bergs" with T_scale==1 (these bergs have now fully-transitioned to non-static,
          !    fully-calved bergs):
          !--interp to grid only from calving bergs (cberg) with T_scale==1--
          !frac_cberg_calved = sum(INTERP(cberg%area * cberg%mass_scaling))/cell_area !fraction of cell pressure
          !from calving bergs
          !--Then, permanently adjust frac_shelf as--:
          !frac_shelf = max(frac_shelf - frac_cberg_calved,0).
          !--If frac_shelf == 0, adjust ice thickness, hmask, etc accordingly.
          !(2) Process calving bergs with T_scale<1:
          !--interp to the grid only from calving bergs (cberg) without T_scale==1 (these bergs have not fully-calved yet)--
          !cberg%mass_scaling=cberg%mass_scaling * T_scale
          !frac_cberg        = sum(INTERP(cberg%area * cberg%mass_scaling))/cell_area !fraction of cell pressure
          !from calving bergs
          !--Then, implement frac_shelf as--:
          !frac_shelf_new = max(frac_shelf - frac_cberg,0)
          !(3) When interpolation time is up, the calving mask and associated ice shelf needs to be eliminated

          min_T_scale=min(T_scale,min_T_scale)

          if (T_scale==1) then !berg is old enough to now evolve as a dynamic berg, but will be deleted if massless
            if (berg%mass_scaling==0) then
              if (bergs%mts) then
                !mark the berg for deletion after halo transfers
                berg%static_berg=0.
                berg%mass_scaling=-1
                berg=>berg%next
              else
                !remove massless particle
                other_berg=>berg
                berg=>berg%next
                call delete_all_bonds(other_berg)
                call delete_iceberg_from_list(bergs%list(grdi,grdj)%first,other_berg)
              endif

              bcount_d=bcount_d+1
            else
              !Particle has fully calved from the ice shelf, and will evolve as a dynamic iceberg
              !The section of the ice shelf from where the particle calved will be eliminated
              bcount_r=bcount_r+1
              berg%static_berg=0.
              !Removes edge particles from a calving iceberg conglomerate so that it can more easily flow away from the ice shelf
              if (bergs%remove_tabular_outer_bonds_when_calve .and. berg%n_bonds<bergs%max_bonds) &
                berg%static_berg=10 !call delete_all_bonds(berg)

              call spread_variable_across_cells(grd, grd%frac_cberg_calved, berg%length * berg%width * berg%mass_scaling, &
                                                grdi, grdj, yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, 1.0)
              berg=>berg%next
            endif

          elseif (T_scale>=0) then !berg is not yet old enough to be released as a dynamic iceberg.
            !Keep gradually increasing berg pressure while decreasing ice-shelf pressure to ensure a smooth transition

            berg%mass_scaling=berg%mass_scaling*T_scale
            call spread_variable_across_cells(grd, grd%frac_cberg, berg%length * berg%width * berg%mass_scaling, &
                                              grdi, grdj, yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, 1.0)
            berg=>berg%next

          elseif (T_scale.lt.0) then
            call error_mesg('KID, new_tabular_bergs_thickness_and_pressure','Berg pressure scaling is negative!', FATAL)
          endif
        endif !end if (berg%static_berg<=3)
      endif ! end if (berg%static_berg<2)
    enddo !do while associated berg
  enddo; enddo

  !You need to delete particles from the list after transferring them between PEs again, so that the other PEs can unbond from
  !them> So maybe you can do a transfer of only those bergs slated to be deleted. Or just transfer the id's and grid cells
  !of the bergs that must be deleted. Or does this work itself out with paralellization anyway?


  ! if (mpp_pe().eq.mpp_root_pe()) print *,'Tscale',T_scale
  if (T_scale/=0) print *,'pe',mpp_pe(),'T_scale',T_scale,'min_T_scale',min_T_scale,'bcount',bcount,'bcount_r',bcount_r,'bcount_d',bcount_d
  call sum_up_spread_fields(bergs, frac_cberg_calved(grd%isc:grd%iec,grd%jsc:grd%jec), 'frac_cberg_calved', ignore_mask_in=.true.)
  call sum_up_spread_fields(bergs, frac_cberg(grd%isc:grd%iec,grd%jsc:grd%jec)       , 'frac_cberg'       , ignore_mask_in=.true.)

  !Adjust frac_cberg_calved and the iceberg mask
  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec

    if (bergs%remove_tabular_outer_bonds_when_calve) then
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg))
        if (berg%static_berg==10) then
          berg%static_berg=0
          call delete_all_bonds(berg)
        else
          berg=>berg%next
        endif
      enddo
    endif
    !In the ice shelf code, cells with frac_cberg_calved == frac_shelf will cause all ice shelf in the cell to be eliminated.
    !Alternatively, all ice shelf in the cell will also be eliminated if frac_cberg_calved == 1.
    !Here, account for potential round-off error so that frac_cberg_calved definitely equals frac_shelf (or 1) where needed -- this should
    !occur wherever the calve mask fraction (saved on the cell's particles hi field) equals 1.
    if (frac_cberg(grdi,grdj)>1) frac_cberg(grdi,grdj)=1
    if (frac_cberg_calved(grdi,grdj)>1) frac_cberg_calved(grdi,grdj)=1
    if (frac_cberg_calved(grdi,grdj)>0 .and. &
      (frac_cberg_calved(grdi,grdj)<frac_shelf(grdi,grdj) .or. frac_shelf(grdi,grdj)/=1)) then
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg))
        if (berg%hi>=1) then
          !All ice in this cell should calve fully. Either set the frac_cberg_calved to frac_shelf, or simply set to 1,
          !Setting it to 1 is a bit more helpful in diagnostic output to easily differentiate between full (1) and partial (<1)
          !calve cells. While it is possible (though unlikely) that this statement could be triggered by a non-calving (already released)
          !berg with hi==1, this is not an issue because any non-zero frac_cberg_calved at the ice front should equal 1 anyway.
          frac_cberg_calved(grdi,grdj)=1
          berg=>null()
        else
          berg=>berg%next
        endif
      enddo
    endif
    !Immediately adjust the iceberg mask to account for fully-calved icebergs
    !Even if the cell is still partially-full of ice shelf after calving, we still unmask it, as any neighboring masked cell
    !will push bergs away using the coastal_drift and tidal_drift features
    ! if (frac_cberg_calved(grdi,grdj)>0) &
    !    grd%msk(grdi,grdj)=1.0-max(frac_shelf(grdi,grdj) - frac_cberg_calved(grdi,grdj),0.)
    if (frac_cberg_calved(grdi,grdj)>0) grd%msk(grdi,grdj)=1.
  enddo; enddo
  call mpp_update_domains(frac_cberg_calved, grd%domain, complete=.false.)
  call mpp_update_domains(frac_cberg,        grd%domain, complete=.false.)
  call mpp_update_domains(grd%msk,           grd%domain, complete=.true.)
  if (allocated(pf_area)) deallocate(pf_area)

  if (.not. bergs%mts) then
    !if bergs%mts, then this is done immediately following process_tabular_calving, within
    !interp_gridded_fields_to_bergs, so no need to do it here...
    do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
      berg=>bergs%list(grdi,grdj)%first
      do while (associated(berg))
        berg%mask_status=grd%msk(grdi,grdj)
        berg=>berg%next
      enddo
    enddo; enddo
  endif
end subroutine new_tabular_bergs_thickness_and_pressure

!> Initialize (begin calving) a tabular iceberg particle from an ice shelf at the given lat/lon coordinates.
!! Save its overlapping area with neighboring cells, which will be used to determine its
!! thickness, mass, and mass scaling in subroutine new_tabular_berg_thickness_and_pressure
!! Interpolation of external fields to the new particle will occur after the berg is released and no longer static
subroutine begin_calving_tabular_iceberg_from_shelf(bergs, grd, lon, lat, calve_mask, frac_shelf, h_shelf, radius)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(icebergs_gridded), pointer :: grd
  real :: lon !< longitude of the new iceberg
  real :: lat !< latitude of the new iceberg
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: calve_mask !< ice shelf calving mask
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: frac_shelf !< The fraction of a grid cell covered by
                                                                             !! the ice shelf [nondim].
  real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: h_shelf !< The ice shelf thickness field (m)
  real :: radius !< radius of the new iceberg
  ! Local variables
  integer :: i,j,k,icnt,icntmax
  real :: orientation
  type(iceberg) :: newberg
  logical :: lret, lres
  real :: xi, yj, calving_to_bergs, calved_to_berg, heat_to_bergs, heat_to_berg
  integer :: stderrunit
  real, pointer :: mass_scaling, initial_thickness, initial_width, initial_length
  logical :: allocations_done
  real :: rx,ry,yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR,c1
  real :: pmask, width
  real :: x,y
  logical :: correct_pmask, overlaps_ocean
  real, dimension(3,3) :: cm_arr, fs_arr
  integer, dimension(3,3) :: overlaps_arr
  real, pointer :: yUxL_overlap, yUxC_overlap, yUxR_overlap
  real, pointer :: yCxL_overlap, yCxC_overlap, yCxR_overlap
  real, pointer :: yDxL_overlap, yDxC_overlap, yDxR_overlap
  ! real, parameter :: rho_seawater=1035.

  ! Get the stderr unit number
  stderrunit = stderr()

  rx = 0.; ry = 0.

!  grd%real_calving(:,:,:)=0.
!  calving_to_bergs=0.
!  heat_to_bergs=0.
  icntmax=0

  ! allocations_done=.false.

  lres=find_cell(grd, lon, lat, i, j)

  if (.not. lres) return
  !Assume that the calve mask spans over the cells that shelf associated with the berg may advect into
  !Get rid of bergs that clearly do not overlap the calve mask
  if (all(frac_shelf(i-1:i+1,j-1:j+1)==0) .and. all(calve_mask(i-1:i+1,j-1:j+1)==0)) return

  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  if (.not.lret) then
    write(stderrunit,*) 'KID, calve_icebergs: something went very wrong!',i,j,xi,yj
    call error_mesg('KID, calve_icebergs', 'berg is not in the correct cell!', FATAL)
  endif
  if (debug.and.(xi<0..or.xi>1..or.yj<0..or.yj>1.)) then
    write(stderrunit,*) 'KID, calve_icebergs: something went very wrong!',i,j,xi,yj
    call error_mesg('KID, calve_icebergs', 'berg xi,yj is not correct!', FATAL)
  endif

  !Ignore bergs on the N and E boundary of the PE, as they will be included in the PEs to the N or E, respectively
  !But the find_cell call should not allow bergs to be found on these boundaries, anyway.
  if ((i==grd%iec .and. xi==1) .or. (j==grd%jec .and. yj==1)) return

  !Do not calve from grounded cells?
  ! if ((bergs%rho_bergs/rho_seawater)*h_shelf(i,j)>grd%ocean_depth(i,j)) return

  ! if (grd%msk(i,j)<0.5) then
  !   write(stderrunit,*) 'KID, calve_icebergs: WARNING!!! Iceberg born in land cell',i,j,newberg%lon,newberg%lat
  !   if (debug) call error_mesg('KID, calve_icebergs', 'Iceberg born in Land Cell!', FATAL)
  ! endif

  if (bergs%hexagonal_icebergs) then
    width=sqrt((radius**2)*2*sqrt(3.))
  else
    width=2*radius
  endif

  ! !interpolate gridded variables to new iceberg
  ! if (grd%tidal_drift>0.) then
  !   call getRandomNumbers(rns, rx)
  !   call getRandomNumbers(rns, ry)
  !   rx = 2.*rx - 1.; ry = 2.*ry - 1.
  ! endif

  call calving_tabular_particle_grid_overlap(bergs, width*width, i, j, xi, yj, &
                                             yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR)

  if (grd%parity_x(i,j)<0.) then
    c1=-1
  else
    c1=1
  endif

  !forget about this particle if it does not overlap any masked cells
  pmask =   yCxC*calve_mask(i   ,j   ) + &
          ((yUxL*calve_mask(i-c1,j+c1) + yDxR*calve_mask(i+c1,j-c1))  + &
           (yUxR*calve_mask(i+c1,j+c1) + yDxL*calve_mask(i-c1,j-c1))) + &
          ((yUxC*calve_mask(i   ,j+c1) + yDxC*calve_mask(i   ,j-c1))  + &
           (yCxL*calve_mask(i-c1,j   ) + yCxR*calve_mask(i+c1,j   )))

  if (pmask<=0) return

  grd%area(i,j)=grd%area_um(i,j)

  cm_arr(1,1)=calve_mask(i-c1,j-c1); cm_arr(1,2)=calve_mask(i-c1,j); cm_arr(1,3)=calve_mask(i-c1,j+c1)
  cm_arr(2,1)=calve_mask(i   ,j-c1); cm_arr(2,2)=calve_mask(i   ,j); cm_arr(2,3)=calve_mask(i   ,j+c1)
  cm_arr(3,1)=calve_mask(i+c1,j-c1); cm_arr(3,2)=calve_mask(i+c1,j); cm_arr(3,3)=calve_mask(i+c1,j+c1)

  fs_arr(1,1)=frac_shelf(i-c1,j-c1); fs_arr(1,2)=frac_shelf(i-c1,j); fs_arr(1,3)=frac_shelf(i-c1,j+c1)
  fs_arr(2,1)=frac_shelf(i   ,j-c1); fs_arr(2,2)=frac_shelf(i   ,j); fs_arr(2,3)=frac_shelf(i   ,j+c1)
  fs_arr(3,1)=frac_shelf(i+c1,j-c1); fs_arr(3,2)=frac_shelf(i+c1,j); fs_arr(3,3)=frac_shelf(i+c1,j+c1)

  overlaps_arr(:,:)=0
  if (yDxL>0) overlaps_arr(1,1)=1; if (yCxL>0) overlaps_arr(1,2)=1; if (yUxL>0) overlaps_arr(1,3)=1
  if (yDxC>0) overlaps_arr(2,1)=1; if (yCxC>0) overlaps_arr(2,2)=1; if (yUxC>0) overlaps_arr(2,3)=1
  if (yDxR>0) overlaps_arr(3,1)=1; if (yCxR>0) overlaps_arr(3,2)=1; if (yUxR>0) overlaps_arr(3,3)=1

  !Fix round-off error that may cause some bergs to erroneously have pmask/=0
  if (pmask/=1 .and. pmask>0.99) then

    if (sum(overlaps_arr*cm_arr)==(sum(overlaps_arr))) pmask=1

    ! correct_pmask=.true.
    ! if (yCxC>0 .and. calve_mask(i   ,j   )/=1) correct_pmask=.false.
    ! if (yUxL>0 .and. calve_mask(i-c1,j+c1)/=1) correct_pmask=.false.
    ! if (yDxR>0 .and. calve_mask(i+c1,j-c1)/=1) correct_pmask=.false.
    ! if (yUxR>0 .and. calve_mask(i+c1,j+c1)/=1) correct_pmask=.false.
    ! if (yDxL>0 .and. calve_mask(i-c1,j-c1)/=1) correct_pmask=.false.
    ! if (yUxC>0 .and. calve_mask(i   ,j+c1)/=1) correct_pmask=.false.
    ! if (yDxC>0 .and. calve_mask(i   ,j-c1)/=1) correct_pmask=.false.
    ! if (yCxL>0 .and. calve_mask(i-c1,j   )/=1) correct_pmask=.false.
    ! if (yCxR>0 .and. calve_mask(i+c1,j   )/=1) correct_pmask=.false.
    ! if (correct_pmask) pmask=1
  endif

  !for debugging:
  newberg%cn=pmask

  !temporarily save on hi the fraction that the particle's cell is calve_masked.
  !Typically, this will equal the calve mask ([0,1]), calve mask could also be >1 (to differentiate
  !between multiple bergs: left of decimal = berg number, right of decimal = calve mask fraction), so
  !always adjust so that hi is between 0 and 1.
  newberg%hi=calve_mask(i,j)
  if (calve_mask(i,j)>1) then
    newberg%hi=calve_mask(i,j)-floor(calve_mask(i,j))
  else
    newberg%hi=calve_mask(i,j)
  endif

  !Full particles get static_berg=-2
  !Partially-full particles that overlap a full cell get static_berg=-2.5
  !Otherwise, particle does overlaps only a non-full but masked cell (static_berg=-3)
  !In all cases, after these bergs receive their bonds, they get static_berg=abs(static_berg)
  newberg%static_berg=-2

  overlaps_ocean=.false.
  if (pmask>1) then
    call error_mesg('KID, begin_calving_tabular_iceberg_from_shelf',&
                    'pmask is somehow greater than 1!', FATAL)
  elseif (pmask==1) then
    !if pmask is 1 (berg only overlaps fully-masked cells), but the berg overlaps ocean (rather than ice shelf),
    !it will be processed as a partially-full or non-full cell
    if (sum(overlaps_arr*fs_arr)/=(sum(overlaps_arr))) overlaps_ocean=.true.
    ! if (yCxC>0 .and. frac_shelf(i   ,j   )<1) overlaps_ocean=.true.
    ! if (yUxL>0 .and. frac_shelf(i-c1,j+c1)<1) overlaps_ocean=.true.
    ! if (yDxR>0 .and. frac_shelf(i+c1,j-c1)<1) overlaps_ocean=.true.
    ! if (yUxR>0 .and. frac_shelf(i+c1,j+c1)<1) overlaps_ocean=.true.
    ! if (yDxL>0 .and. frac_shelf(i-c1,j-c1)<1) overlaps_ocean=.true.
    ! if (yUxC>0 .and. frac_shelf(i   ,j+c1)<1) overlaps_ocean=.true.
    ! if (yDxC>0 .and. frac_shelf(i   ,j-c1)<1) overlaps_ocean=.true.
    ! if (yCxL>0 .and. frac_shelf(i-c1,j   )<1) overlaps_ocean=.true.
    ! if (yCxR>0 .and. frac_shelf(i+c1,j   )<1) overlaps_ocean=.true.
  endif

  if (pmask<1 .or. overlaps_ocean) then
    if (any(overlaps_arr + ceiling(cm_arr) + fs_arr == 3)) then

    ! if ((yCxC>0 .and. calve_mask(i   ,j   )>0 .and. frac_shelf(i   ,j   )==1) .or. &
    !     (yUxL>0 .and. calve_mask(i-c1,j+c1)>0 .and. frac_shelf(i-c1,j+c1)==1) .or. &
    !     (yDxR>0 .and. calve_mask(i+c1,j-c1)>0 .and. frac_shelf(i+c1,j-c1)==1) .or. &
    !     (yUxR>0 .and. calve_mask(i+c1,j+c1)>0 .and. frac_shelf(i+c1,j+c1)==1) .or. &
    !     (yDxL>0 .and. calve_mask(i-c1,j-c1)>0 .and. frac_shelf(i-c1,j-c1)==1) .or. &
    !     (yUxC>0 .and. calve_mask(i   ,j+c1)>0 .and. frac_shelf(i   ,j+c1)==1) .or. &
    !     (yDxC>0 .and. calve_mask(i   ,j-c1)>0 .and. frac_shelf(i   ,j-c1)==1) .or. &
    !     (yCxL>0 .and. calve_mask(i-c1,j   )>0 .and. frac_shelf(i-c1,j   )==1) .or. &
    !     (yCxR>0 .and. calve_mask(i+c1,j   )>0 .and. frac_shelf(i+c1,j   )==1)) then

      newberg%static_berg=-2.5
    else
      newberg%static_berg=-3
    endif
  endif

  !otherwise, save its overlap with neighboring cells, making use of some otherwise unneeded icenewberg memory
  newberg%sst   = yUxL
  newberg%uo    = yUxC
  newberg%vo    = yUxR
  newberg%ui    = yCxL
  newberg%vi    = yCxC
  newberg%ua    = yCxR
  newberg%va    = yDxL
  newberg%ssh_x = yDxC
  newberg%ssh_y = yDxR

  !TODO within this routine, save a list of soon-to-calve tabular bergs for each grid cell, with their area in the cell?
  !or save area in a cell on the particle...
  ! call spread_grid_var_to_particle(bergs, newberg, h_shelf, i, j, xi, yj, newberg%thickness, var_frac=frac_shelf)

  !Do not calve grounded particles?
  !if ((bergs%rho_bergs/rho_seawater)*berg%thickness>newberg%od) return

  newberg%lon=lon; newberg%lat=lat
  newberg%ine=i;   newberg%jne=j
  newberg%xi=xi;   newberg%yj=yj
  newberg%width=width
  newberg%length=newberg%width

  newberg%uvel=0.; newberg%vvel=0.
  if (bergs%interactive_icebergs_on .or. footloose) then
    newberg%uvel_prev=0.;        newberg%vvel_prev=0.
    newberg%uvel_old=0.;         newberg%vvel_old=0.
    newberg%lon_old=newberg%lon; newberg%lat_old=newberg%lat
  endif
  newberg%fl_k=0.
  newberg%axn=0.; newberg%ayn=0.
  newberg%bxn=0.; newberg%byn=0.

  newberg%start_lon=newberg%lon
  newberg%start_lat=newberg%lat
  newberg%start_year=bergs%current_year
  newberg%id = generate_id(grd, i, j)
  newberg%start_day=bergs%current_yearday
!  newberg%start_mass=initial_mass
!  newberg%mass_scaling=mass_scaling
  newberg%mass_of_bits=0.
  newberg%mass_of_fl_bits=0.
  newberg%mass_of_fl_bergy_bits=0.
  newberg%halo_berg=0.

  newberg%heat_density=grd%stored_heat(i,j)/grd%stored_ice(i,j,k) ! This is in J/kg

  if (bergs%mts) then
    if (.not. allocations_done) then
      if (.not. allocated(newberg%axn_fast)) allocate(newberg%axn_fast)
      if (.not. allocated(newberg%ayn_fast)) allocate(newberg%ayn_fast)
      if (.not. allocated(newberg%bxn_fast)) allocate(newberg%bxn_fast)
      if (.not. allocated(newberg%byn_fast)) allocate(newberg%byn_fast)
      if (.not. allocated(newberg%conglom_id)) allocate(newberg%conglom_id)
    endif
    newberg%axn_fast=0.; newberg%ayn_fast=0.; newberg%bxn_fast=0.; newberg%byn_fast=0.; newberg%conglom_id=0
  endif

  if (bergs%iceberg_bonds_on) then
    if (.not. allocations_done) then
      if (.not. allocated(newberg%n_bonds))  allocate(newberg%n_bonds)
    endif
    newberg%n_bonds=0
  endif

  if (bergs%dem) then
    if (.not. allocations_done) then
      if (.not. allocated(newberg%ang_vel)) allocate(newberg%ang_vel)
      if (.not. allocated(newberg%ang_accel)) allocate(newberg%ang_accel)
      if (.not. allocated(newberg%rot)) allocate(newberg%rot)
    endif
    newberg%ang_vel=0.; newberg%ang_accel=0.; newberg%rot=0.
  endif

  if (bergs%tabular_calving) then !(obviously this is true)
    if (.not. allocated(newberg%mask_status)) allocate(newberg%mask_status)
    newberg%mask_status=grd%msk(i,j)
  endif

  call add_new_berg_to_list(bergs%list(i,j)%first, newberg)
  ! calved_to_berg=initial_mass*mass_scaling ! Units of kg
  ! Heat content TODO
  ! heat_to_berg=calved_to_berg*newberg%heat_density ! Units of J
  ! grd%stored_heat(i,j)=grd%stored_heat(i,j)-heat_to_berg
  ! heat_to_bergs=heat_to_bergs+heat_to_berg
  ! ! Stored mass
  ! grd%stored_ice(i,j,k)=grd%stored_ice(i,j,k)-calved_to_berg
  ! calving_to_bergs=calving_to_bergs+calved_to_berg
  ! grd%real_calving(i,j,k)=grd%real_calving(i,j,k)+calved_to_berg/bergs%dt

  bergs%nbergs_calved=bergs%nbergs_calved+1

  ! allocations_done=.true.

!  bergs%net_calving_to_bergs=bergs%net_calving_to_bergs+calving_to_bergs
!  bergs%net_heat_to_bergs=bergs%net_heat_to_bergs+heat_to_bergs

end subroutine begin_calving_tabular_iceberg_from_shelf

!> For calving tabular bergs from the ice shelf. Save the area of overlap the particles have with each surrounding cell.
!! If the particle does not overlap
subroutine calving_tabular_particle_grid_overlap(bergs, Area, i, j, x, y, &
                                                 yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  real :: area !< Area of the iceberg that is calving
  integer, intent(in) :: i !< i-index of cell contained center of berg
  integer, intent(in) :: j !< j-index of cell contained center of berg
  real, intent(in) :: x !< Nondimensional x-position within cell [0,1]
  real, intent(in) :: y !< Nondimensional y-position within cell [0,1]

  ! Local variables
  type(icebergs_gridded), pointer :: grd
  real :: xL, xC, xR, yD, yC, yU, L
  real :: yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR
  real :: S, H, origin_x, origin_y, x0, y0
  real :: Area_Q1,Area_Q2 , Area_Q3,Area_Q4, Area_hex, Area_square
  real :: tol, orientation
  real :: Dn, Hocean
  ! real, parameter :: rho_seawater=1035.
  integer :: stderrunit
  logical :: zero_fill

  ! Get the stderr unit number
  stderrunit = stderr()

!  tol=1.e-10
  grd=>bergs%grd

  !Initialize weights for each cell
  yDxL=0.  ; yDxC=0. ; yDxR=0. ; yCxL=0. ; yCxR=0.
  yUxL=0.  ; yUxC=0. ; yUxR=0. ; yCxC=1.

  orientation=bergs%initial_orientation*(pi/180)

  if (.not. bergs%hexagonal_icebergs) then ! Treat icebergs as squares during spreading to cells, rectangles during thermodynamics

    ! L is the non dimensional length of the iceberg [ L=(Area of berg/ Area of grid cell)^0.5 ].
    if (grd%area_um(i,j)>0) then
      L=min( sqrt(Area / grd%area_um(i,j)),1.0)
    else
      L=1.
    endif

    if (bergs%rotate_icebergs_for_mass_spreading) then

      !No bonds yet for call this:
      !call find_orientation_using_iceberg_bonds(grd,berg,orientation)

      !Subtracting the position of the nearest corner from x,y  (The mass will then be spread over the 4 cells connected to that corner)
      origin_x=1. ; origin_y=1.
      if (x<0.5) origin_x=0.
      if (y<0.5) origin_y=0.

      !Position of the square center, relative to origin at the nearest vertex
      x0=(x-origin_x)
      y0=(y-origin_y)

      call Square_into_quadrants_using_triangles(x0,y0,L,orientation,Area_square, Area_Q1, Area_Q2, Area_Q3, Area_Q4)

      if (min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4)) <-tol) then
        call error_mesg('KID, square spreading', 'Intersection with square should not be negative!!!', WARNING)
        write(stderrunit,*) 'KID, yU,yC,yD', Area_Q1, Area_Q2, Area_Q3, Area_Q4
      endif

      Area_Q1=Area_Q1/Area_square
      Area_Q2=Area_Q2/Area_square
      Area_Q3=Area_Q3/Area_square
      Area_Q4=Area_Q4/Area_square

      !Now, you decide which quadrant belongs to which mass on ocean cell.
      if ((x.ge. 0.5) .and. (y.ge. 0.5)) then !Top right vertex
        yUxR=Area_Q1
        yUxC=Area_Q2
        yCxC=Area_Q3
        yCxR=Area_Q4
      elseif ((x .lt. 0.5) .and. (y.ge. 0.5)) then  !Top left vertex
        yUxC=Area_Q1
        yUxL=Area_Q2
        yCxL=Area_Q3
        yCxC=Area_Q4
      elseif ((x.lt.0.5) .and. (y.lt. 0.5)) then !Bottom left vertex
        yCxC=Area_Q1
        yCxL=Area_Q2
        yDxL=Area_Q3
        yDxC=Area_Q4
      elseif ((x.ge.0.5) .and. (y.lt. 0.5)) then!Bottom right vertex
        yCxR=Area_Q1
        yCxC=Area_Q2
        yDxC=Area_Q3
        yDxR=Area_Q4
      endif

    else
      !no rotation for mass spreading. Given that bergs have not yet rotated anyway, this should be identical to the above,
      !as long as initial orientation is 0. Useful for debugging.

      xL=min(0.5, max(0., 0.5-(x/L)))
      xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
      xC=max(0., 1.-(xL+xR))
      yD=min(0.5, max(0., 0.5-(y/L)))
      yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
      yC=max(0., 1.-(yD+yU))

      yDxL=yD*xL
      yDxC=yD*xC
      yDxR=yD*xR
      yCxL=yC*xL
      yCxR=yC*xR
      yUxL=yU*xL
      yUxC=yU*xC
      yUxR=yU*xR
      yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

      !TODO: do you need to account for bergs that might overlap the edges of the model domain (if any)?
      !fraction_used=1. ! rectangular bergs do share mass with boundaries (all mass is included in cells)
    endif
  else ! hexagonal

    if (grd%area_um(i,j)>0) then
      ! Non-dimensionalize element length by grid area. (This gives the non-dim Apothem of the hexagon)
      H=min(( (sqrt(Area/(2.*sqrt(3.))) / sqrt(grd%area_um(i,j)))),1.)
    else
      ! Largest allowable H, since this makes S=0.49, and S has to be less than 0.5
      H=(sqrt(3.)/2)*(0.49)
    endif
    S=(2/sqrt(3.))*H !Side of the hexagon

    if (S>0.5) then
      ! The width of an iceberg should not be greater than half the grid cell, or else it can spread over 3 cells
      ! (i.e. S must be less than 0.5 non-dimensionally)
      !print 'Elements must be smaller than a whole grid cell', 'i.e.: S= ' , S , '>=0.5'
      call error_mesg('KID, hexagonal spreading', &
                      'Diameter of the iceberg is larger than a grid cell. Use smaller icebergs', WARNING)
    endif

    !Subtracting the position of the nearest corner from x,y
    !(The mass will then be spread over the 4 cells connected to that corner)
    origin_x=1. ; origin_y=1.
    if (x<0.5) origin_x=0.
    if (y<0.5) origin_y=0.

    !Position of the hexagon center, relative to origin at the nearest vertex
    x0=(x-origin_x)
    y0=(y-origin_y)

    !no bonds yet to call this:
    !if (bergs%rotate_icebergs_for_mass_spreading) call find_orientation_using_iceberg_bonds(grd,berg,orientation)

    call Hexagon_into_quadrants_using_triangles(x0,y0,H,orientation,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)

    if (min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4)) <-tol) then
      call error_mesg('KID, hexagonal spreading', 'Intersection with hexagons should not be negative!!!', WARNING)
      write(stderrunit,*) 'KID, yU,yC,yD', Area_Q1, Area_Q2, Area_Q3, Area_Q4
    endif

    Area_Q1=Area_Q1/Area_hex
    Area_Q2=Area_Q2/Area_hex
    Area_Q3=Area_Q3/Area_hex
    Area_Q4=Area_Q4/Area_hex

    !Now, you decide which quadrant belongs to which mass on ocean cell.
    if ((x.ge. 0.5) .and. (y.ge. 0.5)) then !Top right vertex
      yUxR=Area_Q1
      yUxC=Area_Q2
      yCxC=Area_Q3
      yCxR=Area_Q4
    elseif ((x .lt. 0.5) .and. (y.ge. 0.5)) then  !Top left vertex
      yUxC=Area_Q1
      yUxL=Area_Q2
      yCxL=Area_Q3
      yCxC=Area_Q4
    elseif ((x.lt.0.5) .and. (y.lt. 0.5)) then !Bottom left vertex
      yCxC=Area_Q1
      yCxL=Area_Q2
      yDxL=Area_Q3
      yDxC=Area_Q4
    elseif ((x.ge.0.5) .and. (y.lt. 0.5)) then!Bottom right vertex
      yCxR=Area_Q1
      yCxC=Area_Q2
      yDxC=Area_Q3
      yDxR=Area_Q4
    endif
  endif
end subroutine calving_tabular_particle_grid_overlap

!> Returns a list of partially-full bergs connected to full bergs
recursive subroutine make_list_of_bonded_to_full(berg,first)
  type(iceberg), pointer :: berg !< Berg to process
  type(iceberg), pointer :: first !< The first berg in the list of tabular bergs
  ! Local variables
  type(bond), pointer :: current_bond
  type(iceberg), pointer :: other_berg

  berg%sss=0
  current_bond=>berg%first_bond
  do while (associated(current_bond))
    if  (associated(current_bond%other_berg)) then
      other_berg=>current_bond%other_berg
      if (other_berg%id>0) then
        !this berg has not been processed
        other_berg%id=-other_berg%id
        if (other_berg%static_berg==2) then
          call make_list_of_bonded_to_full(other_berg, first)
        else
          call insert_tabular_particle_into_list(first, other_berg)
          other_berg%sss=1 !Marks the berg as 1 bond away from a filled berg
        endif
      endif
    endif
    current_bond=>current_bond%next_bond
  enddo
end subroutine make_list_of_bonded_to_full

!> This is the same as make_list_of_bonded_to_full except the list is of bergs bonded to
!! bergs with static_berg=2.5 (partially-full, and not eventually connecting to a full berg,
!! but overlapping a full ice shelf grid cell).
recursive subroutine make_list_of_bonded_to_full2(berg,first)
  type(iceberg), pointer :: berg !< Berg to process
  type(iceberg), pointer :: first !< The first berg in the list of tabular bergs
  ! Local variables
  type(bond), pointer :: current_bond
  type(iceberg), pointer :: other_berg
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  berg%sss=1
  current_bond=>berg%first_bond
  do while (associated(current_bond))
    if  (associated(current_bond%other_berg)) then
      other_berg=>current_bond%other_berg
      if (other_berg%id>0) then
        !this berg has not been processed
        other_berg%id=-other_berg%id
        if (other_berg%static_berg==2.5) then
          call make_list_of_bonded_to_full2(other_berg, first)
        elseif (other_berg%static_berg==3) then
          call insert_tabular_particle_into_list(first, other_berg)
          !The parent berg with static_berg=2.5 is treated as if it is 1 bond away from a filled berg
          !(even though it is not). Bergs with static_berg=3 that are bonded to the parent berg are
          !marked so that they are treated as if they are 2 bonds away from a filled berg.
          other_berg%sss=2
        else
          write(stderrunit,*) 'KID, make_list_of_bonded_to_full2: something went very wrong!', other_berg%static_berg
          call error_mesg('KID, make_list_of_bonded_to_full2',&
            'Error in determining bonds away from the active ice front!!', FATAL)
        endif
      endif
    endif
    current_bond=>current_bond%next_bond
  enddo
end subroutine make_list_of_bonded_to_full2

!> Determine how many bonds away from a full particle each partially-full particle is.
subroutine assign_bonds_from_full(first,count)
  type(iceberg), pointer :: first !< The first berg to add to the tabular list
  integer :: count !> tracks number of bonds away from a "full" particle
  ! Local variables
  type(iceberg), pointer :: berg ! Berg to process
  type(bond), pointer :: current_bond
  type(iceberg), pointer :: other_berg, prev_berg

  berg=>first
  do while (associated(berg))
    current_bond=>berg%first_bond
    do while (associated(current_bond))
      if  (associated(current_bond%other_berg)) then
        other_berg=>current_bond%other_berg
        if (other_berg%id>0) then
          !this berg has not been processed
          other_berg%id=-other_berg%id
          !Add the berg to the front of the list of
          !partially-full tabular berg particles
          call insert_tabular_particle_into_list(first, other_berg)
          other_berg%sss=count
        endif
      endif
      current_bond=>current_bond%next_bond
    enddo
    prev_berg=>berg

    if (associated(berg%next_t)) then
      berg=>berg%next_t
      call delete_tabular_particle_from_list(first, prev_berg)
    else
      !Delete the berg that was just processed from the list
      !If it is at the start of the list, the the whole list is nullified.
      call delete_tabular_particle_from_list(first, prev_berg)
      !the end of the list has been reached. Restart from the beginning of the list,
      !which may have new bergs added
      berg=>first
      count=count+1
    endif
  enddo
end subroutine assign_bonds_from_full

!> Inserts a berg into the front of a list of tabular iceberg particles
subroutine insert_tabular_particle_into_list(first, newberg)
  ! Arguments
  type(iceberg), pointer :: first !< The first berg in the list of tabular bergs
  type(iceberg), pointer :: newberg !< New berg to insert

  if (associated(first)) then
    !must be inserted at front of the list
    newberg%next_t=>first
    newberg%prev_t=>null()
    first%prev_t=>newberg
    first=>newberg
  else
    ! list is empty so create it
    first=>newberg
    first%next_t=>null()
    first%prev_t=>null()
  endif

end subroutine insert_tabular_particle_into_list

!> Remove a berg from the list of tabular iceberg particles
subroutine delete_tabular_particle_from_list(first, this)
  ! Arguments
  type(iceberg), pointer :: first !< List of tabular bergs
  type(iceberg), pointer :: this !< Berg to be deleted


  ! Connect neighbors to each other
  if (associated(this%prev_t)) this%prev_t%next_t=>this%next_t
  if (associated(this%next_t)) this%next_t%prev_t=>this%prev_t

  ! If deleting the first berg particle, need to reassign the next particle
  ! as the front of the list
  if (this%id.eq.first%id) first=>this%next_t

  this%prev_t=>NULL()
  this%next_t=>NULL()

end subroutine delete_tabular_particle_from_list

end module ice_shelf_tabular_calving
