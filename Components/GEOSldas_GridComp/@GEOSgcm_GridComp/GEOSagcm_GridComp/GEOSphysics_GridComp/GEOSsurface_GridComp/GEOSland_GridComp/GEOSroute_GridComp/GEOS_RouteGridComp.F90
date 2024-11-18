
#include "MAPL_Generic.h"

#ifndef RUN_FOR_REAL
#define MAPL_DimsCatchOnly MAPL_DimsTileOnly
#endif

!=============================================================================
module GEOS_RouteGridCompMod

!BOP
! !MODULE: GEOS_Route -- child to the "Land" gridded component.  

! !DESCRIPTION:
!   {\tt GEOS\_Route} is a gridded component to route total runoff produced in
!   {\tt GEOS\_Catch} (RUNOFF in {\tt GEOScatch\_GridComp} or {\tt GEOScatchCN\_GridComp}) through a 290,188   
!   watersheds globally (excluding submerged catchments (watersheds) from the global list of 291,284. 
!   All of its calculations are done on Pfafstetter watershed space. {\tt GEOS\_Route} has no children. \\
!
!   IMPORTS   : RUNOFF \\
!   INTERNALS : AREACAT, LENGSC2, DNSTR, WSTREAM, WRIVER, LRIVERMOUTH, ORIVERMOUTH \\
!   EXPORTS   : QSFLOW, QINFLOW, QOUTFLOW \\

! !USES: 

  use ESMF
  use MAPL_Mod
  use MAPL_ConstantsMod
  use ROUTING_MODEL,          ONLY:     &
       river_routing, ROUTE_DT
#if 0
  USE catch_constants, ONLY:          &
       N_CatG => N_Pfaf_Catchs  
#endif
  
  implicit none
  integer, parameter :: N_CatG = 291284
  integer,parameter :: nmax=150  
  integer,parameter :: upmax=34
  !integer, parameter :: nt_all = 112573
  private

  type T_RROUTE_STATE
     private
     type (ESMF_RouteHandle) :: routeHandle
     type (ESMF_Field)       :: field
     integer :: nTiles
     integer :: nt_global
     integer :: nt_local
     integer :: comm
     integer :: nDes
     integer :: myPe
     integer :: minCatch
     integer :: maxCatch
     integer, pointer :: pfaf(:) => NULL()
     real,    pointer :: tile_area(:) => NULL() !m2
     integer, pointer :: nsub(:) => NULL()
     integer, pointer :: subi(:,:) => NULL()
     real,    pointer :: subarea(:,:) => NULL() !m2

     integer, pointer :: scounts_global(:) => NULL()
     integer, pointer :: rdispls_global(:) => NULL()
     integer, pointer :: scounts_cat(:) => NULL()
     integer, pointer :: rdispls_cat(:) => NULL()
 
     real,    pointer :: runoff_save(:) => NULL()
     real,    pointer :: areacat(:) => NULL() !m2
     real,    pointer :: lengsc(:) => NULL() !m

     real,    pointer :: wstream(:) => NULL() !m3
     real,    pointer :: wriver(:)  => NULL() !m3
     integer, pointer :: downid(:) => NULL()
     integer, pointer :: upid(:,:) => NULL()

  end type T_RROUTE_STATE

  ! Wrapper for extracting internal state
  ! -------------------------------------
  type RROUTE_WRAP
     type (T_RROUTE_STATE), pointer :: PTR => null()
  end type RROUTE_WRAP

  include "mpif.h"

  
! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (ESMF_Config          )            :: CF
    
    type (T_RROUTE_STATE), pointer         :: route_internal_state => null()
    type (RROUTE_wrap)                     :: wrap

    integer      :: RUN_DT
    real         :: DT

!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC                                 ,&
                          NAME=COMP_NAME			   ,&
                          RC=STATUS )

    VERIFY_(STATUS)

    Iam = trim(COMP_NAME) // 'SetServices'

! -----------------------------------------------------------
! Set the Initialize and Run entry points
! -----------------------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
!    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, RC=STATUS)
!    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, RUN1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, RUN2, RC=STATUS )
    VERIFY_(STATUS)

!------------------------------------------------------------
! Set generic final method 
!------------------------------------------------------------

    
! -----------------------------------------------------------
! Get the configuration
! -----------------------------------------------------------
! 
    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)
!
! -----------------------------------------------------------
! Get the intervals
! -----------------------------------------------------------
!
    call ESMF_ConfigGetAttribute (CF, DT                         ,&
                                  Label="RUN_DT:"                ,&
                                  RC=STATUS)

    VERIFY_(STATUS)

    RUN_DT = nint(DT)

! -----------------------------------------------------------
! At the moment, this will refresh when the land parent 
! needs to refresh.

    call ESMF_ConfigGetAttribute ( CF, DT, Label=trim(COMP_NAME)//"_DT:", &
         default=DT, RC=STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Set the state variable specs.
! -----------------------------------------------------------

!BOS

! -----------------------------------------------------------
!   Import States
! -----------------------------------------------------------

    call MAPL_AddImportSpec(GC,                          &
         LONG_NAME          = 'runoff_flux'               ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RUNOFF'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

! -----------------------------------------------------------
!   INTERNAL STATE
! -----------------------------------------------------------

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'area_of_catchment'        ,&
         UNITS              = 'km+2'                     ,&
         SHORT_NAME         = 'AREACAT'                  ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 
 
    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'length_of_channel_segment',&
         UNITS              = 'km+2'                     ,&
         SHORT_NAME         = 'LENGSC'                   ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'index_of_downtream_catchment',&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'DNSTR'                    ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'volume_of_water_in_local_stream',&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WSTREAM'                  ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  )
    
    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'volume_of_water_in_river' ,&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WRIVER'                   ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'TileID_of_the_lake_tile_at_the_river_mouth' ,&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'LRIVERMOUTH'              ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'TileID_of_the_ocean_tile_at_the_river_mouth' ,&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'ORIVERMOUTH'              ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   
! -----------------------------------------------------------
!  EXPORT STATE:
! -----------------------------------------------------------

    call MAPL_AddExportSpec(GC,                        &
         LONG_NAME          = 'transfer_of_moisture_from_stream_variable_to_river_variable' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QSFLOW'                   ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
                                         RC=STATUS  ) 

    call MAPL_AddExportSpec(GC,                    &
         LONG_NAME          = 'transfer_of_river_water_from_upstream_catchments' ,&
         UNITS              = 'm+3 s-1'                   ,&
         SHORT_NAME         = 'QINFLOW'                   ,&
         DIMS               = MAPL_DimsCatchOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
                                          RC=STATUS  ) 

    call MAPL_AddExportSpec(GC,                    &
         LONG_NAME          = 'transfer_of_river_water_to_downstream_catchments' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QOUTFLOW'                 ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
                                         RC=STATUS  ) 

!EOS

    !call MAPL_TimerAdd(GC,    name="RUN1"  ,RC=STATUS)
    !VERIFY_(STATUS)
    !call MAPL_TimerAdd(GC,    name="RUN2"  ,RC=STATUS)
    !VERIFY_(STATUS)    
    call MAPL_TimerAdd(GC,    name="-RRM" ,RC=STATUS)
    VERIFY_(STATUS)



! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( route_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%ptr => route_internal_state
    
! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)

! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN1"           ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN2"           ,RC=STATUS)
    VERIFY_(STATUS)

! All done
!---------
    
    call MAPL_GenericSetServices(GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

! -----------------------------------------------------------
! INITIALIZE -- Initialize method for the route component
! -----------------------------------------------------------

  subroutine INITIALIZE (GC,IMPORT, EXPORT, CLOCK, RC )

! -----------------------------------------------------------
! !ARGUMENTS:
! -----------------------------------------------------------

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

! -----------------------------------------------------------
! ErrLog Variables
! -----------------------------------------------------------

    character(len=ESMF_MAXSTR)          :: IAm="Initialize"
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! -----------------------------------------------------------
! Locals
! -----------------------------------------------------------

    type (ESMF_VM) :: VM
    integer        :: comm
    integer        :: nDEs
    integer        :: myPE
    integer        :: beforeMe, minCatch, maxCatch, pf, i
    integer        :: ntiles, nt_global

    type(ESMF_Grid) :: tileGrid
    type(ESMF_Grid) :: newTileGrid
    type(ESMF_Grid) :: catchGrid
    type(ESMF_DistGrid) :: distGrid
    type(ESMF_Field) :: field, field0

    type(MAPL_MetaComp), pointer   :: MAPL
    type(MAPL_LocStream) :: locstream
    
    integer, pointer :: ims(:) => NULL()
    integer, pointer :: pfaf(:) => NULL()
    integer, pointer :: arbSeq(:) => NULL()
    integer, pointer :: arbSeq_pf(:) => NULL()    
    integer, pointer :: arbSeq_ori(:) => NULL()    
    integer, allocatable :: arbIndex(:,:)
    real, pointer :: tile_area_src(:) => NULL()
    integer,pointer :: local_id(:)  => NULL()
    real, pointer :: tile_area_local(:) => NULL(), tile_area_global(:) => NULL()
    real, pointer :: tile_area(:) => NULL()    
    real, pointer :: ptr2(:) => NULL()

    real,pointer :: subarea_global(:,:)=> NULL(),subarea(:,:)=> NULL() ! Arrays for sub-area and fractions
    integer,pointer :: subi_global(:,:)=> NULL(),subi(:,:)=> NULL()
    integer,pointer :: nsub_global(:)=> NULL(),nsub(:)=> NULL()
    real,pointer :: area_cat_global(:)=> NULL(),area_cat(:)=> NULL()
    integer,pointer :: scounts(:)=>NULL()
    integer,pointer :: scounts_global(:)=>NULL(),rdispls_global(:)=>NULL()
    integer,pointer :: scounts_cat(:)=>NULL(),rdispls_cat(:)=>NULL()    

    real,pointer :: runoff_save(:)=>NULL(), areacat(:)=>NULL()
    real,pointer :: lengsc_global(:)=>NULL(), lengsc(:)=>NULL()
    integer,pointer :: downid_global(:)=>NULL(), downid(:)=>NULL()
    integer,pointer :: upid_global(:,:)=>NULL(), upid(:,:)=>NULL()    

    real,pointer :: wstream(:)=>NULL(),wriver(:)=>NULL()
    
    type (T_RROUTE_STATE), pointer         :: route => null()
    type (RROUTE_wrap)                     :: wrap


    real, pointer :: dataPtr(:)
    integer :: j,nt_local,mpierr,it
    ! ------------------
    ! begin

    
    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)

    route => wrap%ptr

    ! get vm
    ! extract comm
    call ESMF_VMGetCurrent(VM,                                RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet       (VM,       mpiCommunicator =comm,   RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet       (VM, localpet=MYPE, petcount=nDEs,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    route%comm = comm
    route%ndes = ndes
    route%mype = mype
 
    allocate(ims(1:ndes))
    ! define minCatch, maxCatch
    call MAPL_DecomposeDim ( n_catg,ims,ndes ) ! ims(mype+1) gives the size of my partition
    ! myPE is 0-based!
    beforeMe = sum(ims(1:mype))
    minCatch = beforeMe + 1
    maxCatch = beforeMe + ims(myPe+1)
    !print *, "my PE is:",mype,", minCatch is:",minCatch,", maxCatch is:",maxCatch
 
    ! get LocStream
    call MAPL_Get(MAPL, LocStream = locstream, RC=status)
    VERIFY_(STATUS) 
    ! extract Pfaf (TILEI on the "other" grid)    
    call MAPL_LocStreamGet(locstream, &
         tileGrid=tilegrid, nt_global=nt_global, RC=status)
    route%nt_global = nt_global
    !if (mapl_am_I_root()) print *, "nt_global=",nt_global           
    allocate(pfaf(nt_global))
    open(77,file="../input/pfaf_input.txt",status="old",action="read")
    read(77,*)pfaf
    close(77)
    ! exchange Pfaf across PEs

    ntiles = 0
    !loop over total_n_tiles
    allocate(arbSeq_ori(1:nt_global))
    do i = 1, nt_global
       pf = pfaf(i)
       if (pf >= minCatch .and. pf <= maxCatch) then ! I want this!
          !print *,"my PE is:",mype,"pf=",pf
          ntiles = ntiles+1
          !realloc if needed
          arbSeq_ori(ntiles) = i
       end if
    end do ! global tile loop
    !if (mapl_am_I_root()) print *, "ntiles:",ntiles
    allocate(arbSeq(ntiles))
    arbSeq=arbSeq_ori(1:ntiles)
    deallocate(arbSeq_ori)
  
    distgrid = ESMF_DistGridCreate(arbSeqIndexList=arbSeq, rc=status)
    VERIFY_(STATUS)

    newTileGRID = ESMF_GridEmptyCreate(rc=status)
    VERIFY_(STATUS)       
    allocate(arbIndex(nTiles,1), stat=status)
    VERIFY_(STATUS)

    arbIndex(:,1) = arbSeq

    call ESMF_GridSet(newTileGrid,  &
         name='redist_tile_grid_for_'//trim(COMP_NAME),    &
         distgrid=distgrid, & 
         gridMemLBound=(/1/), &
         indexFlag=ESMF_INDEX_USER, &
         distDim = (/1/), &
         localArbIndexCount=ntiles, &
         localArbIndex=arbIndex, &
         minIndex=(/1/), &
         maxIndex=(/NT_GLOBAL/), &
         rc=status)
    VERIFY_(STATUS)

    deallocate(arbIndex)

    call ESMF_GridCommit(newTileGrid, rc=status)
    VERIFY_(STATUS)   

    ! now create a "catch" grid to be the "native" grid for this component
    distgrid = ESMF_DistGridCreate(arbSeqIndexList=(/minCatch:maxCatch/), &
         rc=status)
    VERIFY_(STATUS)  
    catchGRID = ESMF_GridEmptyCreate(rc=status)
    VERIFY_(STATUS)

    allocate(arbIndex(ims(myPE+1),1), stat=status)
    VERIFY_(STATUS)

    arbIndex(:,1) = (/minCatch:maxCatch/)
         
    call ESMF_GridSet(catchGrid,  &
         name='catch_grid_for_'//trim(COMP_NAME),    &
         distgrid=distgrid, & 
         gridMemLBound=(/1/), &
         indexFlag=ESMF_INDEX_USER, &
         distDim = (/1/), &
         localArbIndexCount=ims(myPE+1), &
         localArbIndex=arbIndex, &
         minIndex=(/1/), &
         maxIndex=(/N_CatG/), &
         rc=status)
    VERIFY_(STATUS)

    deallocate(arbIndex)
    call ESMF_GridCommit(catchGrid, rc=status)
    VERIFY_(STATUS)
    call ESMF_GridCompSet(gc, grid=catchGrid, RC=status)
    VERIFY_(STATUS)   
    call MAPL_LocStreamGet(locstream, TILEAREA = tile_area_src, LOCAL_ID=local_id, RC=status)

    VERIFY_(STATUS) 
    field0 = ESMF_FieldCreate(grid=tilegrid, datacopyflag=ESMF_DATACOPY_VALUE, &
         farrayPtr=tile_area_src, name='TILE_AREA_SRC', RC=STATUS)
    VERIFY_(STATUS)
    ! create field on the "new" tile grid
    allocate(tile_area(ntiles), stat=status)
    VERIFY_(STATUS)
    field = ESMF_FieldCreate(grid=newtilegrid, datacopyflag=ESMF_DATACOPY_VALUE, &
         farrayPtr=tile_area, name='TILE_AREA', RC=STATUS)  
  
    VERIFY_(STATUS)
    ! create routehandle
    call ESMF_FieldRedistStore(srcField=field0, dstField=field, &
                routehandle=route%routehandle, rc=status)
    VERIFY_(STATUS)  
    ! redist tile_area
    call ESMF_FieldRedist(srcField=FIELD0, dstField=FIELD, &
         routehandle=route%routehandle, rc=status)
    VERIFY_(STATUS)


    ntiles = 0
    !loop over total_n_tiles
    allocate(arbSeq_ori(1:nt_global))
    do i = 1, nt_global
       pf = pfaf(i)
       if (pf >= minCatch .and. pf <= maxCatch) then ! I want this!
          !print *,"my PE is:",mype,"pf=",pf
          ntiles = ntiles+1
          !realloc if needed
          arbSeq_ori(ntiles) = pf
       end if
    end do ! global tile loop
    ntiles = maxCatch-minCatch+1
    !if (mapl_am_I_root()) print *, "ntiles:",ntiles
!    allocate(arbSeq_pf(ntiles))
    allocate(arbSeq_pf(maxCatch-minCatch+1))
    arbSeq_pf = [(i, i = minCatch, maxCatch)]
    deallocate(arbSeq_ori)
   
    ! redist pfaf (NOTE: me might need a second routehandle for integers)

    route%pfaf => arbSeq_pf
    route%ntiles = ntiles  
    route%minCatch = minCatch
    route%maxCatch = maxCatch 
    allocate(ptr2(ntiles), stat=status)
    VERIFY_(STATUS)
    route%field = ESMF_FieldCreate(grid=catchGrid, datacopyflag=ESMF_DATACOPY_VALUE, &
        farrayPtr=ptr2, name='RUNOFF', RC=STATUS)
    VERIFY_(STATUS)
  ! Read sub-area data from text files
    allocate(nsub_global(N_CatG),subarea_global(nmax,N_CatG))
    open(77,file="../input/Pfaf_nsub_M36.txt",status="old",action="read"); read(77,*)nsub_global; close(77)
    open(77,file="../input/Pfaf_asub_M36.txt",status="old",action="read"); read(77,*)subarea_global; close(77)       
    allocate(nsub(ntiles),subarea(nmax,ntiles))
    nsub=nsub_global(minCatch:maxCatch)
    subarea=subarea_global(:,minCatch:maxCatch)
    subarea=subarea*1.e6 !km2->m2
    deallocate(nsub_global,subarea_global)

    route%nsub => nsub
    route%subarea => subarea
    
    allocate(subi_global(nmax,N_CatG),subi(nmax,ntiles))
    open(77,file="../input/Pfaf_isub_M36.txt",status="old",action="read");read(77,*)subi_global;close(77)
    subi=subi_global(:,minCatch:maxCatch)
    route%subi => subi
    deallocate(subi_global)

    nt_local=size(tile_area_src,1)

    allocate(scounts(ndes),scounts_global(ndes),rdispls_global(ndes))
    scounts=0
    scounts(mype+1)=nt_local  
    call MPI_Allgather(scounts(mype+1), 1, MPI_INTEGER, scounts_global, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr) 
    rdispls_global(1)=0
    do i=2,nDes
      rdispls_global(i)=rdispls_global(i-1)+scounts_global(i-1)
    enddo
    deallocate(scounts)
    route%scounts_global=>scounts_global
    route%rdispls_global=>rdispls_global

    allocate(scounts(ndes),scounts_cat(ndes),rdispls_cat(ndes))
    scounts=0
    scounts(mype+1)=ntiles  
    call MPI_Allgather(scounts(mype+1), 1, MPI_INTEGER, scounts_cat, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr) 
    rdispls_cat(1)=0
    do i=2,nDes
      rdispls_cat(i)=rdispls_cat(i-1)+scounts_cat(i-1)
    enddo
    deallocate(scounts)
    route%scounts_cat=>scounts_cat
    route%rdispls_cat=>rdispls_cat

    route%nt_local=nt_local
    allocate(runoff_save(1:nt_local))
    route%runoff_save => runoff_save
    route%runoff_save=0.

    allocate(tile_area_local(nt_local),tile_area_global(nt_global))  
    open(77,file="../input/area_m36_1d.txt",status="old",action="read");read(77,*)tile_area_global;close(77)
    tile_area_local=tile_area_global(rdispls_global(mype+1)+1:rdispls_global(mype+1)+nt_local)*1.e6 !km2->m2
    route%tile_area => tile_area_local
    deallocate(tile_area_global)

    allocate(areacat(1:ntiles))
    areacat=0. 
    do i=1,ntiles
      do j=1,nmax
        it=route%subi(j,i) 
        if(it>0)then 
          areacat(i)=areacat(i)+route%subarea(j,i)
        endif
        if(it==0)exit
      enddo
    enddo  
    route%areacat=>areacat

    allocate(lengsc_global(n_catg),lengsc(ntiles))   
    open(77,file="../input/Pfaf_lriv_PR.txt",status="old",action="read");read(77,*)lengsc_global;close(77)
    lengsc=lengsc_global(minCatch:maxCatch)*1.e3 !km->m
    route%lengsc=>lengsc
    deallocate(lengsc_global)

    allocate(downid_global(n_catg),downid(ntiles))
    open(77,file="../input/downstream_1D_new_noadj.txt",status="old",action="read");read(77,*)downid_global;close(77)    
    downid=downid_global(minCatch:maxCatch)
    route%downid=>downid
    deallocate(downid_global)

    allocate(upid_global(upmax,n_catg),upid(upmax,ntiles))   
    open(77,file="../input/upstream_1D.txt",status="old",action="read");read(77,*)upid_global;close(77)  
    upid=upid_global(:,minCatch:maxCatch)   
    route%upid=>upid
    deallocate(upid_global)

    allocate(wriver(ntiles),wstream(ntiles))
!    allocate(wriver_global(n_catg),wstream_global(n_catg))
!    open(77,file="../input/restart/wriver_global_.txt",status="old",action="read");read(77,*)wriver_global;close(77)    
!    open(77,file="../input/restart/wstream_global_.txt",status="old",action="read");read(77,*)wstream_global;close(77) 
!    wriver=wriver_global(minCatch:maxCatch)
!    wstream=wstream_global(minCatch:maxCatch)
!    deallocate(wriver_global,wstream_global)
    wriver=0.
    wstream=0.
    route%wstream=>wstream
    route%wriver=>wriver
    !This should be read from restart file
!    route%wstream=0.
!    route%wriver=0.

    !if (mapl_am_I_root())then
    !  open(88,file="nsub.txt",action="write")
    !  open(89,file="subarea.txt",action="write")
    !  open(90,file="subi.txt",action="write")
    !  open(91,file="tile_area.txt",action="write")
    !  do i=1,nTiles
    !    write(88,*)route%nsub(i)
    !    write(89,'(150(1x,f10.4))')route%subarea(:,i)
    !    write(90,'(150(i7))')route%subi(:,i)
    !    write(91,*)route%tile_area(i)
    !  enddo
    !  stop
    !endif
   
    deallocate(ims)
    call MAPL_GenericInitialize ( GC, import, export, clock, rc=status )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)      
  end subroutine INITIALIZE
  
! -----------------------------------------------------------
! RUN -- Run method for the route component
! -----------------------------------------------------------
  subroutine RUN1 (GC,IMPORT, EXPORT, CLOCK, RC )

! -----------------------------------------------------------
! !ARGUMENTS:
! -----------------------------------------------------------

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC
  end subroutine RUN1


  subroutine RUN2 (GC,IMPORT, EXPORT, CLOCK, RC )

! -----------------------------------------------------------
! !ARGUMENTS:
! -----------------------------------------------------------

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

! -----------------------------------------------------------
! ErrLog Variables
! -----------------------------------------------------------

    character(len=ESMF_MAXSTR)          :: IAm="Run2"
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! -----------------------------------------------------------
! Locals
! -----------------------------------------------------------

    type (MAPL_MetaComp),     pointer   :: MAPL
    type (ESMF_State       )            :: INTERNAL
!    type(ESMF_Alarm)                    :: ALARM
    type (ESMF_Config )                 :: CF
    type(ESMF_VM)                       :: VM

! -----------------------------------------------------
! IMPORT pointers
! ----------------------------------------------------- 

    real, dimension(:), pointer :: RUNOFF 
    real, dimension(:), pointer :: RUNOFF_SRC0   

! -----------------------------------------------------
! INTERNAL pointers
! ----------------------------------------------------- 

    real, dimension(:), pointer :: AREACAT
    real, dimension(:), pointer :: LENGSC
    real, dimension(:), pointer :: DNSTR
    real, dimension(:), pointer :: WSTREAM
    real, dimension(:), pointer :: WRIVER
    real, dimension(:), pointer :: LRIVERMOUTH
    real, dimension(:), pointer :: ORIVERMOUTH

! -----------------------------------------------------
! EXPORT pointers 
! -----------------------------------------------------

    real, dimension(:), pointer :: QSFLOW
    real, dimension(:), pointer :: QINFLOW
    real, dimension(:), pointer :: QOUTFLOW
  
! Time attributes and placeholders

!    type(ESMF_Time) :: CURRENT_TIME

! Others

    type(ESMF_Grid)                    :: TILEGRID
    type (MAPL_LocStream)              :: LOCSTREAM
 
    integer                            :: NTILES, N_CatL, N_CYC
    logical, save                      :: FirstTime=.true.
    real, pointer, dimension(:)    :: tile_area
    integer, pointer, dimension(:) :: pfaf_code

    INTEGER, DIMENSION(:,:), POINTER, SAVE   :: AllActive,DstCatchID 
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: srcProcsID, LocDstCatchID  
    integer, dimension (:),allocatable, SAVE :: GlbActive
    INTEGER, SAVE                            :: N_Active, ThisCycle   
    INTEGER                                  :: Local_Min, Local_Max
    integer                                  :: K, N, I, req
    REAL                                     :: mm2m3, rbuff, HEARTBEAT 
    REAL, ALLOCATABLE, DIMENSION(:)          :: RUNOFF_CATCH, RUNOFF_ACT,AREACAT_ACT,& 
         LENGSC_ACT, QSFLOW_ACT,QOUTFLOW_ACT
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_index
    type(ESMF_Field) :: runoff_src

    integer                                :: ndes, mype
    type (T_RROUTE_STATE), pointer         :: route => null()
    type (RROUTE_wrap)                     :: wrap
    INTEGER, DIMENSION(:)  ,ALLOCATABLE  :: scounts, scounts_global,rdispls, rcounts  
    real, dimension(:), pointer :: runoff_global,runoff_local,area_local,runoff_cat_global    

    integer :: mpierr, nt_global,nt_local, it, j, upid,cid,temp(1),tid,istat

    type(ESMF_Time) :: CurrentTime   
    integer :: YY,MM,DD,HH,MMM,SS
    character(len=4) :: yr_s
    character(len=2) :: mon_s,day_s

    real,pointer :: runoff_save(:)=>NULL()
    real,pointer :: WSTREAM_ACT(:)=>NULL()
    real,pointer :: WRIVER_ACT(:)=>NULL()
    real,allocatable :: runoff_save_m3(:),runoff_global_m3(:),QOUTFLOW_GLOBAL(:)
    real,allocatable :: WTOT_BEFORE(:),WTOT_AFTER(:),QINFLOW_LOCAL(:),UNBALANCE(:),UNBALANCE_GLOBAL(:),ERROR(:),ERROR_GLOBAL(:)
    real,allocatable :: QFLOW_SINK(:),QFLOW_SINK_GLOBAL(:),WTOT_BEFORE_GLOBAL(:),WTOT_AFTER_GLOBAL(:)
    real,allocatable :: wriver_global(:),wstream_global(:),qsflow_global(:)
    

    ! ------------------
    ! begin    
    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)
    route => wrap%ptr

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS) 
    Iam = trim(COMP_NAME) // "RUN2"

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)
    call MAPL_Get(MAPL, HEARTBEAT = HEARTBEAT, RC=STATUS)
    VERIFY_(STATUS)
    !if (mapl_am_I_root()) print *, "HEARTBEAT=",HEARTBEAT 
! Start timers
! ------------

    call MAPL_TimerOn(MAPL,"RUN2")
! Get parameters from generic state
! ---------------------------------

    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
    VERIFY_(STATUS) 
! get pointers to inputs variables
! ----------------------------------

    ndes = route%ndes
    mype = route%mype  
    ntiles = route%ntiles  
    nt_global = route%nt_global  
    runoff_save => route%runoff_save
    nt_local = route%nt_local

    ! get the field from IMPORT
    call ESMF_StateGet(IMPORT, 'RUNOFF', field=runoff_src, RC=STATUS)
    VERIFY_(STATUS)    
    call ESMF_FieldGet(runoff_src, farrayPtr=RUNOFF_SRC0, rc=status)   
    VERIFY_(STATUS) 


! get pointers to internal variables
! ---------------------------------- 
    call MAPL_GetPointer(INTERNAL, AREACAT , 'AREACAT', RC=STATUS)
    VERIFY_(STATUS)    
    call MAPL_GetPointer(INTERNAL, LENGSC  , 'LENGSC',  RC=STATUS)
    VERIFY_(STATUS)        
    call MAPL_GetPointer(INTERNAL, DNSTR   , 'DNSTR'  , RC=STATUS)
    VERIFY_(STATUS)         
    call MAPL_GetPointer(INTERNAL, WSTREAM , 'WSTREAM', RC=STATUS)
    VERIFY_(STATUS)     
    call MAPL_GetPointer(INTERNAL, WRIVER  , 'WRIVER' , RC=STATUS)
    VERIFY_(STATUS)         
    call MAPL_GetPointer(INTERNAL, LRIVERMOUTH, 'LRIVERMOUTH' , RC=STATUS)
    VERIFY_(STATUS)        
    call MAPL_GetPointer(INTERNAL, ORIVERMOUTH, 'ORIVERMOUTH' , RC=STATUS)
    VERIFY_(STATUS)  
! get pointers to EXPORTS
! -----------------------

    call MAPL_GetPointer(EXPORT, QSFLOW,   'QSFLOW'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QINFLOW,  'QINFLOW' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QOUTFLOW, 'QOUTFLOW', RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_Get(MAPL, LocStream=LOCSTREAM, RC=STATUS)
    VERIFY_(STATUS)   
    call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS)
    VERIFY_(STATUS)    
    call MAPL_TimerOn  ( MAPL, "-RRM" )

    ! For efficiency, the time step to call the river routing model is set at ROUTE_DT 

    N_CYC = ROUTE_DT/HEARTBEAT    
    RUN_MODEL : if (ThisCycle == N_CYC) then   

       runoff_save = runoff_save + RUNOFF_SRC0/real (N_CYC)

       call ESMF_ClockGet(clock, currTime=CurrentTime, rc=status)
       VERIFY_(status)
       call ESMF_TimeGet(CurrentTime, yy=YY, mm=MM, dd=DD, h=HH, m=MMM, s=SS, rc=rc)  

       allocate(runoff_global(nt_global))
       call MPI_allgatherv  (                          &
          runoff_save,  route%scounts_global(mype+1)      ,MPI_REAL, &
          runoff_global, route%scounts_global, route%rdispls_global,MPI_REAL, &
          MPI_COMM_WORLD, mpierr) 

       allocate(RUNOFF_ACT(ntiles))
       RUNOFF_ACT=0.
       do i=1,ntiles
         do j=1,nmax
           it=route%subi(j,i) 
           if(it>0)then
             RUNOFF_ACT(i)=RUNOFF_ACT(i)+route%subarea(j,i)*runoff_global(it)/1000.   
           endif
           if(it==0)exit
         enddo
       enddo  

       deallocate(runoff_global) 

      !----check runoff balance----------------------------------------
       IF(1==0)THEN
         allocate(runoff_save_m3(nt_local),runoff_global_m3(nt_global))
         runoff_save_m3=runoff_save*route%tile_area/1000. 
         call MPI_allgatherv  (                          &
              runoff_save_m3,  route%scounts_global(mype+1)      ,MPI_REAL, &
              runoff_global_m3, route%scounts_global, route%rdispls_global,MPI_REAL, &
              MPI_COMM_WORLD, mpierr) 
         allocate(runoff_cat_global(n_catg) )  
         call MPI_allgatherv  (                          &
              RUNOFF_ACT,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              runoff_cat_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)     
         if(mapl_am_I_root())then 
           open(88,file="runoff_global_m3.txt",status="unknown", position="append")
           write(88,*)sum(runoff_global_m3)
           close(88)
           open(88,file="runoff_cat_global.txt",status="unknown", position="append")
           write(88,*)sum(runoff_cat_global)
           close(88)  
           print *,"sum(runoff_global_m3)=",sum(runoff_global_m3)
           print *,"sum(runoff_cat_global)",sum(runoff_cat_global)   
         endif   
         deallocate(runoff_save_m3,runoff_global_m3,runoff_cat_global)
       ENDIF 
       !--------------------------------------------

       allocate (AREACAT_ACT (1:ntiles))       
       allocate (LENGSC_ACT  (1:ntiles))
       allocate (QSFLOW_ACT  (1:ntiles))
       allocate (QOUTFLOW_ACT(1:ntiles))     

       LENGSC_ACT=route%lengsc/1.e3 !m->km
       AREACAT_ACT=route%areacat/1.e6 !m2->km2

       WSTREAM_ACT => route%wstream
       WRIVER_ACT => route%wriver

       !---check water balance------    
       IF(1==0)THEN  
         allocate(WTOT_BEFORE(ntiles),WTOT_AFTER(ntiles),UNBALANCE(ntiles),UNBALANCE_GLOBAL(n_catg))
         allocate(QFLOW_SINK(ntiles),QFLOW_SINK_GLOBAL(n_catg),WTOT_BEFORE_GLOBAL(n_catg),WTOT_AFTER_GLOBAL(n_catg))
         allocate(runoff_save_m3(nt_local),runoff_global_m3(nt_global),ERROR(ntiles),ERROR_GLOBAL(n_catg))
         WTOT_BEFORE=WSTREAM_ACT+WRIVER_ACT
       ENDIF
       !----------------------------

       ! Call river_routing_model
       ! ------------------------     
       CALL RIVER_ROUTING  (ntiles, RUNOFF_ACT,AREACAT_ACT,LENGSC_ACT,  &
            WSTREAM_ACT,WRIVER_ACT, QSFLOW_ACT,QOUTFLOW_ACT) 

       allocate(QOUTFLOW_GLOBAL(n_catg))
       call MPI_allgatherv  (                          &
            QOUTFLOW_ACT,  route%scounts_cat(mype+1)      ,MPI_REAL, &
            QOUTFLOW_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
            MPI_COMM_WORLD, mpierr) 

       allocate(QINFLOW_LOCAL(ntiles))
       QINFLOW_LOCAL=0.
       do i=1,nTiles
         do j=1,upmax
           if(route%upid(j,i)>0)then
             upid=route%upid(j,i)
             WRIVER_ACT(i)=WRIVER_ACT(i)+QOUTFLOW_GLOBAL(upid)*real(route_dt)
             QINFLOW_LOCAL(i)=QINFLOW_LOCAL(i)+QOUTFLOW_GLOBAL(upid)
           else
             exit
           endif
         enddo
       enddo

      !---check water balance------
       IF(1==0)THEN

         WTOT_AFTER=WRIVER_ACT+WSTREAM_ACT
         ERROR = WTOT_AFTER - (WTOT_BEFORE + RUNOFF_ACT*route_dt + QINFLOW_LOCAL*route_dt - QOUTFLOW_ACT*route_dt)
         where(QOUTFLOW_ACT>0.) UNBALANCE = abs(ERROR)/(QOUTFLOW_ACT*route_dt)
         where(QOUTFLOW_ACT<=0.) UNBALANCE = 0.
         call MPI_allgatherv  (                          &
              UNBALANCE,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              UNBALANCE_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)           
         QFLOW_SINK=0.
         do i=1,nTiles
           if(route%downid(i)==-1)then
              QFLOW_SINK(i) = QOUTFLOW_ACT(i)
           endif
         enddo
         call MPI_allgatherv  (                          &
              QFLOW_SINK,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              QFLOW_SINK_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)
         call MPI_allgatherv  (                          &
              WTOT_BEFORE,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              WTOT_BEFORE_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)   
         call MPI_allgatherv  (                          &
              WTOT_AFTER,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              WTOT_AFTER_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr) 
         runoff_save_m3=runoff_save*route%tile_area/1000. 
         call MPI_allgatherv  (                          &
              runoff_save_m3,  route%scounts_global(mype+1)      ,MPI_REAL, &
              runoff_global_m3, route%scounts_global, route%rdispls_global,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)        
         if(mapl_am_I_root())then 
           open(88,file="WTOT_AFTER.txt",status="unknown", position="append")
           write(88,*)sum(WTOT_AFTER_GLOBAL)
           close(88)
           open(88,file="WTOT_BEFORE_RUNOFF_QSINK.txt",status="unknown", position="append")
           write(88,*) sum(WTOT_BEFORE_GLOBAL)+sum(runoff_global_m3)*route_dt-sum(QFLOW_SINK_GLOBAL)*route_dt
           close(88)  
           open(88,file="WTOT_ERROR_2_RUNOFF.txt",status="unknown", position="append")
           write(88,*) (sum(WTOT_AFTER_GLOBAL)-(sum(WTOT_BEFORE_GLOBAL)+sum(runoff_global_m3)*route_dt-sum(QFLOW_SINK_GLOBAL)*route_dt))/(sum(runoff_global_m3)*route_dt)
           close(88)              
         endif                     

         call MPI_allgatherv  (                          &
              ERROR,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              ERROR_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)
         temp = maxloc(abs(ERROR_GLOBAL))
         cid = temp(1)
         if(cid>=route%minCatch.and.cid<=route%maxCatch)then
           tid=cid-route%minCatch+1
           print *,"my PE is:",mype,", max abs value of ERROR=", ERROR(tid)," at pfafid: ",route%minCatch+tid-1,", W_BEFORE=",WTOT_BEFORE(tid),", RUNOFF=",RUNOFF_ACT(tid)*route_dt,", QINFLOW=",QINFLOW_LOCAL(tid)*route_dt,", QOUTFLOW=",QOUTFLOW_ACT(tid)*route_dt,", W_AFTER=",WTOT_AFTER(tid)
         endif  
         if(FirstTime)then     
           if(mapl_am_I_root())then  
             open(88,file="ERROR_TOTAL.txt",action="write")
             do i=1,n_catg
                write(88,*)ERROR_GLOBAL(i)
             enddo
           endif
         endif
    
         if(mapl_am_I_root())print *, "The clock's final current time is ", YY, "/", MM, "/", DD, " ", HH, ":", MMM, ":", SS

         deallocate(WTOT_BEFORE,WTOT_AFTER,UNBALANCE,UNBALANCE_GLOBAL,ERROR,QFLOW_SINK,QFLOW_SINK_GLOBAL,WTOT_BEFORE_GLOBAL,WTOT_AFTER_GLOBAL)
         deallocate(runoff_save_m3,runoff_global_m3,ERROR_GLOBAL)

       ENDIF 
      !----------------------------

       deallocate(RUNOFF_ACT,AREACAT_ACT,LENGSC_ACT,QOUTFLOW_ACT,QINFLOW_LOCAL)

       ! output
       if(HH==0)then
         allocate(wriver_global(n_catg),wstream_global(n_catg),qsflow_global(n_catg))       
         call MPI_allgatherv  (                          &
              WSTREAM_ACT,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              wstream_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)
         call MPI_allgatherv  (                          &
              WRIVER_ACT,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              wriver_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)       
         call MPI_allgatherv  (                          &
              QSFLOW_ACT,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              qsflow_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)         
         if(mapl_am_I_root())then
              if(FirstTime) call system("mkdir -p ../river/river_storage ../river/stream_storage ../river/river_flow ../river/stream_flow", istat)
              write(yr_s,'(I4.4)')YY
              write(mon_s,'(I2.2)')MM
              write(day_s,'(I2.2)')DD        
              open(88,file="../river/river_storage/river_storage_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")
              open(89,file="../river/stream_storage/stream_storage_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")
              open(90,file="../river/river_flow/river_flow_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")              
              open(91,file="../river/stream_flow/stream_flow_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")
              do i=1,n_catg
                write(88,*)wriver_global(i)
                write(89,*)wstream_global(i)
                write(90,*)QOUTFLOW_GLOBAL(i)
                write(91,*)qsflow_global(i)
              enddo
              close(88);close(89);close(90);close(91)
         endif           
         deallocate(wriver_global,wstream_global,qsflow_global)
       endif

       deallocate(QOUTFLOW_GLOBAL,QSFLOW_ACT)

    ! initialize the cycle counter and sum (runoff_tile)       
       WSTREAM_ACT=>NULL()
       WRIVER_ACT=>NULL()      

       runoff_save = 0.
       ThisCycle   = 1           

       if(FirstTime) FirstTime=.False.

    else
       
       runoff_save = runoff_save + RUNOFF_SRC0/real (N_CYC)
       
       ThisCycle = ThisCycle + 1

    endif RUN_MODEL 

    runoff_save => NULL()

! All done
! --------
    call MAPL_TimerOff ( MAPL, "-RRM" ) 
    call MAPL_TimerOff(MAPL,"RUN2")
    !call MPI_Barrier(MPI_COMM_WORLD, mpierr)


    RETURN_(ESMF_SUCCESS)
  end subroutine RUN2

! ---------------------------------------------------------------------------

  subroutine InitializeRiverRouting(MYPE, numprocs, root_proc,           &
       pfaf_code, AllActive, AlldstCatchID, srcProcsID, LocDstCatchID, rc)
    
    implicit none
    INTEGER, INTENT (IN)                             :: MYPE, numprocs
    LOGICAL, INTENT (IN)                             :: root_proc
    INTEGER, DIMENSION (:),  INTENT (IN)             :: pfaf_code
    INTEGER, DIMENSION (N_CatG),          INTENT (INOUT) :: srcProcsID, LocDstCatchID
    INTEGER, DIMENSION (N_CatG,numprocs), INTENT (INOUT) :: Allactive,  AlldstCatchID

    INTEGER, DIMENSION(:)  ,ALLOCATABLE  :: global_buff, scounts, rdispls, rcounts, LocalActive
    INTEGER                              :: N_active, I,J,K,N,i1,i2,NProcs, Local_Min, Local_Max

    integer, optional,    intent(OUT):: rc

    integer :: mpierr
    character(len=ESMF_MAXSTR), parameter :: Iam='InitializeRiverRouting'

    ! STEP 1: Identify active catchments within the local processor. If the catchment is active in 
    !         more than 1 processor, choose an owner.
    ! --------------------------------------------------------------------------------------------

    allocate (LocalActive (1:N_CatG))
    LocalActive = -9999
 
    Local_Min = minval (pfaf_code)
    Local_Max = maxval (pfaf_code)
 
    do N = 1, size (pfaf_code)               
       LocalActive(pfaf_code(n)) = pfaf_code(n) 
    end do 

    allocate (global_buff (N_CatG * numprocs))
    allocate (scounts(numprocs),rdispls(numprocs),rcounts(numprocs))  

    scounts = N_CatG
    rcounts = N_CatG
    
    rdispls(1) = 0
    global_buff= 0
    
    do i=2,numprocs
       rdispls(i)=rdispls(i-1)+rcounts(i-1)
    enddo   
    
    call MPI_allgatherv  (                          &
         LocalActive, scounts         ,MPI_INTEGER, &
         global_buff, rcounts, rdispls,MPI_INTEGER, &
         MPI_COMM_WORLD, mpierr) 
    
    do i=1,numprocs
       Allactive (:,i) = global_buff((i-1)*N_CatG+1:i*N_CatG)
    enddo

    if (root_proc) then

       DO N = 1, N_CatG
          NPROCS = count(Allactive(N,:) >= 1)
          if(NPROCS > 0)then
             if (NPROCS == 1) then
                srcProcsID (N) = maxloc(Allactive(N,:),dim=1) - 1
             else
                i1 = MAX(N - 5,1)
                i2 = MIN(N + 5, N_CatG)
                N_active = 0
                do I = 1,numprocs
                   if(Allactive (N,I) >= 1) then
                      if(count (Allactive(I1:I2,I) > 0) > N_active) then
                         N_active = count (Allactive(I1:I2,I) > 0)
                         J        = I
                      endif
                   endif
                end do
                srcProcsID (N) = J - 1
             endif
          endif
       END DO

    endif

    call MPI_BCAST (srcProcsID, N_CatG, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)

    ! STEP 2: reset downstream catchment indeces (from -1 OR 1:291284) of catchments that are
    !            in the local processor to full domain indeces.
    ! ------------------------------------------------------------------------------------------

    do N = Local_Min, Local_Max
       
       if(LocalActive (N) >=1) then 
          
          if (LocDstCatchID (N) == -1) then
             ! (a) DNST Catch is a sink hole, ocean or lake so water drains to self 
             LocDstCatchID (N)    = N 
             
          endif
          
       else
          
          LocDstCatchID (N) = -9999 ! is inactive
          
       endif
    end do

    global_buff= 0
   
    call MPI_allgatherv  (                          &
         LocDstCatchID,  scounts      ,MPI_INTEGER, &
         global_buff, rcounts, rdispls,MPI_INTEGER, &
         MPI_COMM_WORLD, mpierr)
    
    do i=1,numprocs
       AlldstCatchID (:,i) = global_buff((i-1)*N_CatG+1:i*N_CatG)
    enddo
    
    deallocate (global_buff, scounts, rdispls, rcounts, LocalActive)

    RETURN_(ESMF_SUCCESS)
  end subroutine InitializeRiverRouting
end module GEOS_RouteGridCompMod
