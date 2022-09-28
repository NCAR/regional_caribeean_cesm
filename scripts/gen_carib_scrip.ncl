; gen_scrip.ncl :

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
setvalues NhlGetWorkspaceObjectId()
        "wsMaximumSize" : 100000000000
end setvalues

begin
;---Input files
    ncGrdFilePath = "CARIB_025_grid.nc"

;---Output (and input) files
    dstGridPath = "CARIB_025_grid_SCRIP.nc"

;----------------------------------------------------------------------
; Convert destination grid to a SCRIP convention file.
;----------------------------------------------------------------------
    mom6_grd_file = addfile(ncGrdFilePath,"r")
    tlon     = mom6_grd_file->tlon(:,:)
    tlat     = mom6_grd_file->tlat(:,:)
    qlon     = mom6_grd_file->qlon(:,:)
    qlat     = mom6_grd_file->qlat(:,:)

    Opt                = True
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@NetCDFType     = "netcdf4"
    Opt@Title          = "Caribbean 1/4 deg. Grid"

    ; Corners' lat lon --------------------------------------------------

    grid_dimsizes = getfiledimsizes(mom6_grd_file)
    nyp = grid_dimsizes(0)
    nxp = grid_dimsizes(1)
    ny = grid_dimsizes(2)
    nx = grid_dimsizes(3)
    n = nx*ny

    GridCornerLat = new( (/ny, nx, 4/), "double")
    GridCornerLon = new( (/ny, nx, 4/), "double")

    ;----- create verticie data -----
    do j=0,ny-1
        GridCornerLat(j,:,0) = mom6_grd_file->qlat(j,0)
        GridCornerLat(j,:,1) = mom6_grd_file->qlat(j,0)
        GridCornerLat(j,:,2) = mom6_grd_file->qlat(j+1,0)
        GridCornerLat(j,:,3) = mom6_grd_file->qlat(j+1,0)
    end do

    do i=0,nx-1
        GridCornerLon(:,i,0)   = mom6_grd_file->qlon(0,i)
        GridCornerLon(:,i,1) = mom6_grd_file->qlon(0,i+1)
        GridCornerLon(:,i,2) = mom6_grd_file->qlon(0,i+1)
        GridCornerLon(:,i,3) = mom6_grd_file->qlon(0,i)
    end do

    ; Mask --------------------------------------------------

    maskt = mom6_grd_file->tmask(:,:)
    Opt@GridMask = toint(maskt)
    Opt@GridCornerLat = GridCornerLat
    Opt@GridCornerLon = GridCornerLon

    ; Generate the script file ------------------------------

    ;curvilinear_to_SCRIP(dstGridPath,tlat,tlon,Opt)
    rectilinear_to_SCRIP(dstGridPath,tlat(:,0),tlon(0,:),Opt)

;---Clean up
    delete(Opt)

    ; Append Area to SCRIP file------------------------------

    scripFile = addfile(dstGridPath,"w")

    grid_size = dimsizes(scripFile->grid_center_lat)
    grid_area = new(grid_size,double)
    grid_area!0 = "grid_size"

    do i=0,grid_size-1
      temp_tlat = (/ scripFile->grid_corner_lat(i,3), \
                scripFile->grid_corner_lat(i,1), \
                scripFile->grid_corner_lat(i,0), \
                scripFile->grid_corner_lat(i,2)    /)
      temp_tlon = (/ scripFile->grid_corner_lon(i,3), \
                scripFile->grid_corner_lon(i,1), \
                scripFile->grid_corner_lon(i,0), \
                scripFile->grid_corner_lon(i,2)    /)

      grid_area(i) = area_poly_sphere(temp_tlat, temp_tlon, 1)

    end do

    scripFile->grid_area = grid_area

end
