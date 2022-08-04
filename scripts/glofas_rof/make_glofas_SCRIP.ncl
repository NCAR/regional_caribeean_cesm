load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file
    srcFileName = "/glade/work/gseijo/GloFAS/GloFAS_2000_2014_cesm.nc"

;---Output (and input) files
    dstGridName = "GloFAS_Carib_SCRIP.nc"
;----------------------------------------------------------------------
; Convert original source NCEP grid to a SCRIP convention file.
;----------------------------------------------------------------------
    src_file = addfile(srcFileName,"r")
    rof     = src_file->runoff(0,:,:)

    Opt                = True
;original code->    Opt@GridMask       = where(.not.ismissing(rof),1,0)    ; use "Mask2D" in NCL V6.2.1 and earlier
    Opt@GridMask = src_file->land_ocean_mask(:,:)
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "NCEP Grid"

    rectilinear_to_SCRIP(dstGridName,rof&latitude,rof&longitude,Opt)

;---Clean up
    delete(Opt)

end