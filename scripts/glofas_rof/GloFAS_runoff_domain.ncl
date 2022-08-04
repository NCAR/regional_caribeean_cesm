;===============================================================================
; create data model domain file from JRA-55 TL319 data file and add a mask
; to represent the domain of the caribeean regional example
;===============================================================================

newline = integertochar(10) ; create a newline character

;----- hard-coded input --------------------------------------------------------

; path_in="/glade/p/cesm/omwg_dev/JRA-55/"
path_in ="./"
fn_in   = "glofas_runoff_masked.nc"
fn_out  = "domain.roff.GLOFAS025m_carib.220520.nc"
fn_out2 = "GLOFAS025m_carib.220520.nc"

title       = "GloFAS runoff grid AND MASK -- for use with GloFAS runoff data in the Caribbean regional configuration"
history     = "Created by " + systemfunc("whoami") + " " + systemfunc("date")
notes       = "GloFAS runoff domain data, formatted for use in CESM"
; notes       = notes + newline + "Includes sort field where 1=>largest runoff value, 2=> 2nd largest runoff value,etc."
Conventions = "CF1.0"
Source      = "Copernicus EU GloFAS "
; Source      = Source + newline + "Tsujino et al. 2018, Ocean Modell., 130, 79-139 "
; Source      = Source + newline + "Suzuki et al. 2017. J. Oceanogr., 74, 421-429"

masked = 1

;-------------------------------------------------------------------------------


ni = 750 ; (long) WARNING: hard-coded 0.25 degree
nj =  410 ; (lat) WARNING: hard-coded 0.25 degree
nv =    4 ; WARNING: hard-coded four verticies

pi = 4.0d * atan(1.0d)
deg2rad = pi/180.0d

dx      = 0.25d    ; WARNING: hard-coded for 0.25 degree
half_dx = 0.125d

;-------------------------------------------------------------------------------
print_clock("reading "+fn_in)
;-------------------------------------------------------------------------------

fh_in = addfile(path_in+fn_in,"r")

longitude=fh_in->longitude
latitude =fh_in->latitude
dmask2 =fh_in->mask

;tmp1 =fh_in->rofl(0,:,:)
;tmp2 =fh_in->rofi(0,:,:)
runoff = fh_in->runoff(:,:) ;tmp1 + tmp2
;delete([/tmp1,tmp2/])

;-------------------------------------------------------------------------------
print_clock("create data")
;-------------------------------------------------------------------------------

xc    = new( (/nj,ni/)   , "double")
yc    = new( (/nj,ni/)   , "double")
xv    = new( (/nj,ni,nv/), "double")
yv    = new( (/nj,ni,nv/), "double")
area  = new( (/nj,ni/)   , "double")
dmask = new( (/nj,ni/)   , "integer")
sort  = new( (/nj,ni/)   , "integer")
fracl = new( (/nj,ni/)   , "double")


xc    = -999.
yc    = -999.
xv    = -999.
yv    = -999.
area  = -999.
dmask = 0
order = -999
fracl = -999.

;----- create mask data -----
do j=0,nj-1
   do i=0,ni-1
     if (dmask2(j,i) .gt. 0) then
       dmask(j,i) = 1
     end if
   end do
end do

;dmask = 1
;if (masked .ne.0) then
;   dmask = where( runoff > 0.0, 1, 0)
;end if

;----- create verticie data -----
    j=0       ; next to south pole
    yc(j,:)   = latitude(j)
    yv(j,:,0) = -90.0d
    yv(j,:,1) = -90.0d
    yv(j,:,2) = (latitude(j)+latitude(j+1))/2.0d
    yv(j,:,3) = (latitude(j)+latitude(j+1))/2.0d

    j=nj-1    ; next to north pole
    yc(j,:)   =  latitude(j)
    yv(j,:,0) = (latitude(j)+latitude(j-1))/2.0d
    yv(j,:,1) = (latitude(j)+latitude(j-1))/2.0d
    yv(j,:,2) =  90.0d
    yv(j,:,3) =  90.0d

do j=1,nj-2
    yc(j,:)   =  latitude(j)
    yv(j,:,0) = (latitude(j)+latitude(j-1))/2.0d
    yv(j,:,1) = (latitude(j)+latitude(j-1))/2.0d
    yv(j,:,2) = (latitude(j)+latitude(j+1))/2.0d
    yv(j,:,3) = (latitude(j)+latitude(j+1))/2.0d
end do

do i=0,ni-1
   xc(:,i)   = longitude(i)
   xv(:,i,0) = longitude(i) - half_dx
   xv(:,i,1) = longitude(i) + half_dx
   xv(:,i,2) = longitude(i) + half_dx
   xv(:,i,3) = longitude(i) - half_dx
end do

;----- create area data using verticie data -----

sum_test = 0.0d
do j=0,nj-1
   y0 = sin(yv(j,0,0)*deg2rad) ; south
   y1 = sin(yv(j,0,3)*deg2rad) ; north
   do i=0,ni-1
      x0 = xv(j,i,0) ; west
      x1 = xv(j,i,1) ; east
      area(j,i) = (y1-y0)*(x1-x0)*deg2rad
      sum_test = sum_test+area(j,i)
   end do
end do
print_clock("sum_test/(4*pi) = "+sum_test/(4.0d*pi)  )

;----- create sort data, value 1 => largest runoff, 2 => 2nd largest runoff -----
print_clock("create sort field")

print("test 2d-1d index converions START")
n = -1
do j=0,nj-1
do i=0,ni-1
   n=n+1
   if (i .ne. mod(n,ni)) then
      print("i index error")
   end if
   if (j .ne. toint(floor(n/ni)) ) then
      print("j index error")
   end if
end do
end do
print("test 2d-1d index converions END")


runoff1d = new( (/nj*ni/)   , "double")

n=-1
do j=0,nj-1
do i=0,ni-1
   n=n+1
   runoff1d(n) = runoff(j,i) * area(j,i)
end do
end do

sortindex = dim_pqsort_n(runoff1d,-1,0)

n=-1
do j=0,nj-1
do i=0,ni-1
   n=n+1
   sort(j,i) = sortindex(n)
end do
end do

if (masked .ne. 0) then
   sort = where( runoff > 0.0, sort, -999)
end if

;rint("123456 123456 123456 123456 12345678901234 123456789 123456789
print("sort k      n      i      j    runoff*area       lon       lat")
do k=0,19
   n = sortindex(k)
   i = mod(n,ni)
   j = toint(floor(n/ni))
   print(sprinti("%6.1i",k) + " " \
       + sprinti("%6.1i",n) + " " \
       + sprinti("%6.1i",i) + " " \
       + sprinti("%6.1i",j) + " " \
       + sprintf("%14.11f",runoff(j,i)*area(j,i)) + " " \
       + sprintf("%9.2f",xc(j,i)) + " " \
       + sprintf("%9.2f",yc(j,i)) )
end do


print("number of grid cells                : " + ni*nj)
print("number of non-zero runoff cells     : " + num(runoff .gt. 0))
print("number of non-negative sort indicies: " + num(sort .ge. 0))
print("number of unique       sort indicies: " + count_unique_values_n(sortindex,0))
print("number of unmasked cells            : " + num(dmask .ne. 0))

;=========================================================================================
; make a CESM data model domain file (fields are 2d)
;=========================================================================================

;-------------------------------------------------------------------------------
print_clock("create meta-data")
;-------------------------------------------------------------------------------

xc@bounds       = "xv"
xc@units        = "degrees_north"
xc@long_name    = "longitude of grid cell center"
xc@cartesian_axis = "X"
xc@axis = "X"
xc@point_spacing = "even"
delete(xc@_FillValue)
delete(xc@standard_name)

yc@bounds       = "yv"
yc@units        = "degrees_north"
yc@long_name    = "latitude of grid cell center"
yc@cartesian_axis = "Y"
yc@axis = "Y"
yc@point_spacing = "even"
delete(yc@_FillValue)
delete(yc@standard_name)

xv@units        = "degrees_north"
xv@long_name    = "longitude of grid cell verticies"
delete(xv@_FillValue)

yv@units        = "degrees_north"
yv@long_name    = "latitude of grid cell verticies"
delete(yv@_FillValue)

dmask@coordinates = "yc xc"
dmask@comment     = "0 value if cell is not active"
dmask@note        = "unitless"
dmask@long_name   = "domain mask"
delete(dmask@_FillValue)

sort@coordinates = "yc xc"
sort@comment     = "1=> largest runoff value, 2=> 2nd largest runoff value, etc"
sort@note        = "unitless"
sort@long_name   = "index of runoff magnitude when sorted"
delete(sort@_FillValue)

area@coordinates = "yc xc"
area@units       = "radian2"
area@long_name   = "area of grid cell in radians squared"
delete(area@_FillValue)

fracl@coordinates = "yc xc"
fracl@note        = "unitless, range = [0,1]"
fracl@long_name   = "fraction of land in cell"
delete(fracl@_FillValue)

gAtt = True
gAtt@notes       = notes
gAtt@history     = history
gAtt@Conventions = Conventions
gAtt@Source      = Source
gAtt@title       = title

;-------------------------------------------------------------------------------
print_clock("writing "+fn_out)
;-------------------------------------------------------------------------------

system("/bin/rm -f "+fn_out)   ; clobber any previously created file
fh_out = addfile(fn_out ,"c")  ; open output netCDF file

setfileoption(fh_out,"DefineMode",True) ;---------------------------------------

dimNames = (/"n"    ,"ni"   ,"nj", "nv"  /)
dimSizes = (/ ni*nj , ni    , nj  , nv   /)
dimUnlim = (/False  ,False, False, False /)
filedimdef(fh_out,dimNames,dimSizes,dimUnlim)

fileattdef    (fh_out,gAtt)           ; add global attributes

filevardef(fh_out,"xc"   ,typeof(xc)   , (/"nj","ni"/) )
filevardef(fh_out,"yc"   ,typeof(yc)   , (/"nj","ni"/) )
filevardef(fh_out,"xv"   ,typeof(xv)   , (/"nj","ni","nv"/) )
filevardef(fh_out,"yv"   ,typeof(yv)   , (/"nj","ni","nv"/) )
filevardef(fh_out,"mask" ,typeof(dmask), (/"nj","ni"/) )
; filevardef(fh_out,"sort" ,typeof(sort) , (/"nj","ni"/) )
filevardef(fh_out,"area" ,typeof(area) , (/"nj","ni"/) )
; filevardef(fh_out,"fracl",typeof(fracl), (/"nj","ni"/) )

filevarattdef(fh_out,"xc"   ,xc)
filevarattdef(fh_out,"yc"   ,yc)
filevarattdef(fh_out,"xv"   ,xv)
filevarattdef(fh_out,"yv"   ,yv)
filevarattdef(fh_out,"mask" ,dmask)
; filevarattdef(fh_out,"sort" ,sort)
filevarattdef(fh_out,"area" ,area)
;filevarattdef(fh_out,"fracl",fracl)

setfileoption(fh_out,"DefineMode",False) ;---------------------------------------

fh_out->xc    = (/xc    /)
fh_out->yc    = (/yc    /)
fh_out->xv    = (/xv    /)
fh_out->yv    = (/yv    /)
fh_out->mask  = (/dmask /)
; fh_out->sort  = (/sort  /)
fh_out->area  = (/area  /)
;fh_out->fracl = (/fracl /)

;=========================================================================================
;=========================================================================================
; make scrip/esmf domain file (fields are 1d vectors)
;=========================================================================================

;-------------------------------------------------------------------------------
print_clock("create data")
;-------------------------------------------------------------------------------
grid_size    = ni*nj
grid_rank    = 2
grid_corners = 4

s_yc   = new( (/grid_size/)             , "double")
s_xc   = new( (/grid_size/)             , "double")
s_yv   = new( (/grid_size,grid_corners/), "double")
s_xv   = new( (/grid_size,grid_corners/), "double")
s_area = new( (/grid_size/)             , "double")
s_mask = new( (/grid_size/)             , "integer")
s_dims = new( (/grid_rank/)             , "integer")

s_dims(0) = ni
s_dims(1) = nj

n=-1
do j=0,nj-1
do i=0,ni-1
   n=n+1
   s_yc  (n)   = yc   (j,i)
   s_yv  (n,0) = yv   (j,i,0)
   s_yv  (n,1) = yv   (j,i,1)
   s_yv  (n,2) = yv   (j,i,2)
   s_yv  (n,3) = yv   (j,i,3)
   s_xc  (n)   = xc   (j,i)
   s_xv  (n,0) = xv   (j,i,0)
   s_xv  (n,1) = xv   (j,i,1)
   s_xv  (n,2) = xv   (j,i,2)
   s_xv  (n,3) = xv   (j,i,3)
   s_area(n)   = area (j,i)
   s_mask(n)   = dmask(j,i)
end do
end do

;-------------------------------------------------------------------------------
print_clock("create meta-data")
;-------------------------------------------------------------------------------

gAtt@Conventions = "CESM/SCRIP/ESMF mapping file domain specification"

s_yc@units = "degrees"
s_xc@units = "degrees"
s_yv@units = "degrees"
s_xv@units = "degrees"

s_area@units= "radians^2"
s_mask@note = "unitless"
s_mask@convention = "0 <=> inactive cell"
s_dims@note = "dimensions of underlying multi-dimensional array"

delete(s_dims@_FillValue)
delete(s_yc@_FillValue)
delete(s_xc@_FillValue)

delete(s_yv@_FillValue)
delete(s_xv@_FillValue)

delete(s_xc@bounds)
delete(s_yc@bounds)
delete(s_mask@_FillValue)
delete(s_mask@coordinates)
delete(s_area@_FillValue)
delete(s_area@coordinates)

;----- esmf doesn't like these??? ----
; delete(s_dims@note)
; delete(s_xc@standard_name)
; delete(s_yc@standard_name)
; delete(s_xv@long_name)
; delete(s_yv@long_name)
; delete(s_area@long_name)
; delete(s_area@units)
; delete(s_mask@comment)
; delete(s_mask@note)
; delete(s_mask@long_name)
; s_mask@units = "unitless"
; delete(gAtt@Source)

;-------------------------------------------------------------------------------
print_clock("writing "+fn_out2)
;-------------------------------------------------------------------------------

system("/bin/rm -f "+fn_out2)   ; clobber any previously created file
fh_out = addfile(fn_out2 ,"c")  ; open output netCDF file

setfileoption(fh_out,"DefineMode",True) ;---------------------------------------

delete(dimNames)
delete(dimSizes)
delete(dimUnlim)
dimNames = (/"grid_size","grid_rank","grid_corners" /)
dimSizes = (/ grid_size , grid_rank , grid_corners  /)
dimUnlim = (/ False     , False     , False         /)
filedimdef(fh_out,dimNames,dimSizes,dimUnlim)

fileattdef    (fh_out,gAtt)           ; add global attributes

filevardef(fh_out,"grid_dims"       ,typeof(s_dims) , (/"grid_rank"/) )
filevardef(fh_out,"grid_center_lon" ,typeof(s_xc)   , (/"grid_size"/) )
filevardef(fh_out,"grid_center_lat" ,typeof(s_yc)   , (/"grid_size"/) )
filevardef(fh_out,"grid_corner_lon" ,typeof(s_xv)   , (/"grid_size","grid_corners"/) )
filevardef(fh_out,"grid_corner_lat" ,typeof(s_yv)   , (/"grid_size","grid_corners"/) )
filevardef(fh_out,"grid_imask"      ,typeof(s_mask) , (/"grid_size"/) )
filevardef(fh_out,"grid_area"       ,typeof(s_area) , (/"grid_size"/) )

filevarattdef(fh_out,"grid_center_lon" ,s_xc   )
filevarattdef(fh_out,"grid_center_lat" ,s_yc   )
filevarattdef(fh_out,"grid_corner_lon" ,s_xv   )
filevarattdef(fh_out,"grid_corner_lat" ,s_yv   )
filevarattdef(fh_out,"grid_imask"      ,s_mask )
filevarattdef(fh_out,"grid_area"       ,s_area )
filevarattdef(fh_out,"grid_dims"       ,s_dims )

setfileoption(fh_out,"DefineMode",False) ;---------------------------------------

fh_out->grid_dims       = (/s_dims  /)
fh_out->grid_corner_lon = (/s_xv    /)
fh_out->grid_corner_lat = (/s_yv    /)
fh_out->grid_center_lon = (/s_xc    /)
fh_out->grid_center_lat = (/s_yc    /)
fh_out->grid_imask      = (/s_mask  /)
fh_out->grid_area       = (/s_area  /)

;=========================================================================================
print_clock("done. ")
;=========================================================================================

quit

