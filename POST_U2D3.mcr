#!MC 1300
# Created by Tecplot 360 build 13.1.1.16309
$!VarSet |MFBD| = 'D:\post\Pitch_test'

$!EXPORTSETUP EXPORTFORMAT = AVI
$!EXPORTSETUP IMAGEWIDTH = 878
$!EXPORTSETUP EXPORTFNAME = 'D:\post\post_result\output_post_2d3.avi'

## Begin Animating
$!VarSet |NumLoop|=150
$!Loop |NumLoop|
$!VarSet |Num|=(|Loop|*100)

$!READDATASET  '"|MFBD|\POST_U_2D3_0001\POST_U_2D3_|Num%010d|_0001.DAT" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"X" "Y" "Z" "U" "V" "W"'
$!PICK ADDATPOSITION
  X = 1.54609475032
  Y = 4.29609475032
  CONSIDERSTYLE = YES
$!TWODAXIS GRIDAREA{DRAWBORDER = YES}
$!TWODAXIS YDETAIL{TICKS{SHOWONGRIDBORDERMAX = YES}}
$!TWODAXIS XDETAIL{TICKS{SHOWONGRIDBORDERMAX = YES}}
$!PICK ADDATPOSITION
  X = 1.54609475032
  Y = 4.31658130602
  CONSIDERSTYLE = YES
$!TWODAXIS DEPXTOYRATIO = 0.25
$!TWODAXIS YDETAIL{VARNUM = 3}
$!TRIANGULATE 
  SOURCEZONES =  [1]
  USEBOUNDARY = NO
  INCLUDEBOUNDARYPTS = NO
  TRIANGLEKEEPFACTOR = 0.25
$!GLOBALCONTOUR 1  VAR = 6
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!FIELDLAYERS SHOWCONTOUR = YES
$!GLOBALCONTOUR 1  VAR = 4
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15

$!TWODAXIS YDETAIL{RANGEMIN = -10}
$!TWODAXIS YDETAIL{RANGEMAX = 100}
$!TWODAXIS XDETAIL{RANGEMIN = 10}
$!TWODAXIS XDETAIL{RANGEMAX = 500}
$!GLOBALCONTOUR 1  LEGEND{SHOW = YES}
$!GLOBALCONTOUR 1  LEGEND{AUTORESIZE = YES}
$!GLOBALCONTOUR 1  LEGEND{OVERLAYBARGRID = NO}

  $!IF |Loop| == 1
    $!EXPORTSTART
    EXPORTREGION = CURRENTFRAME
  $!ENDIF
  $!IF |Loop| != 1
    $!EXPORTNEXTFRAME
  $!ENDIF

$!EndLoop

$!EXPORTFINISH 

$!RemoveVar |MFBD|
