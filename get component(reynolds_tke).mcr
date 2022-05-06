#!MC 1300
# Created by Tecplot 360 build 13.1.1.16309
## Set up  Export file type and file name.

$!VarSet |MFBD| = 'D:\post\Pitch_test'

$!VarSet |NumLoop|=25
$!VarSet |Num|=(|NumLoop|*200)
$!VarSet |L0|=1.0
$!VarSet |U0|=11.5258407161
$!VarSet |T0|=(|L0|/|U0|)
$!VarSet |T_turbine| = 42.84
$!VarSet |time|=((|Num|*0.68543297937)/|T_turbine|)

$!READDATASET  '"-F" "1" "|MFBD|\mean_field_3di.h5" "-D" "12" "/dsgs" "/dsgsu" "/nacfx" "/nut" "/pp" "/tke" "/u" "/uw" "/v" "/w" "/wtfx" "/z" "-I" "-K" "1" "1" "1"'
#$!READDATASET  '"-F" "1" "D:\post\Project\Fixed_Turbine\mean_tke_3d.h5" "-D" "9" "/Advection" "/Dissipation" "/Production" "/TKE" "/TransportByPressure" "/TransportByTurbulence" "/TransportByViscousStress" "/u'w'" "/z" "-I" "-K" "1" "1" "1"'
  DATASETREADER = 'HDF5 Loader'
  $!ALTERDATA 
    EQUATION = '{x1}=({Z})*6.2831853071795862/0.00280499344071/192'
  $!ALTERDATA 
    EQUATION = '{y1}=({Y})*6.2831853071795862/0.00280499344071/192'
  $!ALTERDATA 
    EQUATION = '{Z}={/z}'
  $!ALTERDATA 
    EQUATION = '{X}={x1}'
  $!ALTERDATA 
    EQUATION = '{Y}={y1}'
  #$!ALTERDATA 
  #  EQUATION = '{u1}={/u}*11.5258407161'
  #$!ALTERDATA 
  #  EQUATION = '{v1}={/v}*11.5258407161'
  #$!ALTERDATA 
  #  EQUATION = '{w1}={/w}*11.5258407161'

  $!VarSet |NumTurbine|=16
  $!VarSet |NumNac|= (|NumTurbine|*3)
  $!Loop |NumNac|

    $!READDATASET  '"|MFBD|\NacelleLocation\surfnac_|Num%06d|_|Loop%03d|_nf.dat"' 
    READDATAOPTION = APPEND
    RESETSTYLE = NO
    INCLUDETEXT = NO
    INCLUDEGEOM = NO
    INCLUDECUSTOMLABELS = NO
    VARLOADMODE = BYNAME
    ASSIGNSTRANDIDS = YES
    VARNAMELIST = '"X" "Y" "Z" "/dsgs" "/dsgsu" "/nacfx" "/nut" "/pp" "/tke" "/u" "/uw" "/v" "/w" "/wtfx" "/z" "x1" "y1" "u1" "v1" "w1"'
    #VARNAMELIST = '"X" "Y" "Z" "/Advection" "/Dissipation" "/Production" "/TKE" "/TransportByPressure" "/TransportByTurbulence" "/TransportByViscousStress" "/z" "x1" "y1"'
    #VARNAMELIST = '"X" "Y" "Z" "/Advection" "/Dissipation" "/Production" "/TKE" "/TransportByPressure" "/TransportByTurbulence" "/TransportByViscousStress" "/u\'w\'" "/z" "x1" "y1"'
  $!EndLoop
    
  $!Loop |NumTurbine|
    
    $!READDATASET  '"|MFBD|\LineLocation\line_|Num%06d|_|Loop%03d|_nf.dat"'
    READDATAOPTION = APPEND
    RESETSTYLE = NO
    INCLUDETEXT = NO
    INCLUDEGEOM = NO
    INCLUDECUSTOMLABELS = NO
    VARLOADMODE = BYNAME
    ASSIGNSTRANDIDS = YES
    VARNAMELIST = '"X" "Y" "Z" "/dsgs" "/dsgsu" "/nacfx" "/nut" "/pp" "/tke" "/u" "/uw" "/v" "/w" "/wtfx" "/z" "x1" "y1" "u1" "v1" "w1"'
    #VARNAMELIST = '"X" "Y" "Z" "/Advection" "/Dissipation" "/Production" "/TKE" "/TransportByPressure" "/TransportByTurbulence" "/TransportByViscousStress" "/z" "x1" "y1"'
    #VARNAMELIST = '"X" "Y" "Z" "/Advection" "/Dissipation" "/Production" "/TKE" "/TransportByPressure" "/TransportByTurbulence" "/TransportByViscousStress" "/u\'w\'" "/z" "x1" "y1"'  
  $!EndLoop

  $!FIELDLAYERS SHOWMESH = YES

  $!FRAMECONTROL ACTIVATEBYNUMBER
    FRAME = 1
  $!SETFRAMEBACKGROUNDCOLOR BLACK
  $!ROTATE3DVIEW THETA
    ANGLE = 180
    ROTATEORIGINLOCATION = DEFINEDORIGIN
  $!SLICELAYERS SHOW = YES
  $!SLICEATTRIBUTES 1  PRIMARYPOSITION{X = 1962}
  $!SLICEATTRIBUTES 2  SHOWGROUP = YES
  $!SLICEATTRIBUTES 2  PRIMARYPOSITION{Y = 1962}
  $!GLOBALCONTOUR 1  VAR = 11
  #$!CONTOURLEVELS RESETTONICE
  #  CONTOURGROUP = 1
  #  APPROXNUMVALUES = 15
  #$!CONTOURLEVELS NEW
  #  CONTOURGROUP = 1

  $!SLICEATTRIBUTES 3  SHOWGROUP = YES
  $!SLICEATTRIBUTES 3  SLICESURFACE = IPLANES
  $!SLICEATTRIBUTES 3  PRIMARYPOSITION{I = 1}
  $!GLOBALCOLORMAP 1  CONTOURCOLORMAP = DKRAINBOW
  $!FIELDMAP [2-65]  MESH{COLOR = BLACK}
  $!VIEW FITSURFACES

  $!THREEDVIEW VIEWERPOSITION{X = -307.97082787509316}
  $!THREEDVIEW VIEWERPOSITION{Y = -16039.580493221587}
  $!THREEDVIEW VIEWERPOSITION{Z = 787.94114748588447}
  $!THREEDVIEW PSIANGLE = 90
  $!THREEDVIEW THETAANGLE = 0

  $!VIEW FITSURFACES

  $!VIEW ZOOM
  X1 = 783.354830485
  Y1 = -146.714611687
  X2 = 868.240708854
  Y2 = -0.0935490494521
  
  $!REDRAWALL