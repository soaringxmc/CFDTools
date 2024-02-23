#!MC 900
$!VarSet |MFBD| = '.'
$!DRAWGRAPHICS TRUE

#Ñ­»·´ÎÊý
$!VarSet |n| = 85
$!VarSet |dt| = 10000
$!VarSet |t|= 490000
$!VarSet |file|="airfoil"



$!VarSet |filename|="|file|.0|t|"


$!LOOP |n|
  $!DRAWGRAPHICS FALSE
     $!IF |t| <100000
					$!VarSet |filename|="|file|.0|t|"
     $!ENDIF
     $!IF |t| >=100000
					$!VarSet |filename|="|file|.|t|"
     $!ENDIF
#####################################################
     $!READDATASET  '"|MFBD|\|filename|.plt" '
		  READDATAOPTION = NEW
		  RESETSTYLE = YES
		  INCLUDETEXT = NO
		  INCLUDEGEOM = NO
		  INCLUDECUSTOMLABELS = NO
		  VARLOADMODE = BYNAME
		  ASSIGNSTRANDIDS = YES
		  INITIALPLOTTYPE = CARTESIAN3D
		  VARNAMELIST = '"X" "Y" "Z" "RHO" "P" "T" "TAU_WALL" "YPLUS" "YPLUS_WM" "XPLUS" "DIL" "P_AVG" "P_RMS" "U-X" "U-Y" "U-Z"'
		$!THREEDVIEW 
		  PSIANGLE = 60
		  THETAANGLE = -44.1355
		  ALPHAANGLE = -0
		  VIEWERPOSITION
		    {
		    X = 0.5469368589740409
		    Y = -0.6313975559653523
		    Z = 0.5222771398384771
		    }
		$!VIEW PUSH
		$!GLOBALCONTOUR 1  VAR = 4
		$!CONTOURLEVELS RESETTONICE
		  CONTOURGROUP = 1
		  APPROXNUMVALUES = 15
		$!FIELDLAYERS SHOWCONTOUR = YES
		$!WRITEDATASET  "|MFBD|\|filename|.dat"
		  INCLUDETEXT = NO
		  INCLUDEGEOM = NO
		  INCLUDECUSTOMLABELS = NO
		  ASSOCIATELAYOUTWITHDATAFILE = NO
		  BINARY = NO
		  USEPOINTFORMAT = YES
		  PRECISION = 15
		  TECPLOTVERSIONTOWRITE = TECPLOTCURRENT

#####################################################     
   $!VarSet |t| += |dt|
$!ENDLOOP

$!Quit

$!EXPORTFINISH
$!RemoveVar |t|
$!RemoveVar |file|
$!RemoveVar |filename|
$!RemoveVar |dt|
$!RemoveVar |n|
$!RemoveVar |MFBD|
