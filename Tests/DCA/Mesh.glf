# $Id$

package require PWI_Glyph 1.6.9

# Delete any existing grids and database entities.  Reset AS/W, defaults, and
# tolerances.
gg::memClear
gg::aswDeleteBC -glob "*"
gg::aswDeleteVC -glob "*"
gg::aswSet GENERIC -dim 3
gg::defReset
gg::tolReset
# Delay screen updates and checking for user input until script is finished.
gg::updatePolicy DELAYED

gg::aswSet "CGNS-Struct" -dim 3

# Geometric Parameters
# FIXME: not passed to CalculatePoints.py
set Cax 1.0
set Pitch 1.0
set LXUP 1.0
set LXDOWN 1.0
set LZ [expr $Cax * 0.01]
set ThetaTE 30.0
set PI [expr 4.0 * atan(1.0)]
set ThetaTERad [expr $ThetaTE * $PI / 180.0]

# Mesh Parameters
set NBLADE 100
set NUP 50
set NWAKE 50
set NInflow 50
set DW 1.0e-5
set DSLE [expr $Cax * 0.002]
set DSTE $DSLE
set NZ 3

set buf [exec python ./CalculatePoints.py]

set lines [split $buf "\n"]
puts $lines
set PLE [lindex $lines 0]
set PTE [lindex $lines 1]
set RSS [lindex [lindex $lines 2] 0]
set PcSS [lrange [lindex $lines 2] 1 3]
set RPS [lindex [lindex $lines 3] 0]
set PcPS [lrange [lindex $lines 3] 1 3]
puts $PLE
puts $PTE
puts $RSS
puts $PcSS
puts $RPS
puts $PcPS

gg::dbCurveBegin -type CIRCULAR_ARC
  gg::dbCurveAddPt $PLE
  gg::dbCurveAddPt $PTE
  gg::dbCurveAddPt -alternate CENTER $PcPS
set db_PS [gg::dbCurveEnd]

gg::dbCurveBegin -type CIRCULAR_ARC
  gg::dbCurveAddPt $PLE
  gg::dbCurveAddPt $PTE
  gg::dbCurveAddPt -alternate CENTER $PcSS
set db_SS [gg::dbCurveEnd]

gg::dbTransformBegin [list $db_SS]
  gg::xformTranslate [list 0 $Pitch 0]
gg::dbTransformEnd

set conPS [gg::conOnDBEnt $db_PS]
gg::conDim $conPS $NBLADE

set conSS [gg::conOnDBEnt $db_SS]
gg::conDim $conSS $NBLADE

gg::conBegin
  gg::segBegin -type 3D_LINE
    gg::segAddControlPt [list -$LXUP 0.0 0.0]
    gg::segAddControlPt -db [list 0 0 $db_PS]
  gg::segEnd
set conUP [gg::conEnd]
gg::conDim $conUP $NUP

gg::conBegin
  gg::segBegin -type 3D_LINE
    gg::segAddControlPt -db [list 1 0 $db_PS]
    gg::segAddControlPt [list [expr [lindex $PTE 0] + $LXDOWN] [expr [lindex $PTE 1] + $LXDOWN * tan($ThetaTERad)] 0.0]
  gg::segEnd
set conWAKE [gg::conEnd]
gg::conDim $conWAKE $NWAKE

set conUPTrans [gg::conPeriodicTrans $conUP "0 $Pitch 0"]
set conWAKETrans [gg::conPeriodicTrans $conWAKE "0 $Pitch 0"]

gg::conBegin
  gg::segBegin -type 3D_LINE
    gg::segAddControlPt [gg::conGetPt $conUP -arc 0]
    gg::segAddControlPt [gg::conGetPt $conUPTrans -arc 0]
  gg::segEnd
set conInflow [gg::conEnd]
gg::conDim $conInflow $NInflow

gg::conBegin
  gg::segBegin -type 3D_LINE
    gg::segAddControlPt [gg::conGetPt $conWAKE -arc 1]
    gg::segAddControlPt [gg::conGetPt $conWAKETrans -arc 1]
  gg::segEnd
set conOutflow [gg::conEnd]
gg::conDim $conOutflow $NInflow

# Domain

gg::domBegin -type STRUCTURED
  gg::edgeBegin
    gg::edgeAddCon $conUP
    gg::edgeAddCon $conPS
    gg::edgeAddCon $conWAKE
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $conOutflow
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $conWAKETrans
    gg::edgeAddCon $conSS
    gg::edgeAddCon $conUPTrans
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $conInflow
  gg::edgeEnd
set dom1 [gg::domEnd]

gg::domEllSolverBegin $dom1
  gg::domEllSolverStep -iterations 200
gg::domEllSolverEnd

set dom2 [gg::domSplit $dom1 -i $NUP]
set dom3 [gg::domSplit $dom2 -i $NBLADE]

set cons [gg::domGetEdge $dom1]
set conLE [lindex $cons 1]
set cons [gg::domGetEdge $dom3]
set conTE [lindex $cons 3]

gg::conBeginSpacing $conSS $DSLE
gg::conEndSpacing $conSS $DSTE
gg::conBeginSpacing $conPS $DSLE
gg::conEndSpacing $conPS $DSTE
gg::conEndSpacing $conUP $DSLE
gg::conBeginSpacing $conWAKE $DSTE

gg::conDistFunc $conLE -function TANH
gg::conBeginSpacing $conLE $DW
gg::conEndSpacing $conLE $DW

gg::conDistFunc $conTE -function TANH
gg::conBeginSpacing $conTE $DW
gg::conEndSpacing $conTE $DW

gg::domTFISolverRun [list $dom1 $dom2 $dom3]

# Block

gg::blkExtrusionBegin [list $dom1 $dom2 $dom3] -type STRUCTURED
  gg::blkExtrusionAtt TRANSLATE \
	-direction {0 0 1} \
	-equal_spacing \
	-distance $LZ
  gg::blkExtrusionStep [expr $NZ - 1]
set blocks [gg::blkExtrusionEnd]

# B.C.

puts "blocks =  $blocks"
set block1 [lindex $blocks 0]
set block2 [lindex $blocks 1]
set block3 [lindex $blocks 2]

gg::aswCreateBC "Hub" -solid 1 -id 6
set hubs(1) [gg::blkGetFace $block1 1]
set hubs(2) [gg::blkGetFace $block2 1]
set hubs(3) [gg::blkGetFace $block3 1]
gg::aswSetBC [list $hubs(1) $hubs(2) $hubs(3)] "Hub"

gg::aswCreateBC "Shroud" -solid 1 -id 7
set shrouds(1) [gg::blkGetFace $block1 2]
set shrouds(2) [gg::blkGetFace $block2 2]
set shrouds(3) [gg::blkGetFace $block3 2]
gg::aswSetBC [list $shrouds(1) $shrouds(2) $shrouds(3)] "Shroud"

gg::aswSetBC [gg::blkGetFace $block1 3] "Inflow"

gg::aswSetBC [gg::blkGetFace $block3 4] "Outflow"

set blade(1) [gg::blkGetFace $block2 5]
set blade(2) [gg::blkGetFace $block2 6]
gg::aswSetBC [list $blade(1) $blade(2)] "Wall"

gg::aswCreateBC "Periodic1" -solid 0 -id 8
set per1(1) [gg::blkGetFace $block1 5]
set per1(2) [gg::blkGetFace $block3 5]
gg::aswSetBC [list $per1(1) $per1(2)] "Periodic1"

gg::aswCreateBC "Periodic2" -solid 0 -id 9
set per1(1) [gg::blkGetFace $block1 6]
set per1(2) [gg::blkGetFace $block3 6]
gg::aswSetBC [list $per1(1) $per1(2)] "Periodic2"

gg::aswExport "DCA.cgns"

