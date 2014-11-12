proc ReadAirfoil {filename} {
	set f [open $filename r]

	gets $f line
	gets $f nSections

	set sections {}
	for {set i 0} {$i < $nSections} {incr i} {
		set section {}
		gets $f line

		gets $f line
		gets $f nPts
		gets $f line
		set pts {}
		for {set j 0} {$j < $nPts} {incr j} {
			gets $f line
			set xyz [split [join $line]]
			lappend pts $xyz
		}
		lappend section $pts

		gets $f line
		gets $f nPts
		gets $f line
		set pts {}
		for {set j 0} {$j < $nPts} {incr j} {
			gets $f line
			set xyz [split [join $line]]
			lappend pts $xyz
		}
		lappend section $pts

		lappend sections $section
	}

	return $sections
}

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

# Geometric parameters
set X1 -40.0
set XLE 0.0
set YLE 0.0
set XTE 36.985
set YTE -52.3
set Ymid [expr 0.5 * $YTE]
set Pitch 60.0
set X2 [expr $XTE + 40.0]
set LZ 1.0
set NZ 2

# Mesh parameters
set DW 1e-3
set DSLE [expr 3 * $DW] 
set DSTE [expr 2 * $DW]
set DSUP 5.0
set DSWAKE 5.0
set NInflow 100
set NUP 50
set NWAKE 50
set NBLADE 100

# Database
set sections [ReadAirfoil "LS89.xyz"]

set section [lindex $sections 0]
set pts1 [lindex $section 0]
set pts2 [lindex $section 1]

gg::dbCurveBegin -type 3D_LINE
	foreach p $pts1 {
		gg::dbCurveAddPt $p
	}
set curveSS [gg::dbCurveEnd]

gg::dbCurveBegin -type 3D_LINE
	foreach p $pts2 {
		gg::dbCurveAddPt $p
	}
set curvePS [gg::dbCurveEnd]
gg::dbTransformBegin [list $curvePS]
  gg::xformTranslate [list 0 $Pitch 0]
gg::dbTransformEnd

set conPS [gg::conOnDBEnt $curvePS]
gg::conDim $conPS $NBLADE

set conSS [gg::conOnDBEnt $curveSS]
gg::conDim $conSS $NBLADE

gg::conBegin
  gg::segBegin -type 3D_LINE
    gg::segAddControlPt [list $X1 $YLE 0]
    gg::segAddControlPt -db [list 0 0 $curveSS]
  gg::segEnd
set conUP [gg::conEnd]
gg::conDim $conUP $NUP

gg::conBegin
  gg::segBegin -type AKIMA
    gg::segAddControlPt -db [list 1 0 $curveSS]
    gg::segAddControlPt [list [expr $XTE + 15] [expr $YTE - 15] 0]
    gg::segAddControlPt [list $X2 [expr $YTE - 25] 0]
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
    gg::edgeAddCon $conSS
    gg::edgeAddCon $conWAKE
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $conOutflow
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $conWAKETrans
    gg::edgeAddCon $conPS
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
gg::conEndSpacing $conUP $DW
gg::conBeginSpacing $conWAKE $DW

gg::conDistFunc $conLE -function TANH
gg::conBeginSpacing $conLE $DW
gg::conEndSpacing $conLE $DW

gg::conDistFunc $conTE -function TANH
gg::conBeginSpacing $conTE $DW
gg::conEndSpacing $conTE $DW

gg::domTFISolverRun [list $dom1 $dom2 $dom3]

gg::domEllSolverBegin [list $dom1 $dom2]
  gg::domEllSolverStep -iterations 1000
gg::domEllSolverEnd

# Block

gg::blkExtrusionBegin [list $dom1 $dom2 $dom3] -type STRUCTURED
  gg::blkExtrusionAtt TRANSLATE \
	-direction {0 0 1} \
	-equal_spacing \
	-distance $LZ
  gg::blkExtrusionStep $NZ
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

gg::aswExport "LS89.cgns"

