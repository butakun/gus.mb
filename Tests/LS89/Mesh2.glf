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

package require PWI_Glyph 2.0

pw::Application reset
pw::Application setCAESolver CGNS

# Geometric parameters
set X1 -40.0
set XLE 0.0
set YLE 0.0
set XTE 36.985
set YTE -52.3
set Ymid [expr 0.5 * $YTE]
set Pitch 60.0
set X2 [expr $XTE + 40.0]
set LZ 0.01
set NZ 2

set UThroat 0.7

# Mesh parameters
set DW 1e-3
set DSLE 0.1
set DSTE 0.05
set DSUP 5.0
set DSWAKE 5.0
set NInflow 50
set NOutflow 30
set NUP 30
set NWAKE $NInflow
set NPS 50

set NPS_ $NPS
set NInflow_ $NInflow
set NOutflow_ $NOutflow
set NUP_ $NUP

# Database
set sections [ReadAirfoil "LS89.xyz"]

set section [lindex $sections 0]
set pts1 [lindex $section 0]
set pts2 [lindex $section 1]

proc CreateCurve {pts} {
	set creator [pw::Application begin Create]
	set seg [pw::SegmentSpline create]
	foreach p $pts {
		$seg addPoint $p
	}
	set curve [pw::Curve create]
	$curve addSegment $seg
	$creator end
	return $curve
}

set curveSS [CreateCurve $pts1]
set curvePS [CreateCurve $pts2]

pw::Entity transform [pwu::Transform translation [list 0 $Pitch 0]] [list $curvePS]

set pLE [$curveSS getPosition 0.0]
set pTE [$curveSS getPosition 1.0]
set pTEPS [$curvePS getPosition 1.0]
set XLE [lindex $pLE 0]
set XTE [lindex $pTE 0]
set YTE [lindex $pTE 1]
set YTEPS [lindex $pTEPS 1]
set CAX [expr $XTE - $XLE]

set creator [pw::Application begin Create]

set seg [pw::SegmentTrim create]
$seg addPoint [list 0 0 $curvePS]
$seg addPoint [list 1 0 $curvePS]
set con [pw::Connector create]
$con setName "PS"
$con addSegment $seg
$con setDimension $NPS_
set conPS $con

set seg [pw::SegmentTrim create]
$seg addPoint [list 0 0 $curveSS]
$seg addPoint [list $UThroat 0 $curveSS]
set con [pw::Connector create]
$con setName "SSFront"
$con addSegment $seg
$con setDimension $NPS_
set conSSFront $con

set seg [pw::SegmentSpline create]
$seg addPoint [$curveSS getPosition 0.0]
$seg addPoint [list [expr $XLE - 0.2 * $Pitch] [expr $YLE + 0.5 * $Pitch] 0]
$seg addPoint [$curvePS getPosition 0.0]
$seg setSlope Akima
set con [pw::Connector create]
$con setName "passageLE"
$con addSegment $seg
$con setDimension $NInflow_
set conPassageLE $con

set seg [pw::SegmentSpline create]
$seg addPoint [$curveSS getPosition $UThroat]
$seg addPoint [$curvePS getPosition 1.0]
set con [pw::Connector create]
$con setName "passageThroat"
$con addSegment $seg
$con setDimension $NInflow_
set conPassageThroat $con

set seg [pw::SegmentTrim create]
$seg addPoint [list $UThroat 0 $curveSS]
$seg addPoint [list 1 0 $curveSS]
set con [pw::Connector create]
$con setName "SSAft"
$con addSegment $seg
$con setDimension [expr $NInflow_ + $NOutflow_ - 1]
set conSSAft $con

set seg [pw::SegmentSpline create]
$seg addPoint [$curveSS getPosition 1.0]
$seg addPoint [list [expr $XTE + 0.5 * $CAX] [expr $YTE - 0.6 * $Pitch] 0]
$seg addPoint [list [expr $XTE + 1.0 * $CAX] [expr $YTE - 0.95 * $Pitch] 0]
$seg setSlope Akima
set con [pw::Connector create]
$con setName "Wake"
$con addSegment $seg
$con setDimension $NInflow_
set conWake $con

set conWakePS [$conWake createPeriodic -translate [list 0 $Pitch 0]]

set seg [pw::SegmentSpline create]
$seg addPoint [$conWake getPosition -parameter 1.0]
$seg addPoint [$conWakePS getPosition -parameter 1.0]
set con [pw::Connector create]
$con setName "Outflow"
$con addSegment $seg
$con setDimension $NOutflow_
set conOutflow $con

set seg [pw::SegmentSpline create]
$seg addPoint [list [expr $XLE - $CAX] $YLE 0]
$seg addPoint $pLE
set con [pw::Connector create]
$con setName "CutUp]"
$con addSegment $seg
$con setDimension $NUP_
set conCutUp $con

set conCutUpPS [$conCutUp createPeriodic -translate [list 0 $Pitch 0]]

set seg [pw::SegmentSpline create]
$seg addPoint [$conCutUp getPosition -parameter 0.0]
$seg addPoint [$conCutUpPS getPosition -parameter 0.0]
set con [pw::Connector create]
$con setName "Inflow"
$con addSegment $seg
$con setDimension $NInflow_
set conInflow $con

[[$conPS getDistribution 1] getBeginSpacing] setValue $DSLE
[[$conPS getDistribution 1] getEndSpacing] setValue $DSTE
[[$conSSFront getDistribution 1] getBeginSpacing] setValue $DSLE
[[$conSSFront getDistribution 1] getEndSpacing] setValue [expr 2.0 * $DSTE]
[[$conSSAft getDistribution 1] getBeginSpacing] setValue [expr 2.0 * $DSTE]
[[$conSSAft getDistribution 1] getEndSpacing] setValue $DSTE
[[$conWake getDistribution 1] getBeginSpacing] setValue $DW
[[$conWake getDistribution 1] getEndSpacing] setValue [expr 5.0 * $DSTE]

[[$conCutUp getDistribution 1] getEndSpacing] setValue $DSLE
[[$conPassageLE getDistribution 1] getBeginSpacing] setValue $DW
[[$conPassageLE getDistribution 1] getEndSpacing] setValue $DW
[[$conPassageThroat getDistribution 1] getBeginSpacing] setValue $DW
[[$conPassageThroat getDistribution 1] getEndSpacing] setValue $DW

set dom1 [pw::DomainStructured createFromConnectors [list $conCutUp $conPassageLE $conCutUpPS $conInflow]]
set dom2 [pw::DomainStructured createFromConnectors [list $conSSFront $conPassageThroat $conPS $conPassageLE]]

set edge1 [pw::Edge createFromConnectors [list $conSSAft]]
set edge2 [pw::Edge createFromConnectors [list $conWake]]
set edge3 [pw::Edge createFromConnectors [list $conOutflow $conWakePS]]
set edge4 [pw::Edge createFromConnectors [list $conPassageThroat]]
set dom3 [pw::DomainStructured create]
$dom3 addEdge $edge1
$dom3 addEdge $edge2
$dom3 addEdge $edge3
$dom3 addEdge $edge4

$creator end

# Elliptic Solve

set solver [pw::Application begin EllipticSolver [list $dom2]]
$dom2 setEllipticSolverAttribute -edge 1 EdgeAngleCalculation Orthogonal
$dom2 setEllipticSolverAttribute -edge 2 EdgeAngleCalculation Adjacent
$dom2 setEllipticSolverAttribute -edge 3 EdgeAngleCalculation Orthogonal
$dom2 setEllipticSolverAttribute -edge 4 EdgeAngleCalculation Current
$solver run Initialize
$solver run 100
$solver end

set solver [pw::Application begin EllipticSolver [list $dom3]]
$dom3 setEllipticSolverAttribute -edge 1 EdgeAngleCalculation Orthogonal
$dom3 setEllipticSolverAttribute -edge 2 EdgeAngleCalculation Orthogonal
$dom3 setEllipticSolverAttribute -edge 3 EdgeAngleCalculation Orthogonal
$dom3 setEllipticSolverAttribute -edge 4 EdgeAngleCalculation Adjacent
$solver run Initialize
$solver run 100
$solver end

set solver [pw::Application begin EllipticSolver [list $dom1]]
$dom1 setEllipticSolverAttribute -edge 1 EdgeAngleCalculation Current
$dom1 setEllipticSolverAttribute -edge 2 EdgeAngleCalculation Orthogonal
$dom1 setEllipticSolverAttribute -edge 3 EdgeAngleCalculation Orthogonal
$dom1 setEllipticSolverAttribute -edge 4 EdgeAngleCalculation Orthogonal
$solver run Initialize
$solver run 100
$solver end

# Extrusion

set creator [pw::Application begin Create]

set face1 [pw::FaceStructured createFromDomains $dom1]
set face2 [pw::FaceStructured createFromDomains $dom2]
set face3 [pw::FaceStructured createFromDomains $dom3]
set block1 [pw::BlockStructured create]
$block1 addFace $face1
set block2 [pw::BlockStructured create]
$block2 addFace $face2
set block3 [pw::BlockStructured create]
$block3 addFace $face3

$creator end

set solver [pw::Application begin ExtrusionSolver [list $block1 $block2 $block3]]
$block1 setExtrusionSolverAttribute Mode Translate
$block2 setExtrusionSolverAttribute Mode Translate
$block3 setExtrusionSolverAttribute Mode Translate
$block1 setExtrusionSolverAttribute TranslateDirection {0 0 1}
$block2 setExtrusionSolverAttribute TranslateDirection {0 0 1}
$block3 setExtrusionSolverAttribute TranslateDirection {0 0 1}
$block1 setExtrusionSolverAttribute TranslateDistance $LZ
$block2 setExtrusionSolverAttribute TranslateDistance $LZ
$block3 setExtrusionSolverAttribute TranslateDistance $LZ
$solver run 2
$solver end

set dom6 [pw::GridEntity getByName "dom-6"]
$dom6 setName "Inflow"
set dom16 [pw::GridEntity getByName "dom-16"]
$dom16 setName "Outflow"
set dom8 [pw::GridEntity getByName "dom-8"]
$dom8 setName "Shroud1"
set dom13 [pw::GridEntity getByName "dom-13"]
$dom13 setName "Shroud2"
set dom19 [pw::GridEntity getByName "dom-19"]
$dom19 setName "Shroud3"
set dom9 [pw::GridEntity getByName "dom-9"]
$dom9 setName "BladePS"
set dom11 [pw::GridEntity getByName "dom-11"]
$dom11 setName "BladeSS1"
set dom14 [pw::GridEntity getByName "dom-14"]
$dom14 setName "BladeSS2"
set dom7 [pw::GridEntity getByName "dom-7"]
$dom7 setName "PeriodicSS1"
set dom15 [pw::GridEntity getByName "dom-15"]
$dom15 setName "PeriodicSS2"
set dom5 [pw::GridEntity getByName "dom-5"]
$dom5 setName "PeriodicPS1"
set dom17 [pw::GridEntity getByName "dom-17"]
$dom17 setName "PeriodicPS2"

$dom1 setName "Hub1"
$dom2 setName "Hub2"
$dom3 setName "Hub3"

# Boundary Conditions

set bcHub [pw::BoundaryCondition create]
$bcHub setName "Hub"
$bcHub setPhysicalType {Symmetry Plane}
$bcHub apply [list [list $block1 $dom1] [list $block2 $dom2] [list $block3 $dom3]]

# Export

pw::Application setCAESolverAttribute "Structured" 1
pw::Application export -precision Single [pw::Entity sort [list $block1 $block2 $block3]] {LS89-2.cgns}

pw::Application exit

