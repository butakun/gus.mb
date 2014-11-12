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

set NARC 100
set NFLAT 50
set NOUT 100
set NFREE [expr $NARC + 2 * $NFLAT - 2]

set DW 1e-5

# Database

gg::dbCurveBegin -type CIRCULAR_ARC
  gg::dbCurveAddPt {0.0 -1.0 0.0}
  gg::dbCurveAddPt {0.0 1.0 0.0}
  gg::dbCurveAddPt {-1.0 0.0 0.0}
set dbArc [gg::dbCurveEnd]

gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 1 0 $dbArc]
  gg::dbCurveAddPt {2.0 1.0 0.0}
set dbUp [gg::dbCurveEnd]

gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 0 0 $dbArc]
  gg::dbCurveAddPt {2.0 -1.0 0.0}
set dbLo [gg::dbCurveEnd]

# Connectors

set conArc [gg::conOnDBEnt $dbArc]
gg::conDim $conArc $NARC

set conUp [gg::conOnDBEnt $dbUp]
gg::conDim $conUp $NFLAT

set conLo [gg::conOnDBEnt $dbLo]
gg::conDim $conLo $NFLAT

gg::conBegin
  gg::segBegin -type AKIMA
    gg::segAddControlPt {2.0 -10.0 0.0}
    gg::segAddControlPt {-2.0 0.0 0.0}
    gg::segAddControlPt {2.0 10.0 0.0}
  gg::segEnd
set conFree [gg::conEnd]
gg::conDim $conFree $NFREE

gg::conBegin
  gg::segBegin -type 3D_LINE
    gg::segAddControlPt -db [list 1 0 $dbUp]
    gg::segAddControlPt [gg::conGetPt $conFree -arc 1]
  gg::segEnd
set conOutUp [gg::conEnd]
gg::conDim $conOutUp $NOUT

gg::conBegin
  gg::segBegin -type 3D_LINE
    gg::segAddControlPt -db [list 1 0 $dbLo]
    gg::segAddControlPt [gg::conGetPt $conFree -arc 0]
  gg::segEnd
set conOutLo [gg::conEnd]
gg::conDim $conOutLo $NOUT

gg::conBeginSpacing $conOutUp $DW
gg::conBeginSpacing $conOutLo $DW


# Domain

gg::domBegin -type STRUCTURED
  gg::edgeBegin
    gg::edgeAddCon $conLo
    gg::edgeAddCon $conArc
    gg::edgeAddCon $conUp
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $conOutUp
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $conFree
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $conOutLo
  gg::edgeEnd
set dom [gg::domEnd]


# Block

gg::blkExtrusionBegin $dom -type STRUCTURED
  gg::blkExtrusionAtt TRANSLATE \
        -direction {0 0 1} \
        -equal_spacing \
        -distance 1.0
  gg::blkExtrusionStep 2
set blocks [gg::blkExtrusionEnd]


# B.C.

set block [lindex $blocks 0]
puts $blocks
puts $block

set walls [lindex [gg::blkGetFace $block 5] 0]
gg::aswSetBC $walls "Wall"
gg::aswSetBC [gg::blkGetFace $block 6] "Inflow"
gg::aswSetBC [gg::blkGetFace $block 3] "Outflow"
gg::aswSetBC [gg::blkGetFace $block 4] "Outflow"
gg::aswSetBC [gg::blkGetFace $block 1] "Symmetry"
gg::aswSetBC [gg::blkGetFace $block 2] "Symmetry"

gg::aswExport "Blunt.cgns"

