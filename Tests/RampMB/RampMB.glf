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

# Database

gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt [list 0 0 0]
  gg::dbCurveAddPt [list 1 0.1 0]
set _DB(1) [gg::dbCurveEnd]
gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 1 0 $_DB(1)]
  gg::dbCurveAddPt [list 1 1 0]
set _DB(2) [gg::dbCurveEnd]
gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 1 0 $_DB(2)]
  gg::dbCurveAddPt [list 0 1 0]
set _DB(3) [gg::dbCurveEnd]
gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 1 0 $_DB(3)]
  gg::dbCurveAddPt -db [list 0 0 $_DB(1)]
set _DB(4) [gg::dbCurveEnd]

gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt [list 0 0 1]
  gg::dbCurveAddPt [list 1 0.1 1]
set _DB(5) [gg::dbCurveEnd]
gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 1 0 $_DB(5)]
  gg::dbCurveAddPt [list 1 1 1]
set _DB(6) [gg::dbCurveEnd]
gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 1 0 $_DB(6)]
  gg::dbCurveAddPt [list 0 1 1]
set _DB(7) [gg::dbCurveEnd]
gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 1 0 $_DB(7)]
  gg::dbCurveAddPt -db [list 0 0 $_DB(5)]
set _DB(8) [gg::dbCurveEnd]

gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 0 0 $_DB(1)]
  gg::dbCurveAddPt -db [list 0 0 $_DB(5)]
set _DB(9) [gg::dbCurveEnd]

gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 0 0 $_DB(2)]
  gg::dbCurveAddPt -db [list 0 0 $_DB(6)]
set _DB(10) [gg::dbCurveEnd]

gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 0 0 $_DB(3)]
  gg::dbCurveAddPt -db [list 0 0 $_DB(7)]
set _DB(11) [gg::dbCurveEnd]

gg::dbCurveBegin -type 3D_LINE
  gg::dbCurveAddPt -db [list 0 0 $_DB(4)]
  gg::dbCurveAddPt -db [list 0 0 $_DB(8)]
set _DB(12) [gg::dbCurveEnd]

set _CN(1) [gg::conOnDBEnt $_DB(1)]
set _CN(2) [gg::conOnDBEnt $_DB(2)]
set _CN(3) [gg::conOnDBEnt $_DB(3)]
set _CN(4) [gg::conOnDBEnt $_DB(4)]
set _CN(5) [gg::conOnDBEnt $_DB(5)]
set _CN(6) [gg::conOnDBEnt $_DB(6)]
set _CN(7) [gg::conOnDBEnt $_DB(7)]
set _CN(8) [gg::conOnDBEnt $_DB(8)]
set _CN(9) [gg::conOnDBEnt $_DB(9)]
set _CN(10) [gg::conOnDBEnt $_DB(10)]
set _CN(11) [gg::conOnDBEnt $_DB(11)]
set _CN(12) [gg::conOnDBEnt $_DB(12)]

gg::conDim $_CN(1) 100
gg::conDim $_CN(2) 50
gg::conDim $_CN(3) 100
gg::conDim $_CN(4) 50
gg::conDim $_CN(5) 100
gg::conDim $_CN(6) 50
gg::conDim $_CN(7) 100
gg::conDim $_CN(8) 50
gg::conDim $_CN(9) 10
gg::conDim $_CN(10) 10
gg::conDim $_CN(11) 10
gg::conDim $_CN(12) 10

#gg::conBeginSpacing $_CN(1) -sub 1 .02
#gg::conEndSpacing $_CN(1) -sub 1 .02

# DOMAINS

gg::domBegin -type STRUCTURED
  gg::edgeBegin
    gg::edgeAddCon $_CN(1)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(2)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(3)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(4)
  gg::edgeEnd
set _DM(1) [gg::domEnd]

gg::domBegin -type STRUCTURED
  gg::edgeBegin
    gg::edgeAddCon $_CN(5)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(6)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(7)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(8)
  gg::edgeEnd
set _DM(2) [gg::domEnd]

gg::domBegin -type STRUCTURED
  gg::edgeBegin
    gg::edgeAddCon $_CN(9)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(8)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(12)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(4)
  gg::edgeEnd
set _DM(3) [gg::domEnd]

gg::domBegin -type STRUCTURED
  gg::edgeBegin
    gg::edgeAddCon $_CN(10)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(6)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(11)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(2)
  gg::edgeEnd
set _DM(4) [gg::domEnd]

gg::domBegin -type STRUCTURED
  gg::edgeBegin
    gg::edgeAddCon $_CN(1)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(10)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(5)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(9)
  gg::edgeEnd
set _DM(5) [gg::domEnd]

gg::domBegin -type STRUCTURED
  gg::edgeBegin
    gg::edgeAddCon $_CN(3)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(11)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(7)
  gg::edgeEnd
  gg::edgeBegin
    gg::edgeAddCon $_CN(12)
  gg::edgeEnd
set _DM(6) [gg::domEnd]

# Blocks

gg::blkBegin -type STRUCTURED
  gg::faceBegin
    gg::faceAddDom $_DM(1)
  gg::faceEnd
  gg::faceBegin
    gg::faceAddDom $_DM(3)
  gg::faceEnd
  gg::faceBegin
    gg::faceAddDom $_DM(2)
  gg::faceEnd
  gg::faceBegin
    gg::faceAddDom $_DM(4)
  gg::faceEnd
  gg::faceBegin
    gg::faceAddDom $_DM(5)
  gg::faceEnd
  gg::faceBegin
    gg::faceAddDom $_DM(6)
  gg::faceEnd
set _BL(1) [gg::blkEnd]

gg::aswSetBC [list \
    $_DM(3) \
    ] \
  "Inflow"
gg::aswSetBC [list \
    $_DM(4) \
    ] \
  "Outflow"
gg::aswSetBC [list $_DM(5)] "Wall"
gg::aswSetBC [list $_DM(6)] "Wall"
gg::aswSetBC [list $_DM(1) $_DM(2)] "Wall"

#gg::blkJoin $_BL(2) $_BL(1)

set _ggTemp_(1) [gg::blkSplit $_BL(1) -i 50]
set _BL(2) [lindex $_ggTemp_(1) 0]
unset _ggTemp_(1)

set _ggTemp_(1) [gg::blkSplit $_BL(1) -j 20]
set _BL(3) [lindex $_ggTemp_(1) 0]
unset _ggTemp_(1)

gg::aswExport "RampMB.cgns"

