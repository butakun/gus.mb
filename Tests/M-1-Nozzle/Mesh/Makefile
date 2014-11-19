PYTHON = python
AUTOGRID5 = /opt/numeca/bin/igg -autogrid5

all::	M-1-Nozzle.cgns

M-1-Nozzle.cgns:	M-1-Nozzle.geomTurbo
	$(AUTOGRID5) -print -batch -script ./AutoGrid5Script.py

M-1-Nozzle-geomRurbo:
	$(PYTHON) ./MakeGeomTurbo.py > M-1-Nozzle.geomTurbo

