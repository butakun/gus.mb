# $Id: Makefile 59 2010-09-01 08:53:01Z kato $

%.pdf:	%.tex
	pdflatex $<
	pdflatex $<

all::	Scheme.pdf

Scheme.pdf:	Scheme.tex GoverningEquations.tex TurbulenceModels.tex
