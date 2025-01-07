ROOT_FILE 		= thesis
BACKUP_FILES 	= *~
CLEAN_FILES 	= $(BACKUP_FILES) *-frn.tex *.fls *.acn *.acr *.alg *.aux *.bcf *.bbl *.blg *.dvi *.fdb_latexmk *.glg *.glo *.gls *.idx *.ilg *.ind *.ist *.lof *.log *.lot *.lol *.maf *.mtc *.mtc0 *.nav *.nlo *.out *.pdfsync *.ps *.snm *.synctex.gz *.toc *.vrb *.xdy *.tdo *.run.xml *-blx.bib Chapters/*.aux
DIST_FILES 		= $(ROOT_FILE).pdf

.PHONY: clean cleanall tex frontispiece biber prepare simplethesis thesis

clean:
	echo $(CLEAN_FILES) | xargs $(RM); \
	find . -name '$(BACKUP_FILES)' | xargs $(RM);

cleanall: clean
	$(RM) $(DIST_FILES)

tex:
	latexmk -silent -f -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make $(ROOT_FILE).tex

biber:
	biber $(ROOT_FILE)

prepare: tex
	$(MAKE) biber

simplethesis: prepare
	$(MAKE) tex

thesis: prepare
	$(MAKE) frontispiece
	$(MAKE) tex
