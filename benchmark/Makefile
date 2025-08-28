# Prerequisite Unix tools: ls, wc, echo, cut, bc
LSD ?= ../Build/muLSD
COMPARE ?= ../Build/compare_lines

# LSD expects .pgm as input, replace if needed
EXT ?= png

# Database for test: all images with extension EXT0
EXT0 = jpg
DB = $(wildcard *.$(EXT0))

all : eval

# Use ImageMagick to convert to required format, if needed
ifneq ($(EXT),$(EXT0))
%.$(EXT) : %.$(EXT0)
	magick $< $@
endif

.PRECIOUS : %_seg.txt
%_seg.txt : %.$(EXT)
	$(LSD) $< $@

%_eval.txt : %_seg.txt
	$(COMPARE) $< $(subst seg,gt,$<) > $@

eval : $(subst .$(EXT0),_eval.txt,$(DB))
	echo Based on $(shell ls *_eval.txt |wc -l) images:

	echo -n \(0 > dummy.txt
	for i in $^; do \
		echo -n + >> dummy.txt; \
		cut -d ' ' -f 2 $$i | tr -d '\n' >> dummy.txt; \
	done
	echo \) / $(shell ls *_eval.txt |wc -l) >> dummy.txt
	echo -n P=
	bc -l < dummy.txt

	echo -n \(0 > dummy.txt
	for i in $^; do \
		echo -n + >> dummy.txt; \
		cut -d ' ' -f 3 $$i | tr -d '\n' >> dummy.txt; \
	done
	echo \) / $(shell ls *_eval.txt |wc -l) >> dummy.txt
	echo -n R=
	bc -l < dummy.txt

	echo -n \(0 > dummy.txt
	for i in $^; do \
		echo -n + >> dummy.txt; \
		cut -d ' ' -f 5 $$i | tr -d '\n' >> dummy.txt; \
	done
	echo \) / $(shell ls *_eval.txt |wc -l) >> dummy.txt
	echo -n IOU=
	bc -l < dummy.txt

	rm dummy.txt

.PHONY : clean
clean :
	rm -f dummy.txt *_eval.txt *_seg.txt
