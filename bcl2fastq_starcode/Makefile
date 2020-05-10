SHELL := /bin/bash
.DELETE_ON_ERROR:

#===============================================================================
# VARIABLES

# sequencing run ids
# RUNS := $(shell find -L ../data/seq-runs -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)
RUNS := 200402_NB552046_0038_AHLHM2BGXF

#===============================================================================
# RECIPIES

all: starcode kb
starcode: $(addprefix pipeline/, $(addsuffix /starcode.csv, $(RUNS)))
kb: $(addprefix pipeline/, $(addsuffix /kb.tsv, $(RUNS)))

# cleanup
clean:
	rm -rf pipeline/* ../data/expected-amplicons.idx

.PRECIOUS: $(addprefix pipeline/, %/conditions.csv %/starcode.csv %/output.bus %/kb.tsv)
.SECONDARY:

#===============================================================================
# GET DATA!

../data/seq-runs/200402_NB552046_0038_AHLHM2BGXF:
	@mkdir -p $(@D)
	@wget --save-cookies cookies.txt \
	    'https://docs.google.com/uc?export=download&id=1GVhErh_LZm2R9m7LOCJouI7pN2bHhE4s'\
	    -O - \
         | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p' \
	 > confirm.txt \
	 && wget --load-cookies cookies.txt \
	 -O $@.tar.gz \
	 'https://docs.google.com/uc?export=download&id=1GVhErh_LZm2R9m7LOCJouI7pN2bHhE4s&confirm='$$(<confirm.txt) \
	 && tar -I pigz -xf $@.tar.gz -C $@ \
	 && rm cookies.txt confirm.txt

#===============================================================================
# STARCODE PIPELINE
# 1) demultiplex all i5/i7 pairs
# 2) run starcode on all resulting fastq's to collapse sequences lev=2 apart

# 0) Directory Prep:
# ------------------
# copy barcode map from data/ folder.
pipeline/%/bc-map.csv:
	@echo "Grabbing $@"
	@mkdir -p $(dir $@)
	@cp ../data/barcode-map.csv $@

# grab relevant section of samplesheet (make sure to catch windows return)
pipeline/%/conditions.csv: ../data/plate-maps/%/SampleSheet.csv pipeline/%/bc-map.csv
	@echo "Parsing $<"
	@python src/strip-windows.py $< \
	    | awk '/index2/,0' \
	    | Rscript src/bc-lengths.R $(word 2, $^) $@ \
	    2> $(@:.csv=.err)

# 1) Demultiplexing:
# ------------------
# starcode compatible demux
pipeline/%/demux-starcode/demux: ../data/seq-runs/% ../data/plate-maps/%/SampleSheet.csv
	@echo "Demultiplexing all i5/i7 combos in $< for starcode"
	@mkdir -p $(@D)
	@bcl2fastq \
	    --runfolder-dir $< \
	    --output-dir $(@D) \
	    --sample-sheet $(lastword $^) \
	    --processing-threads 8 \
	    --no-lane-splitting \
	    2> $(@D)/bcl2fastq.out \
	    && touch $@

# 2) Starcode:
# ------------
# Run starcode on each fastq, using gnu parallel to parse the conditions for barcode length
pipeline/%/starcode.csv: pipeline/%/demux-starcode/demux pipeline/%/conditions.csv
	@echo "Counting BCs for all fastq's in $(<D)"
	@parallel --header : --colsep "," \
		zcat $(<D)/"{Sample_ID}"_S*_R1_001.fastq.gz \
		\| awk -v bc_len="{bc_len}" -f src/count-bcs.awk \
		\| starcode -d2 -t1 --sphere --print-clusters 2> /dev/null \
		\| python src/tidy-star.py \
		\| awk -v name="{Sample_ID}" \''{print name, $$1, $$2, $$3}'\' OFS="," \
	:::: $(word 2, $^) \
	2> $(@:.csv=.err) \
	| sort -t ',' -k 1,1 -k 2,2 -k 3,3nr -k 4,4 \
	> $(@:.csv=.tmp) \
	&& echo "Sample_ID,Centroid,Count,barcode" \
	| cat - $(@:.csv=.tmp) > $@ \
	&& rm $(@:.csv=.tmp)

