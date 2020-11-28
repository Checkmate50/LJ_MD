F=simple
SOURCE=$(F).xyz
DIR=results/$(F)
SEARCH=results/$(F)*

CC=gcc

build: 
	$(CC) lj_md.c -o lj_md.out

run:
	rm -rf results/*
	./lj_md.out $(SOURCE) $(DIR)

vmd:
	find $(SEARCH) -type f -print0 -exec vmd {} +
	
clean:
	rm -f *.out *.stackdump
	rm -r results/*