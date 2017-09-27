all:
	gcc -Wall src/main.c src/util.c src/matrix.c -g -lm -o invmat
doc:
	doxygen Doxyfile
clean:
	rm invmat