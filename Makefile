all: bin/processo-raytracer bin/open-mp-raytracer

bin:
	mkdir bin

bin/processo-raytracer: processo/processo-raytracer.c bin
	gcc processo/processo-raytracer.c -lm -o bin/processo-raytracer

bin/open-mp-raytracer: open-mp/open-mp-raytracer.c bin
	gcc open-mp/open-mp-raytracer.c -lm -o bin/open-mp-raytracer
clean:

mrproper: clean
	rm -rf bin/*