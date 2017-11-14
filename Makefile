all: bin/processo-raytracer bin/open-mp-raytracer bin/threads-raytracer bin/bag-threads-raytracer

bin:
	mkdir bin

bin/processo-raytracer: processo/processo-raytracer.c bin
	gcc processo/processo-raytracer.c -lm -o bin/processo-raytracer

bin/open-mp-raytracer: open-mp/open-mp-raytracer.c bin
	gcc open-mp/open-mp-raytracer.c -lm -o bin/open-mp-raytracer

bin/threads-raytracer: threads/threads-raytracer.c bin
	gcc threads/threads-raytracer.c -lm -pthread -o bin/threads-raytracer

bin/bag-threads-raytracer: bag-of-tasks/bag-of-tasks-threads-raytracer.c bin
	gcc bag-of-tasks/bag-of-tasks-threads-raytracer.c -lm -pthread -o bin/bag-threads-raytracer

clean:

mrproper: clean
	rm -rf bin
	rm -rf *.bmp