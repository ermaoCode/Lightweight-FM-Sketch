objects = main.o LFM.o murmur3.o hyperloglog.o

ELDA: $(objects)
	g++ -o ELDA $(objects)

main.o: main.cc hyperloglog.h LFM.h
	g++ -c $<
LFM.o: LFM.cc LFM.h
	g++ -c $<
hyperloglog.o: hyperloglog.cc hyperloglog.h murmur3.h
	g++ -c $<
murmur3.o: murmur3.c murmur3.h
	g++ -c $< 

.PHONY: clean
clean:
	rm ELDA $(objects)

