CC = g++

CXXFLAGS = -c -Wall -std=c++11
LIBS = -lboost_iostreams -lboost_system -lboost_filesystem
OBJS = FilterGenerator.o kiss_fft.o wavelet2d.o AudioFile.o main.o

filter: $(OBJS)
	$(CC) $(OBJS) -o filter $(LIBS)

FilterGenerator.o: FilterGenerator.cpp FilterGenerator.h AudioFile.h wavelet2d.h
	$(CC) $(CXXFLAGS) FilterGenerator.cpp $(LIBS)

kiss_fft.o: kiss_fft.c _kiss_fft_guts.h
	$(CC) $(CXXFLAGS) kiss_fft.c

wavelet2d.o: wavelet2d.cpp wavelet2d.h kiss_fft.h
	$(CC) $(CXXFLAGS) wavelet2d.cpp

AudioFile.o: AudioFile.cpp AudioFile.h
	$(CC) $(CXXFLAGS) AudioFile.cpp

main.o: main.cpp FilterGenerator.h wavelet2d.h AudioFile.h
	$(CC) $(CXXFLAGS) main.cpp
