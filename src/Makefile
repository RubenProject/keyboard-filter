CC = g++

CXXFLAGS = -c -Wall -std=c++11
LDFLAGS = ./libwavelib.a
OBJS = FilterGenerator.o AudioFile.o main.o

filter: $(OBJS)
	$(CC) $(OBJS) -o filter $(LDFLAGS)

FilterGenerator.o: FilterGenerator.cpp FilterGenerator.h AudioFile.h
	$(CC) $(CXXFLAGS) FilterGenerator.cpp

AudioFile.o: AudioFile.cpp AudioFile.h
	$(CC) $(CXXFLAGS) AudioFile.cpp

main.o: main.cpp FilterGenerator.h AudioFile.h
	$(CC) $(CXXFLAGS) main.cpp
