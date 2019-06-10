CC = g++

CFLAGS = -g -Wall -w

TARGET = raytracer

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp

clean:
	rm $(TARGET)
