CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra

CFLAGS  += $(shell pkg-config --cflags raylib)
LDLIBS  += $(shell pkg-config --libs raylib) -lm

TARGET  = plot_widget_rl
SRC     = main.c

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(TARGET)
