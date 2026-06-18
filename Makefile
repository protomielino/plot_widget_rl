CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra

CFLAGS  += $(shell pkg-config --cflags raylib)
LDLIBS  += $(shell pkg-config --libs raylib) -lm

TARGET  = plot_widget_rl
SRCS    = main.c plot_widget.c izhikevich.c
OBJS    = $(SRCS:.c=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) $(OBJS)
