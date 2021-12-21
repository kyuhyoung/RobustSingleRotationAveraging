############################################################
####	Makefile template for opencv use
############################################################
RM 			= @rm -rfv
CC = g++
CFLAGS = -g -Wall
SRCS = $(wildcard ./*.cpp)
#PROG = HelloWorld
PROG      = $(notdir $(shell pwd))	## current foldername is target name

OPENCV = `pkg-config opencv --cflags --libs`
LIBS = $(OPENCV) -lstdc++fs

$(PROG):$(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

#	clean
clean:
#	$(RM) $(BUILD_PATH)
	$(RM) $(PROG)
	@echo "Done."


