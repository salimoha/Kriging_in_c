# created by Shahrouz Alimohammadi
# Note: it is easy to modify this Makefile
# only change the paths to the directory where you have boost and suitsparse and LAPACK library 
# replace only this part "/home/shahrouz/software/" in the file
# 09-2013
# shahoruz.alm@gmail.com

# How to use it-> type in the command line:
#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = g++

# define any compile-time flags
#CFLAGS = -Wall -g

# define any directories containing header files other than /usr/include
# /home/shahrouz/
INCLUDES = -I/home/shahrouz/libc/SuiteSparse/SPQR/Include/ -I/home/shahrouz/libc/SuiteSparse/CHOLMOD/Include/ -I/home/shahrouz/libc/SuiteSparse/SuiteSparse_config/ -I/home/shahrouz/libc/boost_1_55_0/

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -p -g -L/home/shahrouz/libc/SuiteSparse/SPQR/Lib/ -L/home/shahrouz/libc/SuiteSparse/CHOLMOD/Lib/ -L/home/shahrouz/libc/SuiteSparse/AMD/Lib/ -L/home/shahrouz/libc/SuiteSparse/COLAMD/Lib/ -L/home/shahrouz/libc/SuiteSparse/SuiteSparse_config/

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lspqr -lcholmod -lamd -lcolamd -lsuitesparseconfig -lblas -llapack -lrt

# define the C source files
SRCS = main.cpp

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file 
MAIN = testApp

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  Simple compiler named mycc has been compiled

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

run:
	./$(MAIN)

clean:
	$(RM) fun.txt *.o *~ $(MAIN) 

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
