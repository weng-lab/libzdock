#
#	On BSD machines, RANLIB should be 'ranlib'
#
#	On System V machines, RANLIB should be ':'
#
SHELL		= /bin/sh
RANLIB		= :

CXX		= gcc -felide-constructors
#CXX		= CC

.SUFFIXES:	.cc

.cc.o:
	$(CXX) $(CCFLAGS) -c $< -o $@

OPT		= -O
DEBUG		=
CCFLAGS		= $(OPT) $(DEBUG)
LIBRARY		= pdb++

LIBARCH		= lib$(LIBRARY).a
OBJS		= pdb_read.o pdb_sprntf.o pdb_sscanf.o pdb_chars.o \
		pdb_type.o pdb++.o pdbinput.o
SRCS		= pdb_read.cc pdb_sprntf.cc pdb_sscanf.cc pdb_chars.cc \
		pdb_type.cc pdb++.cc pdbinput.cc

all:		$(LIBARCH)

install:	$(LIBARCH)
		install -F /usr/local/lib $(LIBARCH)

$(LIBARCH):     $(OBJS)
		@echo "Loading $(LIBARCH) ... "
		@ar cru $(LIBARCH) $(OBJS)
		@$(RANLIB) $(LIBARCH)
		@echo "done"

clean:;		@rm -f $(OBJS)

spotless:;	@rm -f $(OBJS) $(LIBARCH)

tags:           $(HDRS) $(SRCS); @ctags -w $(HDRS) $(SRCS)
