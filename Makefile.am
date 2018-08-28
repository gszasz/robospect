AUTOMAKE_OPTIONS = foreign
SUBDIRS = src

archive : VERSION := $(shell cvs status src/main.c | grep Working | awk '{print($$3)}')
archive : 
#	cvs2cl -r -t
	tar -zchf /tmp/robospect-${VERSION}.tar.gz --transform 's%^%./robospect-${VERSION}/%' --transform 's%manual/%%' INSTALL README LICENSE doc/QUICKSTART doc/manual/robospect_user_guide.pdf  src/*.c src/fileio/*.c src/models/*.c src/*.h src/math/*.c utilities/lines2moog.pl utilities/robo_decode_flags.pl utilities/robo_simulate_spectra.pl src/Makefile.in src/Makefile.am src/config.h.in Makefile.in Makefile.am configure depcomp install-sh missing compile configure.ac 
	tar -zchf /tmp/robospect-${VERSION}-compat.tar.gz --transform 's%.compat%%' --transform 's%^%./robospect-${VERSION}-compat/%' --transform 's%manual/%%'  INSTALL README LICENSE doc/QUICKSTART doc/manual/robospect_user_guide.pdf  src/*.c src/fileio/*.c src/models/*.c src/*.h src/math/*.c utilities/lines2moog.pl utilities/robo_decode_flags.pl utilities/robo_simulate_spectra.pl src/Makefile.in src/Makefile.am.compat src/config.h.in Makefile.in Makefile.am configure depcomp install-sh missing compile configure.ac 
