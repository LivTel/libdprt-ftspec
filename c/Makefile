# $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/Makefile,v 0.9 2002-05-20 10:40:04 cjm Exp $

include	../../Makefile.common
include	../Makefile.common

INCDIR 		= $(LIBDPRT_SRC_HOME)/include
DOCSDIR 	= $(LIBDPRT_DOC_HOME)/cdocs
CFLAGS 		= -g $(CCHECKFLAG) -I$(INCDIR) -I$(CFITSIOINCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR) -L$(LT_LIB_HOME)
LINTFLAGS 	= -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
DOCFLAGS 	= -static
SRCS 		= DpRtLibrary.c dprt.c
HEADERS		= $(SRCS:%.c=%.h)
OBJS		= $(SRCS:%.c=%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)

top: shared docs

shared: $(INCDIR)/DpRtLibrary.h $(LT_LIB_HOME)/libdprt.so

$(INCDIR)/DpRtLibrary.h: $(DPRT_BIN_HOME)/java/DpRtLibrary.class
	$(JAVAH) -jni -d $(INCDIR) -classpath $(DPRT_BIN_HOME)/java:$(CLASSPATH) DpRtLibrary
	$(TOUCH) $@

$(LT_LIB_HOME)/libdprt.so: $(OBJS)
	cc $(CCSHAREDFLAG) $(CFLAGS) $(OBJS) -o $@ $(TIMELIB) -lcfitsio $(SOCKETLIB)

static: $(LT_LIB_HOME)/libdprt.a

$(LT_LIB_HOME)/libdprt.a: $(OBJS)
	ar rcv $@ $?

docs: $(DOCS)

$(DOCS) : $(SRCS)
	-$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

# dont checkout DpRtLibrary.h - it is a machine built header
checkout:
	$(CO) $(CO_OPTIONS) $(SRCS)
	cd $(INCDIR); $(CO) $(CO_OPTIONS) dprt.h;

# dont checkin DpRtLibrary.h - it is a machine built header
checkin:
	-$(CI) $(CI_OPTIONS) $(SRCS)
	-(cd $(INCDIR); $(CI) $(CI_OPTIONS) dprt.h;)

staticdepend:
	makedepend -- $(CFLAGS) -- $(SRCS)

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

lint:
	$(LINT)	$(LINTFLAGS) $(SRCS)

clean:
	-$(RM) $(RM_OPTIONS) $(OBJS) $(LT_LIB_HOME)/libdprt.so $(LT_LIB_HOME)/libdprt.a 
	-$(RM) $(RM_OPTIONS) $(INCDIR)/DpRtLibrary.h $(TIDY_OPTIONS)

tidy:
	-$(RM)  $(RM_OPTIONS) $(TIDY_OPTIONS)
	-cd $(INCDIR); $(RM) $(RM_OPTIONS) $(TIDY_OPTIONS);

backup: tidy
	-$(RM) $(RM_OPTIONS) $(OBJS)

#
# $Log: not supported by cvs2svn $
# Revision 0.8  2001/05/15 18:38:38  cjm
# Split dprt into libdprt.
#
# Revision 0.7  1999/12/10 12:21:44  cjm
# Added backup target.
#
# Revision 0.6  1999/09/08 13:56:34  cjm
# Changed checkin target to use CI_OPTIONS.
#
# Revision 0.5  1999/08/31 12:50:15  cjm
# Changed LIBDIR to LT_LIB_HOME.
#
# Revision 0.4  1999/08/24 16:05:29  cjm
# Changed for libccd1 directory structure.
# Improved docs rule.
#
# Revision 0.3  1999/07/12 16:05:00  cjm
# Changes to ignore CDoc return value.
#
# Revision 0.2  1999/06/25 09:50:59  dev
# changed tidy and checkin rules so that errors are ignored
# (e.g. files already deleted/files already checked in)
#
# Revision 0.1  1999/06/24 11:05:01  dev
# initial revision
#
#
# DO NOT DELETE

DpRtLibrary.o: /usr/include/stdio.h /usr/include/features.h
DpRtLibrary.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
DpRtLibrary.o: /usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66/include/stddef.h
DpRtLibrary.o: /usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66/include/stdarg.h
DpRtLibrary.o: /usr/include/bits/types.h /usr/include/libio.h
DpRtLibrary.o: /usr/include/_G_config.h /usr/include/bits/stdio_lim.h
DpRtLibrary.o: /usr/include/stdlib.h /usr/include/sys/types.h
DpRtLibrary.o: /usr/include/time.h /usr/include/endian.h
DpRtLibrary.o: /usr/include/bits/endian.h /usr/include/sys/select.h
DpRtLibrary.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
DpRtLibrary.o: /usr/include/sys/sysmacros.h /usr/include/alloca.h
DpRtLibrary.o: /usr/local/jdk1.2.2/include/jni.h
DpRtLibrary.o: /usr/local/jdk1.2.2/include/linux/jni_md.h
DpRtLibrary.o: /home/dev/src/libdprt/include/dprt.h
dprt.o: /usr/include/stdio.h /usr/include/features.h /usr/include/sys/cdefs.h
dprt.o: /usr/include/gnu/stubs.h
dprt.o: /usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66/include/stddef.h
dprt.o: /usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66/include/stdarg.h
dprt.o: /usr/include/bits/types.h /usr/include/libio.h
dprt.o: /usr/include/_G_config.h /usr/include/bits/stdio_lim.h
dprt.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
dprt.o: /usr/include/endian.h /usr/include/bits/endian.h
dprt.o: /usr/include/sys/select.h /usr/include/bits/select.h
dprt.o: /usr/include/bits/sigset.h /usr/include/sys/sysmacros.h
dprt.o: /usr/include/alloca.h /usr/include/string.h /usr/include/math.h
dprt.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
dprt.o: /usr/include/bits/mathcalls.h
dprt.o: /usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66/include/float.h
dprt.o: /home/dev/src/cfitsio/fitsio.h /home/dev/src/cfitsio/longnam.h
dprt.o: /home/dev/src/libdprt/include/dprt.h
