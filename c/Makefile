# $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/Makefile,v 0.12 2004-01-30 16:44:45 cjm Exp $

include	../../Makefile.common
include	../Makefile.common

INCDIR 		= $(LIBDPRT_SRC_HOME)/include
RJSINCDIR	= $(LIBDPRT_SRC_HOME)/rjs
OBJECTINCDIR	= $(LIBDPRT_SRC_HOME)/object/include
DOCSDIR 	= $(LIBDPRT_DOC_HOME)/cdocs
BINDIR		= $(LIBDPRT_BIN_HOME)/c/$(HOSTTYPE)
CFLAGS 		= -g $(CCHECKFLAG) -I$(INCDIR) -I$(CFITSIOINCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR) -I$(RJSINCDIR) \
		-I$(OBJECTINCDIR) -L$(LT_LIB_HOME)
LINTFLAGS 	= -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
DOCFLAGS 	= -static
SRCS 		= dprt.c DpRtLibrary.c
HEADERS		= $(SRCS:%.c=%.h)
OBJS		= $(SRCS:%.c=$(BINDIR)/%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)
LIBS		= -lcfitsio -lrjs_dprt -ldprt_object

top: shared docs

shared: $(INCDIR)/DpRtLibrary.h $(LT_LIB_HOME)/libdprt.so

$(INCDIR)/DpRtLibrary.h: $(DPRT_BIN_HOME)/java/DpRtLibrary.class
	$(JAVAH) -jni -d $(INCDIR) -classpath $(DPRT_BIN_HOME)/java:$(CLASSPATH) DpRtLibrary
	$(TOUCH) $@

$(LT_LIB_HOME)/libdprt.so: $(OBJS)
	$(CC) $(CCSHAREDFLAG) $(CFLAGS) $(OBJS) -o $@ $(TIMELIB) $(LIBS) $(SOCKETLIB)

$(BINDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

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
# Revision 0.11  2002/11/26 18:13:11  cjm
# Addd lib rjs_dprt linkage.
#
# Revision 0.10  2002/05/20 11:34:45  cjm
# Moved object files to binary tree.
#
# Revision 0.9  2002/05/20 10:40:04  cjm
# Added libraries to link against.
#
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

DpRtLibrary.o: /usr/include/stdio.h /usr/include/sys/feature_tests.h
DpRtLibrary.o: /usr/include/sys/isa_defs.h /usr/include/sys/va_list.h
DpRtLibrary.o: /usr/include/stdio_tag.h /usr/include/stdio_impl.h
DpRtLibrary.o: /usr/include/stdlib.h /usr/java1.2/include/jni.h
DpRtLibrary.o: /usr/include/stdarg.h /usr/java1.2/include/solaris/jni_md.h
DpRtLibrary.o: /usr/include/sys/types.h /usr/include/sys/machtypes.h
DpRtLibrary.o: /usr/include/sys/int_types.h /usr/include/sys/select.h
DpRtLibrary.o: /usr/include/sys/time.h /usr/include/sys/time.h
DpRtLibrary.o: /usr/include/limits.h
DpRtLibrary.o: /space/home/dev/src/libdprt/include/dprt.h
dprt.o: /usr/include/stdio.h /usr/include/sys/feature_tests.h
dprt.o: /usr/include/sys/isa_defs.h /usr/include/sys/va_list.h
dprt.o: /usr/include/stdio_tag.h /usr/include/stdio_impl.h
dprt.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/math.h
dprt.o: /space/home/dev/src/cfitsio/fitsio.h /usr/include/stddef.h
dprt.o: /space/home/dev/src/cfitsio/longnam.h
dprt.o: /space/home/dev/src/libdprt/include/dprt.h
