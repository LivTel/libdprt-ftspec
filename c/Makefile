# $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/Makefile,v 0.14 2004-03-31 16:46:43 cjm Exp $

include	../../../Makefile.common
include	../../Makefile.common
include	../Makefile.common

LIBNAME		= $(LIBDPRT_HOME)_$(CCS_HOME)
INCDIR 		= $(LIBDPRT_CCS_SRC_HOME)/include
RJSINCDIR	= $(LIBDPRT_SRC_HOME)/rjs
OBJECTINCDIR	= $(LIBDPRT_SRC_HOME)/object/include
DOCSDIR 	= $(LIBDPRT_CCS_DOC_HOME)/cdocs
BINDIR		= $(LIBDPRT_CCS_BIN_HOME)/c/$(HOSTTYPE)
CFLAGS 		= -g $(CCHECKFLAG) -I$(INCDIR) -I$(CFITSIOINCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR) \
		-I$(JNIGENERALINCDIR) -I$(RJSINCDIR) -I$(OBJECTINCDIR) -L$(LT_LIB_HOME)
LINTFLAGS 	= -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
DOCFLAGS 	= -static
SRCS 		= dprt.c ngat_dprt_ccs_DpRtLibrary.c
HEADERS		= $(SRCS:%.c=%.h)
OBJS		= $(SRCS:%.c=$(BINDIR)/%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)
LIBS		= -lcfitsio -lrjs_dprt -ldprt_object -ldprt_jni_general

top: shared docs

shared: $(INCDIR)/ngat_dprt_ccs_DpRtLibrary.h $(LT_LIB_HOME)/$(LIBNAME).so

$(INCDIR)/ngat_dprt_ccs_DpRtLibrary.h: $(DPRT_BIN_HOME)/java/ngat/dprt/ccs/DpRtLibrary.class
	$(JAVAH) -jni -d $(INCDIR) -classpath $(DPRT_BIN_HOME)/java:$(CLASSPATH) ngat.dprt.ccs.DpRtLibrary
	$(TOUCH) $@

$(LT_LIB_HOME)/$(LIBNAME).so: $(OBJS)
	$(CC) $(CCSHAREDFLAG) $(CFLAGS) $(OBJS) -o $@ $(TIMELIB) $(LIBS) $(SOCKETLIB)

$(BINDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

static: $(LT_LIB_HOME)/$(LIBNAME).a

$(LT_LIB_HOME)/$(LIBNAME).a: $(OBJS)
	ar rcv $@ $?

docs: $(DOCS)

$(DOCS) : $(SRCS)
	-$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

# dont checkout ngat_dprt_ccs_DpRtLibrary.h - it is a machine built header
checkout:
	$(CO) $(CO_OPTIONS) $(SRCS)
	cd $(INCDIR); $(CO) $(CO_OPTIONS) dprt.h;

# dont checkin ngat_dprt_ccs_DpRtLibrary.h - it is a machine built header
checkin:
	-$(CI) $(CI_OPTIONS) $(SRCS)
	-(cd $(INCDIR); $(CI) $(CI_OPTIONS) dprt.h;)

staticdepend:
	makedepend -p$(BINDIR)/ -- $(CFLAGS)  -- $(SRCS)

depend:
	makedepend -p$(BINDIR)/ -- $(CFLAGS)  -- $(SRCS)

lint:
	$(LINT)	$(LINTFLAGS) $(SRCS)

clean:
	-$(RM) $(RM_OPTIONS) $(OBJS) $(LT_LIB_HOME)/$(LIBNAME).so $(LT_LIB_HOME)/$(LIBNAME).a 
	-$(RM) $(RM_OPTIONS) $(INCDIR)/ngat_dprt_ccs_DpRtLibrary.h $(TIDY_OPTIONS)

tidy:
	-$(RM)  $(RM_OPTIONS) $(TIDY_OPTIONS)
	-cd $(INCDIR); $(RM) $(RM_OPTIONS) $(TIDY_OPTIONS);

backup: tidy
	-$(RM) $(RM_OPTIONS) $(OBJS)

#
# $Log: not supported by cvs2svn $
# Revision 0.13  2004/03/31 08:56:16  cjm
# Library now repackaged as libdprt_ccs. Directory changed to ccs specific one.
#
# Revision 0.12  2004/01/30 16:44:45  cjm
# Added dprt_object library and include.
#
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

/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/stdio.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/features.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/sys/cdefs.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/gnu/stubs.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2/include/stddef.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/types.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/wordsize.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/typesizes.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/libio.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/_G_config.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/wchar.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/wchar.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/gconv.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2/include/stdarg.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/stdio_lim.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/sys_errlist.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/stdlib.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/string.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/math.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/huge_val.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/mathdef.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/mathcalls.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /home/dev/src/cfitsio/fitsio.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /home/dev/src/cfitsio/longnam.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /home/dev/src/libdprt/jni_general/include/dprt_jni_general.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /home/dev/src/libdprt/rjs/rjs_dprt.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/dirent.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /usr/include/bits/dirent.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /home/dev/src/libdprt/rjs/lt_filenames.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /home/dev/src/libdprt/rjs/dprt_libfits.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /home/dev/src/libdprt/object/include/object.h
/home/dev/bin/libdprt/ccs/c/i386-linux/dprt.o: /home/dev/src/libdprt/ccs/include/dprt.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/stdio.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/features.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/sys/cdefs.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/gnu/stubs.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2/include/stddef.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/bits/types.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/bits/wordsize.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/bits/typesizes.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/libio.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/_G_config.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/wchar.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/bits/wchar.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/gconv.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2/include/stdarg.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/bits/stdio_lim.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/bits/sys_errlist.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/stdlib.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/include/string.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/java/j2sdk1.4.1_01/include/jni.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /usr/java/j2sdk1.4.1_01/include/linux/jni_md.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /home/dev/src/libdprt/ccs/include/ngat_dprt_ccs_DpRtLibrary.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /home/dev/src/libdprt/object/include/object.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /home/dev/src/libdprt/ccs/include/dprt.h
/home/dev/bin/libdprt/ccs/c/i386-linux/ngat_dprt_ccs_DpRtLibrary.o: /home/dev/src/libdprt/jni_general/include/dprt_jni_general.h
