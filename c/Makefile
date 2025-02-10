# $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/Makefile,v 0.15 2025-02-10 11:52:00 cjm Exp $

include	../../../Makefile.common
include	../../Makefile.common
include	../Makefile.common

LIBNAME		= $(LIBDPRT_HOME)_$(FTSPEC_HOME)
INCDIR 		= $(LIBDPRT_FTSPEC_SRC_HOME)/include
DOCSDIR 	= $(LIBDPRT_FTSPEC_DOC_HOME)/cdocs
BINDIR		= $(LIBDPRT_FTSPEC_BIN_HOME)/c/$(HOSTTYPE)
CFLAGS 		= -g $(CCHECKFLAG) -I$(INCDIR) -I$(CFITSIOINCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR) \
		-I$(JNIGENERALINCDIR) -L$(LT_LIB_HOME)
LINTFLAGS 	= -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
DOCFLAGS 	= -static
SRCS 		= dprt.c ngat_dprt_ftspec_DpRtLibrary.c
HEADERS		= $(SRCS:%.c=%.h)
OBJS		= $(SRCS:%.c=$(BINDIR)/%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)
LIBS		= -lcfitsio -ldprt_jni_general

top: shared docs

shared: $(INCDIR)/ngat_dprt_ftspec_DpRtLibrary.h $(LT_LIB_HOME)/$(LIBNAME).so

$(INCDIR)/ngat_dprt_ftspec_DpRtLibrary.h: $(DPRT_BIN_HOME)/java/ngat/dprt/ftspec/DpRtLibrary.class
	$(JAVAH) -jni -d $(INCDIR) -classpath $(DPRT_BIN_HOME)/java:$(CLASSPATH) ngat.dprt.ftspec.DpRtLibrary
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

# dont checkout ngat_dprt_ftspec_DpRtLibrary.h - it is a machine built header
checkout:
	$(CO) $(CO_OPTIONS) $(SRCS)
	cd $(INCDIR); $(CO) $(CO_OPTIONS) dprt.h;

# dont checkin ngat_dprt_ftspec_DpRtLibrary.h - it is a machine built header
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
	-$(RM) $(RM_OPTIONS) $(INCDIR)/ngat_dprt_ftspec_DpRtLibrary.h $(TIDY_OPTIONS)

tidy:
	-$(RM)  $(RM_OPTIONS) $(TIDY_OPTIONS)
	-cd $(INCDIR); $(RM) $(RM_OPTIONS) $(TIDY_OPTIONS);

backup: tidy
	-$(RM) $(RM_OPTIONS) $(OBJS)

#
# $Log: not supported by cvs2svn $
#

