# $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/Makefile,v 0.4 1999-08-24 16:05:29 cjm Exp $

include	../../Makefile.common
include	../Makefile.common

LIBDIR 		= $(LT_BIN_HOME)/lib
INCDIR 		= $(DPRT_SRC_HOME)/include
DOCSDIR 	= $(DPRT_DOC_HOME)/cdocs
CFLAGS 		= -g $(CCHECKFLAG) -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
LINTFLAGS 	= -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
DOCFLAGS 	= -static
SRCS 		= DpRtLibrary.c dprt.c
HEADERS		= $(SRCS:%.c=%.h)
OBJS		= $(SRCS:%.c=%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)

top: shared docs

shared: $(LIBDIR)/libdprt.so

$(LIBDIR)/libdprt.so: $(OBJS)
	cc $(CCSHAREDFLAG) $(CFLAGS) $(OBJS) -o $@ $(TIMELIB)

static: $(LIBDIR)/libdprt.a

$(LIBDIR)/libdprt.a: $(OBJS)
	ar rcv $@ $?

docs: $(DOCS)

$(DOCS) : $(SRCS)
	-$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

checkout:
	$(CO) $(CO_OPTIONS) $(SRCS)
	cd $(INCDIR); $(CO) $(CO_OPTIONS) $(HEADERS);

checkin:
	-$(CI) $(SRCS)
	-(cd $(INCDIR); $(CI) $(HEADERS);)

staticdepend:
	makedepend -- $(CFLAGS) -- $(SRCS)

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

lint:
	$(LINT)	$(LINTFLAGS) $(SRCS)

clean:
	-$(RM) $(RM_OPTIONS) $(OBJS) $(LIBDIR)/libdprt.so $(LIBDIR)/libdprt.a $(TIDY_OPTIONS)

tidy:
	-$(RM)  $(RM_OPTIONS) $(TIDY_OPTIONS)
	-cd $(INCDIR); $(RM) $(RM_OPTIONS) $(TIDY_OPTIONS);

#
# $Log: not supported by cvs2svn $
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

