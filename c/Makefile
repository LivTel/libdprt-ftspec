# $Header: /space/home/eng/cjm/cvs/libdprt-ftspec/c/Makefile,v 0.2 1999-06-25 09:50:59 dev Exp $

include			../../Makefile.common

INCDIR = ../include
DOCSDIR = ../cdocs
JAVADIR = ../java
CFLAGS = -g $(CCHECKFLAG) -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
LINTFLAGS = -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
DOCFLAGS = -static
SRCS 		= 	DpRtLibrary.c dprt.c
HEADERS		=	$(SRCS:%.c=%.h)
OBJS		=	$(SRCS:%.c=%.o)
DOCS 		= 	$(SRCS:%.c=$(DOCSDIR)/%.html)

top: shared

shared: $(JAVADIR)/libdprt.so $(DOCS)

$(JAVADIR)/libdprt.so: $(OBJS)
	cc $(CCSHAREDFLAG) $(CFLAGS) $(OBJS) -o $@ $(TIMELIB)

static: libdprt.a $(DOCS)

libdprt.a: $(OBJS)
	ar rcv $@ $?

$(DOCS) : $(OBJS)
	$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

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
	-$(RM) $(RM_OPTIONS) $(OBJS) $(JAVADIR)/libdprt.so libdprt.a $(TIDY_OPTIONS)

tidy:
	-$(RM)  $(RM_OPTIONS) $(TIDY_OPTIONS)
	-cd $(INCDIR); $(RM) $(RM_OPTIONS) $(TIDY_OPTIONS);

#
# $Log: not supported by cvs2svn $
# Revision 0.1  1999/06/24 11:05:01  dev
# initial revision
#
#
# DO NOT DELETE

