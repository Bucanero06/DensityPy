SUBDIRS=libGBTO
PWD=$(shell pwd)

install:
	echo $(PWD)
	@for i in $(SUBDIRS) ; \
	do \
		cd $$i ; \
		echo "......Entering: $$i" ; \
		make $@ ; \
		cd .. ; \
		echo "......Leaving: $$i" ; echo; echo;\
		cp $$PWD/$$i/$$i.a ${LIBDIR}/$$i.a; \
		cp $$i/*.mod ../inc/ ; \
	done

clean:
	@for i in $(SUBDIRS); \
        do \
        cd $$i ; \
        echo "Current directory: $$i" ; \
                $(MAKE) $@ ; \
                cd .. ; \
        done

