all: whamrun mbarrun

whamrun:
	./syncin_whamrun

mbarrun:
	./syncin_mbarrun


lonestar:
	rsync -avzL --exclude=.*.swp --exclude=.*.un~ --exclude=*~ \
	  *.pbs *.sh init Makefile \
	  oo1@lonestar.tacc.utexas.edu:scratch/wham/data/1VII/

