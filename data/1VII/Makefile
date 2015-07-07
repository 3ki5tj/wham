all: whamrun mbarrun

whamrun:
	cat whamrun1/xvg_wham.log         >  xvg_wham.log
	cat whamrun2/xvg_wham.log         >> xvg_wham.log
	cat whamrun3/xvg_wham.log         >> xvg_wham.log
	cat whamrun4/xvg_wham.log         >> xvg_wham.log
	cat whamrun1/xvg_whamtm.log       >  xvg_whamtm.log
	cat whamrun2/xvg_whamtm.log       >> xvg_whamtm.log
	cat whamrun3/xvg_whamtm.log       >> xvg_whamtm.log
	cat whamrun4/xvg_whamtm.log       >> xvg_whamtm.log
	./stat.py xvg_wham*.log

mbarrun:
	cat mbarrun1/xvg_mbar.log         >  xvg_mbar.log
	cat mbarrun2/xvg_mbar.log         >> xvg_mbar.log
	cat mbarrun3/xvg_mbar.log         >> xvg_mbar.log
	cat mbarrun4/xvg_mbar.log         >> xvg_mbar.log
	cat mbarrun5/xvg_mbar.log         >  xvg_mbar.log
	cat mbarrun6/xvg_mbar.log         >> xvg_mbar.log
	cat mbarrun7/xvg_mbar.log         >> xvg_mbar.log
	cat mbarrun8/xvg_mbar.log         >> xvg_mbar.log
	cat mbarrun1/xvg_mbartm.log       >  xvg_mbartm.log
	cat mbarrun2/xvg_mbartm.log       >> xvg_mbartm.log
	cat mbarrun3/xvg_mbartm.log       >> xvg_mbartm.log
	cat mbarrun4/xvg_mbartm.log       >> xvg_mbartm.log
	cat mbarrun5/xvg_mbartm.log       >  xvg_mbartm.log
	cat mbarrun6/xvg_mbartm.log       >> xvg_mbartm.log
	cat mbarrun7/xvg_mbartm.log       >> xvg_mbartm.log
	cat mbarrun8/xvg_mbartm.log       >> xvg_mbartm.log
	./stat.py xvg_mbar*.log


lonestar:
	rsync -avzL --exclude=.*.swp --exclude=.*.un~ --exclude=*~ \
	  *.pbs *.sh init \
	  oo1@lonestar.tacc.utexas.edu:scratch/wham/data/1VII/

