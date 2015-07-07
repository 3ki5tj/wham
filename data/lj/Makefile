all: mkrunev mkrunevkth mkrunevhp mkrunevhpl

mkrunev:
	cat runev1/lj2.log   >  lj2.log
	cat runev2/lj2.log   >> lj2.log
	cat runev3/lj2.log   >> lj2.log
	cat runev4/lj2.log   >> lj2.log
	cat runev1/lj2tm.log >  lj2tm.log
	cat runev2/lj2tm.log >> lj2tm.log
	cat runev3/lj2tm.log >> lj2tm.log
	cat runev4/lj2tm.log >> lj2tm.log
	./stat.py lj2.log lj2tm.log

mkrunevkth:
	cat runevkth1/lj2kth.log   >  lj2kth.log
	cat runevkth2/lj2kth.log   >> lj2kth.log
	cat runevkth3/lj2kth.log   >> lj2kth.log
	cat runevkth4/lj2kth.log   >> lj2kth.log
	cat runevkth1/lj2kthtm.log >  lj2kthtm.log
	cat runevkth2/lj2kthtm.log >> lj2kthtm.log
	cat runevkth3/lj2kthtm.log >> lj2kthtm.log
	cat runevkth4/lj2kthtm.log >> lj2kthtm.log
	./stat.py lj2kth.log lj2kthtm.log

mkrunevhp:
	cat runevhp1/lj2hp.log   >  lj2hp.log
	cat runevhp2/lj2hp.log   >> lj2hp.log
	cat runevhp3/lj2hp.log   >> lj2hp.log
	cat runevhp4/lj2hp.log   >> lj2hp.log
	cat runevhp1/lj2hptm.log >  lj2hptm.log
	cat runevhp2/lj2hptm.log >> lj2hptm.log
	cat runevhp3/lj2hptm.log >> lj2hptm.log
	cat runevhp4/lj2hptm.log >> lj2hptm.log
	./stat.py lj2hp.log lj2hptm.log

mkrunevhpl:
	cat runevhpl1/lj2hpl.log   >  lj2hpl.log
	cat runevhpl2/lj2hpl.log   >> lj2hpl.log
	cat runevhpl3/lj2hpl.log   >> lj2hpl.log
	cat runevhpl4/lj2hpl.log   >> lj2hpl.log
	cat runevhpl1/lj2hpltm.log >  lj2hpltm.log
	cat runevhpl2/lj2hpltm.log >> lj2hpltm.log
	cat runevhpl3/lj2hpltm.log >> lj2hpltm.log
	cat runevhpl4/lj2hpltm.log >> lj2hpltm.log
	./stat.py lj2hpl.log lj2hpltm.log

lonestar:
	make -C .. lonestar
