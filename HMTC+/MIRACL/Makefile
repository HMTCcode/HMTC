all:
	cp mirdef.hpp mirdef.h
	g++ -c -m64 -g mrcore.c
	g++ -c -m64 -g mrarth0.c
	g++ -c -m64 -g mrarth1.c
	g++ -c -m64 -g mrarth2.c
	g++ -c -m64 -g mralloc.c
	g++ -c -m64 -g mrsmall.c
	g++ -c -m64 -g mrio1.c
	g++ -c -m64 -g mrio2.c
	g++ -c -m64 -g mrgcd.c
	g++ -c -m64 -g mrjack.c
	g++ -c -m64 -g mrxgcd.c
	g++ -c -m64 -g mrarth3.c
	g++ -c -m64 -g mrbits.c
	g++ -c -m64 -g mrrand.c
	g++ -c -m64 -g mrprime.c
	g++ -c -m64 -g mrcrt.c
	g++ -c -m64 -g mrscrt.c
	g++ -c -m64 -g mrmonty.c
	g++ -c -m64 -g mrpower.c
	g++ -c -m64 -g mrsroot.c
	g++ -c -m64 -g mrcurve.c
	g++ -c -m64 -g mrfast.c
	g++ -c -m64 -g mrshs.c
	g++ -c -m64 -g mrshs256.c
	g++ -c -m64 -g mrshs512.c
	g++ -c -m64 -g mrsha3.c
	g++ -c -m64 -g mrfpe.c
	g++ -c -m64 -g mraes.c
	g++ -c -m64 -g mrgcm.c
	g++ -c -m64 -g mrlucas.c
	g++ -c -m64 -g mrzzn2.c
	g++ -c -m64 -g mrzzn2b.c
	g++ -c -m64 -g mrzzn3.c
	g++ -c -m64 -g mrzzn4.c
	g++ -c -m64 -g mrecn2.c
	g++ -c -m64 -g mrstrong.c
	g++ -c -m64 -g mrbrick.c
	g++ -c -m64 -g mrebrick.c
	g++ -c -m64 -g mrec2m.c
	g++ -c -m64 -g mrgf2m.c
	g++ -c -m64 -g mrflash.c
	g++ -c -m64 -g mrfrnd.c
	g++ -c -m64 -g mrdouble.c
	g++ -c -m64 -g mrround.c
	g++ -c -m64 -g mrbuild.c
	g++ -c -m64 -g mrflsh1.c
	g++ -c -m64 -g mrpi.c
	g++ -c -m64 -g mrflsh2.c
	g++ -c -m64 -g mrflsh3.c
	g++ -c -m64 -g mrflsh4.c
	cp mrmuldv.g64 mrmuldv.c
	g++ -c -m64 -g mrmuldv.c
	g++ -c -m64 -g big.cpp
	g++ -c -m64 -g zzn.cpp
	g++ -c -m64 -g ecn.cpp
	g++ -c -m64 -g ec2.cpp
	g++ -c -m64 -g ecn3.cpp
	g++ -c -m64 -g ecn4.cpp
	g++ -c -m64 -g ecn6.cpp
	g++ -c -m64 -g ecn8.cpp
	g++ -c -m64 -g flash.cpp
	g++ -c -m64 -g crt.cpp
	g++ -c -m64 -g zzn2.cpp
	g++ -c -m64 -g zzn3.cpp
	g++ -c -m64 -g zzn4.cpp
	#g++ -c -m64 -g zzn6.cpp
	g++ -c -m64 -g zzn8.cpp
	g++ -c -m64 -g zzn6a.cpp
	g++ -c -m64 -g zzn12a.cpp
	g++ -c -m64 -g zzn18.cpp
	g++ -c -m64 -g zzn24.cpp
	g++ -c -m64 -g ecn2.cpp
	ar rc miracl.a mrcore.o mrarth0.o mrarth1.o mrarth2.o 	mralloc.o mrsmall.o mrzzn2.o mrzzn3.o
	ar r miracl.a mrio1.o mrio2.o mrjack.o mrgcd.o mrxgcd.o 	mrarth3.o mrbits.o mrecn2.o mrzzn4.o
	ar r miracl.a mrrand.o mrprime.o mrcrt.o mrscrt.o mrmonty.o 	mrcurve.o mrsroot.o mrzzn2b.o
	ar r miracl.a mrpower.o mrfast.o mrshs.o mrshs256.o mraes.o 	mrlucas.o mrstrong.o mrgcm.o    
	ar r miracl.a mrflash.o mrfrnd.o mrdouble.o mrround.o 	mrbuild.o
	ar r miracl.a mrflsh1.o mrpi.o mrflsh2.o mrflsh3.o mrflsh4.o 
	ar r miracl.a mrbrick.o mrebrick.o mrec2m.o mrgf2m.o 	mrmuldv.o mrshs512.o  mrsha3.o mrfpe.o
	ar r miracl.a big.o  zzn.o ecn.o ecn3.o ecn4.o ecn6.o ecn8.o ec2.o flash.o crt.o zzn2.o zzn4.o zzn3.o zzn6a.o zzn8.o zzn12a.o zzn18.o zzn24.o ecn2.o
	rm -f *.o	
	cp miracl.a ../
	rm -f *.a	

 
