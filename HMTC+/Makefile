export LD_LIBRARY_PATH=./
CXX = g++
CC= gcc
LIB_DIR = lib

all:
	make -C MIRACL
	${CXX} -o test_prihmtc test_prihmtc.cpp priHMTC.cpp ${LIB_DIR}/bn_pair.cpp ${LIB_DIR}/bn_transfer.cpp ${LIB_DIR}/aes_ctr.cpp ${LIB_DIR}/WjCryptLib_Aes.c ${LIB_DIR}/WjCryptLib_AesCtr.c miracl.a -fpermissive -O2 -w
	

clean:
	rm -f *.a *.so 
	rm -f test_prihmtc


test:
	./test_prihmtc