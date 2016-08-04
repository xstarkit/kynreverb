SRC_KYNREFREV=xside_threads.c xside.c libxspec.c fft_reverberation.c xskynrefrev.c

CC=gcc
INCLUDE=-I/usr/local/share/heasoft/x86_64-unknown-linux-gnu-libc2.19/include
LIBRARY_PATH=-L/usr/local/share/heasoft/x86_64-unknown-linux-gnu-libc2.19/lib
LIBRARY= -lcfitsio_3.37 -lpthread
CFLAGS=-fPIC -O3 -Wall -DOUTSIDE_XSPEC -lm

kynrefrev: $(SRC_KYNREFREV)
		$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC_KYNREFREV) $(LIBRARY_PATH) $(LIBRARY)
