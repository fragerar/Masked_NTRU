
CC=/usr/bin/cc
CFLAGS  =  
LDFLAGS=-lcrypto
CFLAGS = -DCRYPTO_NAMESPACE\(s\)=ntru_\#\#s -O3 -fomit-frame-pointer -march=native

SOURCES = cmov.c crypto_sort_int32.c fips202.c kem.c owcpa.c pack3.c packq.c poly.c poly_lift.c poly_mod.c poly_r2_inv.c poly_rq_mul.c poly_s3_inv.c rng.c sample.c sample_iid.c 
SOURCES_MASKING = gadgets.c debug.c masked_owcpa.c random.c CGV14.c masked_kem.c masked_fips202.c masked_keygen.c
SOURCES_BENCH = cpucycles.c
HEADERS = api_bytes.h api.h cmov.h crypto_hash_sha3256.h crypto_sort_int32.h fips202.h kem.h owcpa.h params.h poly.h rng.h sample.h random.h
HEADERS_MASKING = masking.h

ORDER=1
RNG=1
MACRO = -D MASKING_ORDER=$(ORDER) -D RNG_MODE=$(RNG)

test:
	$(CC) $(CFLAGS) $(SOURCES) $(SOURCES_MASKING) test.c $(LDFLAGS) -o test

bench:
	$(CC) $(CFLAGS) $(SOURCES) $(SOURCES_MASKING) $(SOURCES_BENCH) $(MACRO) benchmarks.c $(LDFLAGS) -o bench

.PHONY: test bench clean 

clean:
	-rm -f test