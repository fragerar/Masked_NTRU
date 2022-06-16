
CC=/usr/bin/cc
LDFLAGS=-lcrypto
CFLAGS = -O2 -march=native -DCRYPTO_NAMESPACE\(s\)=ntru_\#\#s

SOURCES = cmov.c crypto_sort_int32.c fips202.c kem.c owcpa.c pack3.c packq.c poly.c poly_lift.c poly_r2_inv.c poly_s3_inv.c rng.c sample.c sample_iid.c
SOURCES_MASKING = gadgets.c debug.c masked_owcpa.c random.c CGV14.c masked_kem.c masked_fips202.c masked_keygen.c
SOURCES_BENCH = cpucycles.c
HEADERS = api_bytes.h api.h cmov.h crypto_hash_sha3256.h crypto_int8.h crypto_sort_int32.h fips202.h kem.h owcpa.h params.h poly.h poly_r2_inv.h rng.h sample.h
OBJS = poly_mod_3_Phi_n.s poly_mod_q_Phi_n.s poly_r2_mul.s poly_rq_mul.s poly_rq_to_s3.s square_126_509_shufbytes.s square_1_509_patience.s square_15_509_shufbytes.s square_252_509_shufbytes.s square_30_509_shufbytes.s square_3_509_patience.s square_63_509_shufbytes.s square_6_509_patience.s vec32_sample_iid.s

ORDER=1
RNG=1
MACRO = -D MASKING_ORDER=$(ORDER) -D RNG_MODE=$(RNG)


test:
	$(CC) $(CFLAGS) $(SOURCES) $(SOURCES_MASKING) $(OBJS) test.c $(LDFLAGS) -o test

bench:
	$(CC) $(CFLAGS) $(SOURCES) $(SOURCES_MASKING) $(SOURCES_BENCH) $(MACRO) $(OBJS) benchmarks.c $(LDFLAGS) -o bench

.PHONY: test bench clean 

clean:
	-rm -f test
	-rm -f bench
