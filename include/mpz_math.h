#ifndef MPZ_MATH
#define MPZ_MATH

#include<gmpxx.h>

void legendre(mpz_t& ret, const mpz_t& a, const mpz_t& p);


//https://github.com/hvy/integer-factorization/blob/80e73f8910b827d7d17221f829426a15088648c7/src/test.c
void tonelli_shanks(mpz_t& x, const mpz_t& n, const mpz_t& p);


std::string str2hex_str(const std::string& op);

std::string mpz2dec_str(const mpz_t& op);

std::string sha256(const std::string& message);

#endif