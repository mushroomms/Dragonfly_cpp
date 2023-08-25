#include <gmpxx.h>
#include <assert.h>
#include <openssl/sha.h>
#include <openssl/evp.h>
#include "../include/mpz_math.h"


void legendre(mpz_t& ret, const mpz_t& a, const mpz_t& p)
{
    mpz_t v;
    mpz_init(v);
    mpz_sub_ui(v, p, 1);    //v = p - 1
    mpz_fdiv_q_ui(v, v, 2); //v = v // 2
    mpz_powm(ret, a, v, p); // ret =  pow(a, (p - 1) // 2, p)
    mpz_clear(v);
}


//https://github.com/hvy/integer-factorization/blob/80e73f8910b827d7d17221f829426a15088648c7/src/test.c
void tonelli_shanks(mpz_t& x, const mpz_t& n, const mpz_t& p)
{
    mpz_t v;
    mpz_init(v);
    legendre(v, n, p);
    bool is_square_mod_p = (mpz_cmp_ui(v, 1) == 0);
    mpz_clear(v);
    if (!is_square_mod_p)
    {
        gmp_printf("square = %Zd and p =  %Zd is not a square (mod p)  \n", n, p);
        assert(is_square_mod_p);
    }

    mpz_t two, minus_one, p_minus_one, q, s, z, c, r, t, m, tmp, tmp2, t_sqrt, b;

    mpz_init_set_ui(two, 2);
    mpz_init_set_ui(minus_one, -1);
    mpz_init(p_minus_one);
    mpz_init(q);
    mpz_init_set_ui(s, 0);
    mpz_init_set_ui(z, 2);
    mpz_init(c);
    mpz_init(r);
    mpz_init(t);
    mpz_init(m);
    mpz_init(tmp);
    mpz_init(tmp2);
    mpz_init(t_sqrt);
    mpz_init(b);

    mpz_sub_ui(p_minus_one, p, 1);
    mpz_set(q, p_minus_one);
    while (0 != mpz_even_p(q) /* while q is even */)
    {
        mpz_add_ui(s, s, 1);
        mpz_divexact_ui(q, q, 2);
    }

    /* Find the quadratic non-residue */
    while (-1 != mpz_legendre(z, p))
    {
        mpz_add_ui(z, z, 1);
    }

    mpz_powm(c, z, q, p);
    mpz_add_ui(tmp, q, 1); /* q is odd, so tmp will be even */
    mpz_divexact_ui(tmp, tmp, 2);
    mpz_powm(r, n, tmp, p);
    mpz_powm(t, n, q, p);
    mpz_set(m, s);

    while (1)
    {
        if (0 == mpz_cmp_ui(t, 1) /* t == 1 */)
        {
            mpz_set(x, r);
            break;
        }
        else
        {
            long msl = mpz_get_si(m);
            for (long i = 0; i < msl; ++i)
            {
                mpz_pow_ui(tmp, two, i);
                mpz_powm(t_sqrt, t, tmp, p);
                if (0 == mpz_cmp_ui(t_sqrt, 1) /* t^(1/2) == 1 */)
                {
                    mpz_pow_ui(tmp, two, msl - i - 1);
                    mpz_powm(b, c, tmp, p);
                    mpz_mul(r, r, b);
                    mpz_mod(r, r, p);
                    mpz_pow_ui(tmp, b, 2);
                    mpz_mul(t, t, tmp);
                    mpz_mod(t, t, p);
                    mpz_mod(c, tmp, p);
                    mpz_set_ui(m, i);
                    break;
                }
            }
        }
    }

    mpz_clear(two);
    mpz_clear(minus_one);
    mpz_clear(p_minus_one);
    mpz_clear(q);
    mpz_clear(s);
    mpz_clear(z);
    mpz_clear(c);
    mpz_clear(r);
    mpz_clear(t);
    mpz_clear(m);
    mpz_clear(tmp);
    mpz_clear(tmp2);
    mpz_clear(t_sqrt);
    mpz_clear(b);
}


std::string str2hex_str(const std::string& op)
{
	unsigned char hexChars[] = "0123456789abcdef";
	unsigned int a, b;
    std::string hex_rop;

	for (int i = 0; i < op.size(); ++i)
	{
        unsigned char c = op[i];
		a = (c >> 4) & 0x0f;
		b = c & 0x0f;
		hex_rop.push_back(hexChars[a]);
        hex_rop.push_back(hexChars[b]);
	}

    return hex_rop;
}

std::string mpz2dec_str(const mpz_t& op)
{
	
    int base = 10;
    unsigned int  len = mpz_sizeinbase(op, base) + 2;
    char* rop = new char[len];
	mpz_get_str(rop, base, op);
    std::string dec_str(rop);
    delete(rop);
    return dec_str;
}

std::string sha256(const std::string& message)
{
    EVP_MD_CTX *ctx;
    unsigned char hash[SHA256_DIGEST_LENGTH] = {0};
    unsigned int msg_size = message.size();

    ctx = EVP_MD_CTX_new();
    EVP_DigestInit_ex(ctx, EVP_sha256(), NULL);
    EVP_DigestUpdate(ctx, message.c_str(), msg_size);
    EVP_DigestFinal_ex(ctx, hash, &msg_size);

    std::string digest((char *)hash, 32);
    return digest;
}