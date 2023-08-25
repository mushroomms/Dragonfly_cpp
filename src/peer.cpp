#include <gmpxx.h>
#include <openssl/sha.h>
#include <string>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <ctime>
#include <sodium.h>
#include "../include/peer.h"
#include "../include/mpz_math.h"

Peer::Peer(const std::string &password, const std::string &mac_address, const std::string &name)
{
    this->password = password;
    this->mac_address = mac_address.c_str();
    this->name = name;

    //Try out Curve-ID: brainpoolP256t1
    mpz_init_set_str(this->p, "A9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377", 16);
    mpz_init_set_str(this->a, "7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9", 16);
    mpz_init_set_str(this->b, "26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6", 16);
    mpz_init_set_str(this->q, "A9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7", 16);

    mpz_init(this->priv);
    mpz_init(this->mask);
    mpz_init(this->scalar);
    mpz_init(this->peer_scalar);

    this->curve = EllipticCurve(this->a, this->b, this->p);

    //gmp_printf("hello:p = %Zd,curve.p = %Zd   \n",this->p,this->curve.p);
}

Peer::~Peer()
{
    mpz_clear(this->p);
    mpz_clear(this->a);
    mpz_clear(this->b);
    mpz_clear(this->q);

    mpz_clear(this->priv);
    mpz_clear(this->mask);
    mpz_clear(this->scalar);
    mpz_clear(this->peer_scalar);
}

void Peer::initiate(std::string peer_mac, int k)
{

    /* 
    """
        See algorithm in https://tools.ietf.org/html/rfc7664
        in section 3.2.1
    """
     */
    this->peer_mac = peer_mac.c_str();

    int found = 0;
    int num_valid_points = 0;
    int n = mpz_sizeinbase(this->p, 2) + 64;
    std::string str_save;
    const std::string str_for_seed = "Dragonfly Hunting And Pecking";

    mpz_t x, y;
    mpz_init(x), mpz_init(y);

    mpz_t temp, seed, val;
    mpz_init(temp), mpz_init(seed), mpz_init(val);
    for (int counter = 1; counter <= k || found == 0; counter++)
    {
        std::string base = this->compute_hashed_password(counter);
        //std::cout << "base = " << str2hex_str(base) << std::endl;
        this->key_derivation_function(temp, n, base, str_for_seed);

       // gmp_printf("temp = %Zd\n", temp);

        mpz_sub_ui(seed, this->p, 1);
        mpz_mod(seed, temp, seed);
        mpz_add_ui(seed, seed, 1);
        //gmp_printf("seed = %Zd\n", seed);

        this->curve.curve_equation(val, seed);
        //gmp_printf("val = %Zd\n", val);
        if (this->curve.is_quadratic_residue(val))
        {
            mpz_set(x, seed);
            str_save = base;
            found = 1;
            num_valid_points++;
            // std::cout << "["  << this->name << "]" << " Got point after " << counter << " iterations" << std::endl;
        }
    }
    assert(found == 1 && "No valid point found");

    this->curve.curve_equation(temp, x);
    tonelli_shanks(y, temp, this->p);

    mpz_t save;
    mpz_init_set_str(save, str2hex_str(str_save).c_str(), 16);
    if ((mpz_odd_p(y) && mpz_odd_p(save)) || (mpz_even_p(y) && mpz_even_p(save)))
    {
        this->PE.set_x(x);
        this->PE.set_y(y);
    }
    else
    {
        mpz_sub(y, this->p, y);
        this->PE.set_x(x);
        this->PE.set_y(y);
    }

    //check valid point
    this->curve.curve_equation(temp, this->PE.x);
    mpz_powm_ui(val, this->PE.y, 2, this->p);
    assert(mpz_cmp(temp, val) == 0);

    gmp_printf("[%s] Using %i-th valid point = (x=%Zd, y=%Zd) \n", this->name.c_str(), num_valid_points, this->PE.x, this->PE.y);

    mpz_clear(temp), mpz_clear(seed), mpz_clear(val);
    mpz_clear(save);
    mpz_clear(x), mpz_clear(y);
}

std::string Peer::compute_hashed_password(const int &counter)
{

    std::string min_message = std::min(this->mac_address, this->peer_mac);
    std::string max_message = std::max(this->mac_address, this->peer_mac);
    std::string message = max_message + min_message + this->password + std::to_string(counter);

    std::string digest = sha256(message);
    return digest;
}

void Peer::key_derivation_function(mpz_t &rop, const int &n, const std::string &base, const std::string &seed)
{

    std::string combined_seed = base + seed;
    std::string hex_combined_seed = str2hex_str(combined_seed);

    mpz_t mpz_hex_combined_seed;
    mpz_init_set_str(mpz_hex_combined_seed, hex_combined_seed.c_str(), 16);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed(state, mpz_hex_combined_seed);
    mpz_urandomb(rop, state, n);

    // gmp_printf("mpz_urandomb :rop = %Zd ,mpz_sizeinbase:%d \n",
    //             rop,mpz_sizeinbase(rop,2));
    mpz_clear(mpz_hex_combined_seed);
}

void Peer::commit_exchange()
{

    /* 
        """
        This is basically Diffie Hellman Key Exchange (or in our case ECCDH)

        In the Commit Exchange, both sides commit to a single guess of the
        password.  The peers generate a scalar and an element, exchange them
        with each other, and process the other's scalar and element to
        generate a common and shared secret.

        If we go back to elliptic curves over the real numbers, there is a nice geometric
        interpretation for the ECDLP: given a starting point P, we compute 2P, 3P, . . .,
        d P = T , effectively hopping back and forth on the elliptic curve. We then publish
        the starting point P (a public parameter) and the final point T (the public key). In
        order to break the cryptosystem, an attacker has to figure out how often we “jumped”
        on the elliptic curve. The number of hops is the secret d, the private key.
        """
        # seed the PBG before picking a new random number
        # random.seed(time.process_time())

        # None or no argument seeds from current time or from an operating
        # system specific randomness source if available. 
    */

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, randombytes_random());

    //# Otherwise, each party chooses two random numbers, private and mask
    mpz_urandomm(this->priv, state, this->p);
    if (mpz_cmp_ui(this->priv, 0) == 0)
    {
        mpz_set_ui(this->priv, 1);
    }
    mpz_urandomm(this->mask, state, this->p);
    if (mpz_cmp_ui(this->mask, 0) == 0)
    {
        mpz_set_ui(this->mask, 1);
    }
    gmp_printf("[%s] Private = %Zd \n", this->name.c_str(), this->priv);
    gmp_printf("[%s] Mask = %Zd \n", this->name.c_str(), this->mask);
    
    /* 
        # These two secrets and the Password Element are then used to construct
        # the scalar and element:

        # what is q?
        # o  A point, G, on the elliptic curve, which serves as a generator for
        #    the ECC group.  G is chosen such that its order, with respect to
        #    elliptic curve addition, is a sufficiently large prime.
        #
        # o  A prime, q, which is the order of G, and thus is also the size of
        #    the cryptographic subgroup that is generated by G.

        # https://math.stackexchange.com/questions/331329/is-it-possible-to-compute-order-of-a-point-over-elliptic-curve
        # In the elliptic Curve cryptography, it is said that the order of base point
        # should be a prime number, and order of a point P is defined as k, where kP=O.

        # Theorem 9.2.1 The points on an elliptic curve together with O
        # have cyclic subgroups. Under certain conditions all points on an
        # elliptic curve form a cyclic group.
        # For this specific curve the group order is a prime and, according to Theo-
        # rem 8.2.4, every element is primitive.

        # Question: What is the order of our PE?
        # the order must be p, since p is a prime
     */

    // self.scalar = (self.private + self.mask) % self.q
    mpz_add(this->scalar, this->priv, this->mask);
    mpz_mod(this->scalar, this->scalar, this->q);


    /* 
        # If the scalar is less than two (2), the private and mask MUST be
        # thrown away and new values generated.  Once a valid scalar and
        # Element are generated, the mask is no longer needed and MUST be
        # irretrievably destroyed. 
    */
    assert(mpz_cmp_ui(this->scalar, 2) >= 0);

    Point res_point;
    this->curve.ec_point_scalar_mul(res_point, this->mask, this->PE);

    /* 
        # get the inverse of res
        # −P = (x_p , p − y_p ).
    */
    this->curve.ec_point_inv(this->element, res_point);

    assert(this->curve.is_point_on_curve(this->element));

    gmp_printf("[%s] Scalar = %Zd \n", this->name.c_str(), this->scalar);
    gmp_printf("[%s] Element.x = %Zd, Element.y = %Zd \n", this->name.c_str(),
               this->element.x, this->element.y);
}

std::string Peer::compute_shared_secret(const Point &peer_element, const mpz_t &peer_scalar, const std::string &peer_mac)
{
    /* 
    """
        ss = F(scalar-op(private,
                         element-op(peer-Element,
                                    scalar-op(peer-scalar, PE))))

        AP1: K = private(AP1) • (scal(AP2) • P(x, y) ◊ new_point(AP2))
               = private(AP1) • private(AP2) • P(x, y)
        AP2: K = private(AP2) • (scal(AP1) • P(x, y) ◊ new_point(AP1))
               = private(AP2) • private(AP1) • P(x, y)

        A shared secret element is computed using one’s rand and
        the other peer’s element and scalar:
        Alice: K = rand A • (scal B • PW + elemB )
        Bob: K = rand B • (scal A • PW + elemA )

        Since scal(APx) • P(x, y) is another point, the scalar multiplied point
        of e.g. scal(AP1) • P(x, y) is added to the new_point(AP2) and afterwards
        multiplied by private(AP1).
    """ 
    */
    this->peer_element = peer_element;
    mpz_set(this->peer_scalar,peer_scalar);
    this->peer_mac = peer_mac.c_str();

    gmp_printf("[%s] Peer Scalar = %Zd \n", this->name.c_str(), this->peer_scalar);
    gmp_printf("[%s] Peer Element.x = %Zd, Peer Element.y = %Zd \n", this->name.c_str(),
               peer_element.x, peer_element.y);

    assert(this->curve.is_point_on_curve(peer_element));
    /* 
        # If both the peer-scalar and Peer-Element are
        # valid, they are used with the Password Element to derive a shared
        # secret, ss:
     */
    Point z,zz;
    this->curve.ec_point_scalar_mul(z,peer_scalar, this->PE);
    this->curve.ec_point_add(zz, peer_element, z);
    this->curve.ec_point_scalar_mul(this->k, this->priv, zz);

    gmp_printf("[%s] Shared Secret ss = %Zd \n", this->name.c_str(), this->k.x);

    std::string own_message;
    own_message += mpz2dec_str(this->k.x);
    own_message += mpz2dec_str(this->scalar);
    own_message += mpz2dec_str(this->peer_scalar);
    own_message += mpz2dec_str(this->element.x);
    own_message += mpz2dec_str(this->peer_element.x);
    own_message += this->mac_address.c_str();

    std::string token = sha256(own_message);
    printf("[%s] Token = %s \n", this->name.c_str(), str2hex_str(token).c_str() );
    return token;
}


std::string Peer::confirm_exchange(const std::string& peer_token){

    std::string peer_message;
    peer_message += mpz2dec_str(this->k.x);
    peer_message += mpz2dec_str(this->peer_scalar);
    peer_message += mpz2dec_str(this->scalar);
    peer_message += mpz2dec_str(this->peer_element.x);
    peer_message += mpz2dec_str(this->element.x);
    peer_message += this->peer_mac.c_str();

    std::string peer_token_computed = sha256(peer_message);

    printf("\n[%s] Computed Token from Peer = %s \n", this->name.c_str(), str2hex_str(peer_token_computed).c_str());
    printf("[%s] Received Token from Peer = %s \n", this->name.c_str(), str2hex_str(peer_token).c_str());

    // Use assert to check if the tokens are the same
    assert(peer_token == peer_token_computed && "Peer token mismatch!");

    std::string pmk_message;
    pmk_message += mpz2dec_str(k.x);

    mpz_t val;
    mpz_init(val);
    mpz_add(val, this->scalar, this->peer_scalar);
    mpz_mod(val, val, this->q);
    pmk_message += mpz2dec_str(val);
    mpz_clear(val);

    std::string pmk = sha256(pmk_message);
    printf("[%s] Pairwise Master Key(PMK) = %s \n", this->name.c_str(), str2hex_str(pmk).c_str());

    return pmk;
}