#ifndef PEER_H
#define PEER_H

#include <string>
#include <gmpxx.h>
#include "curve.h"

/* 
"""
    Implements https://wlan1nde.wordpress.com/2018/09/14/wpa3-improving-your-wlan-security/
    Take a ECC curve from here: https://safecurves.cr.yp.to/

    Example: NIST P-384
    y^2 = x^3-3x+27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575
    modulo p = 2^384 - 2^128 - 2^96 + 2^32 - 1
    2000 NIST; also in SEC 2 and NSA Suite B

    See here: https://www.rfc-editor.org/rfc/rfc5639.txt

   Curve-ID: brainpoolP256r1
      p =
      A9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377
      A =
      7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9
      B =
      26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6
      x =
      8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262
      y =
      547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997
      q =
      A9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7
      h = 1
    """
 */
class Peer{
    public:
        
        std::string password;
        std::string mac_address;
        std::string name;
        

        mpz_t a;
        mpz_t b;
        mpz_t p;
        mpz_t q;

        mpz_t priv;
        mpz_t mask;
        mpz_t scalar;
        

        EllipticCurve curve;
        Point PE; //Password Element
        Point element;
        Point k;

        Point peer_element;
        mpz_t peer_scalar;
        std::string peer_mac;
        
        
        Peer(const std::string& password,const std::string& mac_address,const std::string& name);
        ~Peer();

        void initiate(std::string peer_mac,int k = 40);
        std::string compute_hashed_password(const int& counter);
        void key_derivation_function(mpz_t& rop,const int& n,const std::string& base,const std::string & seed);
        void commit_exchange();
        std::string compute_shared_secret(const Point& peer_element,const mpz_t& peer_scalar,const std::string& peer_mac);
        std::string confirm_exchange(const std::string& peer_token);
};

#endif