#ifndef CURVE_H
#define CURVE_H


#include"point.h"
#include<gmpxx.h>

// elliptic this over a prime field
class EllipticCurve{
    public:
        mpz_t a, b, p;



        EllipticCurve()
        {
            mpz_init_set_ui(this->a, 0);
            mpz_init_set_ui(this->b, 0);
            mpz_init_set_ui(this->p, 0);
        }


        EllipticCurve(const char * a, const char * b, const char * p,int base)
        {
            mpz_init_set_str(this->a, a, base);
            mpz_init_set_str(this->b, b, base);
            mpz_init_set_str(this->p, p, base);
        }

        
        
        EllipticCurve(const mpz_t& a, const mpz_t& b, const mpz_t& p)
        {
            
            mpz_init_set(this->a, a);
            mpz_init_set(this->b, b);
            mpz_init_set(this->p, p);
            
        }

        EllipticCurve(const  EllipticCurve& ec)
        {
            
            mpz_init_set(this->a, ec.a);
            mpz_init_set(this->b, ec.b);
            mpz_init_set(this->p, ec.p);
        }

        EllipticCurve& operator=(const  EllipticCurve& ec)
        {
            if (this == &ec) {
                return *this; // Self-assignment check
            }
            mpz_init_set(this->a, ec.a);
            mpz_init_set(this->b, ec.b);
            mpz_init_set(this->p, ec.p);
            
            return *this;
        }

        ~EllipticCurve()
        {
             mpz_clear(this->a);
             mpz_clear(this->b);
             mpz_clear(this->p);
        }

        //We currently use the elliptic this
        //NIST P-384
        void curve_equation(mpz_t& ret, const mpz_t& x);

        // https://en.wikipedia.org/wiki/Euler%27s_criterion
        // Computes Legendre Symbol.
        bool is_quadratic_residue(const mpz_t& x);

        bool is_point_on_curve(const Point &point);

        void inv_mod_p(mpz_t& ret, const mpz_t& x);

        void ec_point_inv(Point &ret,const Point &point);

        void ec_point_double_add(Point &ret,const Point &point);

        void ec_point_add(Point &rop,const Point &op1,const Point &op2);

        void ec_point_scalar_mul(Point &rop,const mpz_t& scalar,const Point &point);
  
};


#endif