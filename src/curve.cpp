#include<gmpxx.h>
#include<assert.h>
#include "../include/mpz_math.h"
#include "../include/curve.h"

//(pow(x, 3) + (self.a * x) + self.b) % self.p
// (x ** 3 + self.a * x + self.b) % self.p
void EllipticCurve::curve_equation(mpz_t& ret, const mpz_t& x){
    mpz_t v1,v2;
    mpz_init(v1), mpz_init(v2); 
    mpz_pow_ui(v1, x, 3); // v1 =  pow(x, 3)
    mpz_mul(v2, this->a, x); // v2 =  this->a * x
    mpz_add(ret, v1, v2); // ret = v1 + v2
    mpz_add(ret, ret, this->b); // ret = ret + this->b
    mpz_mod(ret, ret, this->p); // ret = ret % this->p
    mpz_clear(v1);mpz_clear(v2);
}

bool EllipticCurve::is_quadratic_residue(const mpz_t& x){
    mpz_t leg_ret;
    mpz_init(leg_ret);
    legendre(leg_ret, x, this->p);
    bool is_QR = (mpz_cmp_ui(leg_ret,1) == 0);
    mpz_clear(leg_ret);

    return is_QR;
}


bool EllipticCurve::is_point_on_curve(const Point &point){

    if(mpz_cmp_ui(point.x,0) == 0 && mpz_cmp_ui(point.y,0) == 0){
        return true;
    }

    /* 
        x, y, = P[0], P[1]
        left = y * y
        right = (x * x * x) + (self.a * x) + self.b
        return (left - right) % self.p == 0 
    */
    mpz_t left, right, v;
    mpz_init(left), mpz_init(right),mpz_init(v);
    mpz_pow_ui(left,point.y,2); // left = y * y

    mpz_pow_ui(right,point.x,3); // right = (x * x * x)
    mpz_mul(v, this->a, point.x);// v = a * x
    mpz_add(right, right, v);// right = right + v
    mpz_add(right, right, this->b);// right = right + v

    mpz_sub(v, left, right); // v = left - right
    mpz_mod(v, v, this->p); // ret = ret % this->p

    bool is_on_curve =  (mpz_cmp_ui(v,0) == 0);

    mpz_clear(left),  mpz_clear(right), mpz_clear(v);
    return is_on_curve;    
}


void EllipticCurve::inv_mod_p(mpz_t& ret, const mpz_t& x)
{ 
    /* 
        Compute an inverse for x modulo p, assuming that x
        is not divisible by p.
     */
    mpz_t v;
    mpz_init(v);
    mpz_mod(v,x,this->p);
    bool is_not_inverse = (mpz_cmp_ui(v, 0) == 0);
    

    if (is_not_inverse)
    {
        gmp_printf("x = %Zd and p =  %Zd is impossible inverse  \n", x, p);
        assert(!is_not_inverse);
    }

    mpz_sub_ui(v, p, 2);    //v = p - 2
    mpz_powm(ret, x, v, p); // ret =  pow(x, self.p-2, self.p)
    
    
    mpz_clear(v);  
}

void EllipticCurve::ec_point_inv(Point &ret,const Point &point)
{
    if(mpz_cmp_ui(point.x,0) == 0 && mpz_cmp_ui(point.y,0) == 0){
        mpz_set(ret.x,point.x);
        mpz_set(ret.y,point.y);
        return ;
    }

    // Point(P.x, (-P.y) % self.p)
    mpz_t v;
    mpz_init(v);
    mpz_neg(v, point.y);
    mpz_mod(v, v, this->p);
    mpz_set(ret.x,point.x);
    mpz_set(ret.y,v);
    mpz_clear(v);  
}


void EllipticCurve::ec_point_double_add(Point &ret,const Point &point){
    mpz_t numer, denom, lambda;
    //mpz_inits(numer, denom, lambda, NULL);
    mpz_init(numer),mpz_init(denom),mpz_init(lambda);
    
    // calculate lambda
    mpz_mul(numer, point.x, point.x);
    mpz_mul_ui(numer, numer, 3);
    mpz_add(numer, numer, this->a);
    mpz_mul_ui(denom, point.y, 2);
    mpz_invert(denom, denom, this->p);  // TODO check status
    mpz_mul(lambda, numer, denom);
    mpz_mod(lambda, lambda, this->p);
    
    // calculate resulting x coord
    mpz_mul(ret.x, lambda, lambda);
    mpz_sub(ret.x, ret.x, point.x);
    mpz_sub(ret.x, ret.x, point.x);
    mpz_mod(ret.x, ret.x, this->p);

    //calculate resulting y coord
    mpz_sub(ret.y, point.x, ret.x);
    mpz_mul(ret.y, lambda, ret.y);
    mpz_sub(ret.y, ret.y, point.y);
    mpz_mod(ret.y, ret.y, this->p);
    
    //mpz_clears(numer, denom, lambda, NULL);
    mpz_clear(numer),mpz_clear(denom),mpz_clear(lambda);
}

void EllipticCurve::ec_point_add(Point &rop,const Point &op1,const Point &op2){
    mpz_t xdiff, ydiff, lambda;
    mpz_init(xdiff), mpz_init(ydiff), mpz_init(lambda);

    // calculate lambda
    mpz_sub(ydiff, op2.y, op1.y);
    mpz_sub(xdiff, op2.x, op1.x);
    mpz_invert(xdiff, xdiff, this->p);  // TODO check status
    mpz_mul(lambda, ydiff, xdiff);
    mpz_mod(lambda, lambda, this->p);

    // calculate resulting x coord
    mpz_mul(rop.x, lambda, lambda);
    mpz_sub(rop.x, rop.x, op1.x);
    mpz_sub(rop.x, rop.x, op2.x);
    mpz_mod(rop.x, rop.x, this->p);

    //calculate resulting y coord
    mpz_sub(rop.y, op1.x, rop.x);
    mpz_mul(rop.y, lambda, rop.y);
    mpz_sub(rop.y, rop.y, op1.y);
    mpz_mod(rop.y, rop.y, this->p);

    
    mpz_clear(xdiff), mpz_clear(ydiff), mpz_clear(lambda);
}

		
void EllipticCurve::ec_point_scalar_mul(Point & rop, const mpz_t& scalar,const Point &point) {
    Point R0, R1, tmp;
    
    mpz_set(R0.x, point.x);
    mpz_set(R0.y, point.y);
    
    this->ec_point_double_add(R1,point);
    
    int dbits = mpz_sizeinbase(scalar, 2);
    int  i;
    
    for(i = dbits - 2; i >= 0; i--) {
        if(mpz_tstbit(scalar, i)) {
            mpz_set(tmp.x, R0.x);
            mpz_set(tmp.y, R0.y);
            this->ec_point_add(R0, R1, tmp);
            mpz_set(tmp.x, R1.x);
            mpz_set(tmp.y, R1.y);
            this->ec_point_double_add(R1,tmp);
        }
        else {
            mpz_set(tmp.x, R1.x);
            mpz_set(tmp.y, R1.y);
            this->ec_point_add(R1, R0, tmp);
            mpz_set(tmp.x, R0.x);
            mpz_set(tmp.y, R0.y);
            this->ec_point_double_add(R0,tmp);
        }
    }
    
    mpz_set(rop.x, R0.x);
    mpz_set(rop.y, R0.y);
}