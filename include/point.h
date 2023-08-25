#ifndef POINT_H
#define POINT_H

#include<gmpxx.h>

// point in a prime field
class Point{
    public:
        mpz_t x, y;
        Point()
        {
            mpz_init(x);
            mpz_init(y);
        }
        Point(const char * x, const char * y, int base)
        {
            mpz_init_set_str(this->x, x, base);
            mpz_init_set_str(this->y, y, base);
        }
        Point(const mpz_t& x, const mpz_t& y)
        {
            mpz_init_set(this->x, x);
            mpz_init_set(this->y, y);
            
        }
        Point(const Point& p):Point(p.x,p.y){
        }

        Point& operator=(const Point& p){
            if (this == &p) {
                return *this; // Self-assignment check
            }
            mpz_init_set(this->x, p.x);
            mpz_init_set(this->y, p.y);

            return *this;
        }


        ~Point(){
             mpz_clear(this->x);
             mpz_clear(this->y);
        }

        void set_x(const char * x, int base)
        {
            mpz_set_str(this->x, x, base);
        }
        void set_x(const mpz_t& x)
        {
            mpz_set(this->x, x);
        }


        void set_y(const char * y, int base)
        {
            mpz_set_str(this->y, y, base);
        }
        void set_y(const mpz_t& y)
        {
            mpz_set(this->y, y);
        }
};



#endif