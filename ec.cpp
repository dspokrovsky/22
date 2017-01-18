#include <ec.hpp>
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <gcrypt.h>
#include <gmp.h>

#define MAX_MPI_BUF 256

using namespace std;

void show_mpi(gcry_mpi_t mpi)
{
    unsigned char *check_buf = (unsigned char *)malloc(MAX_MPI_BUF);
    memset(check_buf, '\0', MAX_MPI_BUF);
    gcry_mpi_print(GCRYMPI_FMT_HEX, check_buf, MAX_MPI_BUF, NULL, mpi);
    printf("%s\n", check_buf);
}

// -----------------------------------------------------
// Методы точки

Point::Point()
{
    x = gcry_mpi_new(0);
    y = gcry_mpi_new(0);
    z = gcry_mpi_new(0);
}

Point::~Point()
{
    gcry_mpi_release(x);
    gcry_mpi_release(y);
    gcry_mpi_release(z);
}

void Point::print()
{
    cout << "x = ";
    show_mpi(x);
    cout << "y = ";
    show_mpi(y);
    cout << "z = ";
    show_mpi(z);
}

// -----------------------------------------------------
// Методы эллиптической кривой



EllipticCurve::EllipticCurve(MontgometyCurveParameters ec_param)
{
    /*cout << "Эллиптическая кривая" << endl;
    cout << "b*y^2 = x^3 + ax^2 + x (mod p)" << endl;
    cout << "Параметры кривой:\n" << endl;*/
    p = gcry_mpi_new(256);
    a = gcry_mpi_new(256);
    b = gcry_mpi_new(256);
    q = gcry_mpi_new(256);
    if (gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, ec_param.p, 0, NULL) != 0)
        cout << "Ошибка записи p.\n";
    if (gcry_mpi_scan(&a, GCRYMPI_FMT_HEX, ec_param.a, 0, NULL) != 0)
        cout << "Ошибка записи a.\n";
    if (gcry_mpi_scan(&b, GCRYMPI_FMT_HEX, ec_param.b, 0, NULL) != 0)
        cout << "Ошибка записи b.\n";
    if (gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, ec_param.q, 0, NULL) != 0)
        cout << "Ошибка записи q.\n";
    //---
    cout << "p = ";
    show_mpi(p);
    cout << "a = ";
    show_mpi(a);
    cout << "b = ";
    show_mpi(b);
    cout << "q = ";
    show_mpi(q);

    // Генерация случайного числа k
    k = gcry_mpi_new(0);
    gcry_mpi_randomize(k, 50, GCRY_WEAK_RANDOM);
    cout << "k = ";
    show_mpi(k);

    // Генерация l
    l = gcry_mpi_new(0);
    gcry_mpi_set_ui(l, 20);
    cout << "l = ";
    show_mpi(l);
    cout << endl;
}

EllipticCurve::~EllipticCurve()
{
    gcry_mpi_release(p);
    gcry_mpi_release(a);
    gcry_mpi_release(b);
    gcry_mpi_release(q);
    gcry_mpi_release(k);
}

gcry_mpi_t EllipticCurve::calculateRightPart(gcry_mpi_t x)
{

    gcry_mpi_t fx = gcry_mpi_new(0);
    gcry_mpi_t buf = gcry_mpi_new(0);
    gcry_mpi_set_ui(buf, 2);
    gcry_mpi_powm(fx, x, buf, p);
    gcry_mpi_mulm(fx, a, fx, p);
    gcry_mpi_set_ui(buf, 3);
    gcry_mpi_powm(buf, x, buf, p);
    gcry_mpi_addm(fx, buf, fx, p);
    gcry_mpi_addm(fx, fx, x, p);
    gcry_mpi_release(buf);
    return fx;
}

int EllipticCurve::eulerCriteria(gcry_mpi_t fx0)
{
    gcry_mpi_t exp = gcry_mpi_new(0);
    gcry_mpi_t div = gcry_mpi_new(0);
    gcry_mpi_t X = gcry_mpi_new(0);
    gcry_mpi_sub_ui(exp, p, 1);
    gcry_mpi_set_ui(div, 2);
    gcry_mpi_div(exp, NULL, exp, div, 0);
    gcry_mpi_powm(X, fx0, exp, p);
    gcry_mpi_release(exp);
    gcry_mpi_release(div);
    int cmp = gcry_mpi_cmp_ui(X,1);
    gcry_mpi_release(X);
    return cmp;
}

gcry_mpi_t EllipticCurve::calculateY(gcry_mpi_t x0)
{
    gcry_mpi_t fx0 = calculateRightPart(x0);
    gcry_mpi_t y0 = gcry_mpi_new(0);
    gcry_mpi_t invers = gcry_mpi_new(0);
    if (!gcry_mpi_invm(invers,b,p)) std::cout << "ATTENTION NO INVERSE ELEMENT compY\n";//TODO
    gcry_mpi_mulm(fx0, invers, fx0, p);
    if (eulerCriteria(fx0)==0)
    {
        cout << "f(x) - квадратичный вычет в р.\n";
        gcry_mpi_t exp = gcry_mpi_new(0);
        gcry_mpi_t div = gcry_mpi_new(0);
        gcry_mpi_add_ui(exp, p, 1);
        gcry_mpi_set_ui(div, 4);
        gcry_mpi_div(exp, NULL, exp, div, 0);
        gcry_mpi_powm(y0, fx0, exp, p);
        gcry_mpi_release(exp);
        gcry_mpi_release(div);
    }
    else
    {
        cout << "f(x) - не квадратичный вычет в р.\n";
        gcry_mpi_set_ui(y0, 0);
    }
    return y0;
}

int EllipticCurve::belongingAffine(Point &point)
{
    // Подстановка в уравнение кривой координат точки point и проверка на принадлежность
    point.print();
    gcry_mpi_t left = gcry_mpi_new(256);
    gcry_mpi_t exp = gcry_mpi_new(256);
    gcry_mpi_set_ui(exp, 2);
    gcry_mpi_powm(left, point.y, exp, p);
    gcry_mpi_mulm(left, left, b, p);
    //show_mpi(left);
    gcry_mpi_t right = calculateRightPart(point.x);
    //show_mpi(right);
    //std::cout << '\n';
    int cmp = gcry_mpi_cmp(left, right);
    gcry_mpi_release(left);
    gcry_mpi_release(right);
    return cmp;
}

int EllipticCurve::buildPoint(Point &pk,int mode)
{
    switch(mode)
    {
    case 0:
    {
        cout << "Будут использованы ГОСТ параметры для точки P(x,y,z).\n";
        if (gcry_mpi_scan(&P.x, GCRYMPI_FMT_HEX, point_param.x, 0, NULL) != 0)
        {
            cout << "Ошибка записи координаты x.\n";
            return 1;
        }
        if (gcry_mpi_scan(&P.y, GCRYMPI_FMT_HEX, point_param.y, 0, NULL) != 0)
        {
            cout << "Ошибка записи координаты x.\n";
            return 1;
        }
        gcry_mpi_set_ui(P.z, 1);
        P.print();
        cout << endl;
        break;
    }
    case 3:
    {
        cout << "Будет сгенерирована точка 9.\n";
        gcry_mpi_set_ui(P.x, 9);
        P.y = calculateY(P.x);
        gcry_mpi_set_ui(P.z, 1);
        P.print();
        cout << endl;
        break;
    }
    case 1:
    {
        //cout << "Будет сгенерирована случайная точка.\n";
        gcry_mpi_randomize(P.x, 256, GCRY_WEAK_RANDOM);
        gcry_mpi_mod(P.x,P.x,p);
        P.y = calculateY(P.x);
        gcry_mpi_set_ui(P.z, 1);
        //P.print();
        cout << endl;
        break;
    }
    default:
    {
        cout << "Неверный режим.\n";
        return 1;
        break;
    }
    }
    pk.x= P.x;
    pk.y= P.y;
    pk.z= P.z;
    if (belongingAffine(P) == 0)
    {
        cout << "--Точка P(x,y) принадлежит кривой.\n";
        cout << endl;
        return 0;
    }
    else
    {
        cout << "--Точка P(x,y) не принадлежит кривой.\n";
        cout << endl;
        return 1;
    }
}

void EllipticCurve::doublePoint(Point &dp, const Point &p1)
{
    gcry_mpi_t buf1 = gcry_mpi_new(0);
    gcry_mpi_t buf2 = gcry_mpi_new(0);
    gcry_mpi_t buf3 = gcry_mpi_new(0);
    gcry_mpi_t two = gcry_mpi_new(0);;
    gcry_mpi_set_ui(two,2);
    gcry_mpi_t three = gcry_mpi_new(0);
    gcry_mpi_set_ui(three,3);
    gcry_mpi_t one = gcry_mpi_new(0);    
    gcry_mpi_set_ui(one,1);
    //start
    gcry_mpi_mulm(buf1, two, b,p);
    gcry_mpi_mulm(buf1, buf1, p1.y,p);
    gcry_mpi_powm(buf1, buf1, two,p);
    if (!gcry_mpi_invm(buf1,buf1,p)) std::cout << "ATTENTION NO INVERSE ELEMENT dp\n";//TODO
    //buf1-save
    gcry_mpi_powm(buf2, p1.x, two,p);       //x1^2
    gcry_mpi_mulm(buf2, buf2, three,p);         //3*x1^2
    gcry_mpi_mulm(buf3, a, two,p);          //2a
    gcry_mpi_mulm(buf3, buf3, p1.x,p);     // 2ax
    gcry_mpi_addm(buf3, buf3, buf2,p);      //3*x1^2+2ax1
    gcry_mpi_addm(buf3, buf3, one,p);       //3*x1^2+2ax1+1
    gcry_mpi_powm(buf3, buf3, two,p);       // ( 3*x1^2+2ax+1)^2
    gcry_mpi_mulm(buf3, buf3, buf1,p);
    gcry_mpi_mulm(buf3, buf3, b,p);
    gcry_mpi_subm(buf3, buf3, a,p);
    gcry_mpi_subm(buf3, buf3, p1.x,p);
    gcry_mpi_subm(dp.x, buf3, p1.x,p);
    //y3 - start

    gcry_mpi_mulm(buf1, two, b,p);
    gcry_mpi_mulm(buf1, buf1, p1.y,p);
    gcry_mpi_powm(buf1, buf1, three,p);     //2by1
    if (!gcry_mpi_invm(buf1,buf1,p)) std::cout << "ATTENTION NO INVERSE ELEMENT dp2\n";//TODO
                                            //TODO
    gcry_mpi_powm(buf2, p1.x, two,p);
    gcry_mpi_mulm(buf2, buf2, three,p);     //(3)
    gcry_mpi_mulm(buf3, a, two,p);
    gcry_mpi_mulm(buf3, buf3, p1.x,p);
    gcry_mpi_addm(buf3, buf3, buf2,p);
    gcry_mpi_addm(buf3, buf3, one,p);

    gcry_mpi_powm(buf3, buf3, three,p);

    gcry_mpi_mulm(buf3, buf3, buf1,p);
    gcry_mpi_mulm(buf1, buf3, b,p);  //b1

                                        //TODO
    gcry_mpi_powm(buf2, p1.x, two,p);
    gcry_mpi_mulm(buf2, buf2, three,p);
    gcry_mpi_mulm(buf3, a, two,p);
    gcry_mpi_mulm(buf3, buf3, p1.x,p);
    gcry_mpi_addm(buf3, buf3, buf2,p);
    gcry_mpi_addm(buf3, buf3, one,p); //b3

    gcry_mpi_mulm(buf2, p1.x, two,p);
    gcry_mpi_addm(buf2, p1.x, buf2,p);
    gcry_mpi_addm(buf2, a, buf2,p);

    gcry_mpi_mulm(buf2, buf2, buf3,p);

    gcry_mpi_t bin = gcry_mpi_new(0);
    gcry_mpi_mulm(bin, b, two,p);
    gcry_mpi_mulm(bin, p1.y, bin,p);
    if (!gcry_mpi_invm(bin,bin,p)) std::cout << "ATTENTION NO INVERSE ELEMENT dp3 \n";//TODO
    gcry_mpi_mulm(buf2, buf2, bin,p);

    gcry_mpi_subm(buf2, buf2, buf1,p);
    gcry_mpi_subm(dp.y, buf2, p1.y,p);
    dp.z = p1.z;
}



void EllipticCurve::addPointsAff(Point &p3, const Point &p1, const Point &p2)
{
    gcry_mpi_t y2y1 = gcry_mpi_new(0);
    gcry_mpi_t x2x1 = gcry_mpi_new(0);
    gcry_mpi_t two = gcry_mpi_new(0);
    gcry_mpi_set_ui(two,2);
    gcry_mpi_t buf = gcry_mpi_new(0);
    gcry_mpi_subm(y2y1, p2.y, p1.y,p); //y2-y1
    gcry_mpi_powm(y2y1,y2y1,two,p);
    gcry_mpi_mulm(y2y1, y2y1, b,p); //b*(y2-y1)^2

    gcry_mpi_subm(x2x1, p2.x, p1.x,p); //()
    gcry_mpi_powm(x2x1,x2x1,two,p);
    if (!gcry_mpi_invm(x2x1,x2x1,p)) std::cout<< "ATTENTION NO INVERSE ELEMENT 1\n";//TODO
    gcry_mpi_mulm(y2y1, y2y1, x2x1,p);
    gcry_mpi_subm(y2y1, y2y1, a,p);
    gcry_mpi_subm(y2y1, y2y1, p1.x,p);
    gcry_mpi_subm(p3.x, y2y1, p2.x,p);
    //x3 done
    gcry_mpi_subm(x2x1, p2.x, p1.x,p);
    //gcry_mpi_mulm(buf, x2x1, x2x1,p);
    gcry_mpi_t three = gcry_mpi_new(0);
    gcry_mpi_set_ui(three,3);
    gcry_mpi_powm(buf,x2x1,three,p);//(x2-x1)^3
    if (!gcry_mpi_invm(x2x1,buf,p)) std::cout<<"ATTENTION NO INVERSE ELEMENT 2\n";//TODO
    gcry_mpi_subm(y2y1, p2.y, p1.y,p);
    gcry_mpi_powm(y2y1,y2y1,three,p);
    gcry_mpi_mulm(y2y1, b, y2y1,p);
    gcry_mpi_mulm(buf, x2x1, y2y1,p); // second b*(y2-y1)3/(x2-x1)3
    gcry_mpi_subm(x2x1, p2.x, p1.x,p);
    if (!gcry_mpi_invm(x2x1,x2x1,p)) std::cout <<"ATTENTION NO INVERSE ELEMENT 3\n";//TODO
    gcry_mpi_subm(y2y1, p2.y, p1.y,p);
    gcry_mpi_mulm(y2y1, x2x1, y2y1,p);
    gcry_mpi_mulm(x2x1, p1.x, two,p);
    gcry_mpi_addm(x2x1,x2x1,p2.x,p);
    gcry_mpi_addm(x2x1,x2x1,a,p);
    gcry_mpi_mulm(x2x1,x2x1,y2y1,p);
    gcry_mpi_subm(y2y1, x2x1, buf,p);
    gcry_mpi_subm(p3.y, y2y1, p1.y,p);
    p3.z=p1.z;
}

void EllipticCurve::calculateMultiplePoint(Point &k_point, const Point &point, const gcry_mpi_t m)
{
    Point temp;
    temp.x = gcry_mpi_copy(point.x);
    temp.y = gcry_mpi_copy(point.y);
    temp.z = gcry_mpi_copy(point.z);

    gcry_mpi_set_ui(k_point.x, 0);
    gcry_mpi_set_ui(k_point.y, 1);
    gcry_mpi_set_ui(k_point.z, 1);

    for(int i = (int)gcry_mpi_get_nbits(m)-1; i > -1;--i)
    {
        doublePoint(k_point, temp);
        if (gcry_mpi_test_bit(m, i))
            addPointsAff(k_point, k_point, temp);
    }
}

void EllipticCurve::doubling_pointP(Point &dp, const Point p1){
    gcry_mpi_t exp =  gcry_mpi_new(0);
    /*  x2 = (x+z)^2 * (x-z)^2 */
    /*  z2 = ((x+z)^2 - (x-z)^2)*((x+z)^2 + ((A-2)/4)((x+z)^2 - (x-z)^2)) */
    gcry_mpi_t A =  gcry_mpi_new(0);     /* A = (x+z) */
    gcry_mpi_addm(A,p1.x,p1.z,p);


    gcry_mpi_t B =  gcry_mpi_new(0);       /* B = (x-z) */
    gcry_mpi_subm(B,p1.x,p1.z,p);

    gcry_mpi_set_ui(exp,2);
                                       /* A = (x+z)^2 */
    gcry_mpi_powm(A,A,exp,p);
                                          /* B = (x-z)^2 */
    gcry_mpi_powm(B,B,exp,p);
                                        /* x2 = (x+z)^2 * (x-z)^2 */

    gcry_mpi_mulm(dp.x,A,B,p);
                                  /* B = (x+z)^2 - (x-z)^2 */

    gcry_mpi_subm(B,A,B,p);

    /* (486662-2)/4 = 121665 */
    gcry_mpi_t bufa =  gcry_mpi_new(0);
    gcry_mpi_subm(bufa,a,exp,p);
    gcry_mpi_set_ui(exp,4);
    gcry_mpi_div(bufa,0,bufa,exp,0);
    gcry_mpi_mulm(bufa,bufa,B,p);
    gcry_mpi_addm(bufa,A,B,p);
                                  /* z2 = (B)*((x+z)^2 + ((A-2)/4)(B)) */
    gcry_mpi_mulm(dp.z,bufa,B,p);
    //dp.y = p1.y;


    gcry_mpi_set_ui(exp,3);
    gcry_mpi_powm(A,dp.x,exp,p);  //x3
    gcry_mpi_set_ui(exp,2);
    gcry_mpi_powm(B,dp.x,exp,p);
    gcry_mpi_mulm(B,dp.z,B,p);
    gcry_mpi_mulm(B,a,B,p);  //ax2z
    gcry_mpi_addm(A,A,B,p);   //x3+ax2z
    gcry_mpi_powm(B,dp.z,exp,p);
    gcry_mpi_mulm(B,B,dp.x,p);    //xz2
    gcry_mpi_addm(B,B,A,p); //right


    bufa = dp.z;
    if (!gcry_mpi_invm(bufa,bufa,p)) std::cout << "ATTENTION NO INVERSE ELEMENT\n";//TODO
    gcry_mpi_mulm(B,bufa,B,p);
    if (!gcry_mpi_invm(bufa,b,p)) std::cout << "ATTENTION NO INVERSE ELEMENT\n";//TODO
    gcry_mpi_mulm(B,bufa,B,p);
    dp.y = calculateY(B);

}

int EllipticCurve::prBelong(const Point &point)
{
    // point в проективных координатах
    // Подстановка в уравнение
    // b*Y^2*Z = X^3 + a*X*Z^2 + X*Z^2 (p)

    // Левая часть
    gcry_mpi_t left = gcry_mpi_new(0);
    gcry_mpi_t exp = gcry_mpi_new(0);
    gcry_mpi_set_ui(exp, 2);
    gcry_mpi_powm(left, point.y, exp, p);
    gcry_mpi_mulm(left, left, point.z, p);
    gcry_mpi_mulm(left, left, b, p);

    // Правая часть
    gcry_mpi_t right = gcry_mpi_new(0);
    gcry_mpi_t first = gcry_mpi_new(0);
    gcry_mpi_t second = gcry_mpi_new(0);
    gcry_mpi_t third = gcry_mpi_new(0);
    gcry_mpi_powm(second, point.z, exp, p);
    gcry_mpi_mulm(second, second, point.x, p);
    gcry_mpi_mulm(second, second, a, p);

    gcry_mpi_set_ui(exp, 3);
    gcry_mpi_powm(first, point.x, exp, p);

    gcry_mpi_set_ui(exp, 2);
    gcry_mpi_powm(third, point.z, exp, p);
    gcry_mpi_mulm(third, third, point.x, p);
    gcry_mpi_addm(right, first, second, p);
    gcry_mpi_addm(right, right, third, p);

    int cmp = gcry_mpi_cmp(left, right);
    gcry_mpi_release(left);
    gcry_mpi_release(right);
    gcry_mpi_release(exp);
    gcry_mpi_release(first);
    gcry_mpi_release(second);
    gcry_mpi_release(third);
    return cmp;
}

/*
void EllipticCurve::addPointsAff(Point &p3, const Point &p1, const Point &p2)
{
    gcry_mpi_t y2y1 = gcry_mpi_new(0);
    gcry_mpi_t x2x1 = gcry_mpi_new(0);
    gcry_mpi_t buf = gcry_mpi_new(0);
    gcry_mpi_subm(y2y1, p2.y, p1.y,p);
    gcry_mpi_mulm(y2y1, y2y1, y2y1,p);
    gcry_mpi_mulm(y2y1, y2y1, b,p);

    gcry_mpi_subm(x2x1, p2.x, p1.x,p);
    gcry_mpi_mulm(x2x1, x2x1, x2x1,p);
    //gcry_mpi_t invers = gcry_mpi_new(0);
    if (!gcry_mpi_invm(x2x1,x2x1,p)) std::cout << "ATTENTION NO INVERSE ELEMENT add1\n";//TODO
    gcry_mpi_mulm(y2y1, y2y1, x2x1,p);
    gcry_mpi_subm(y2y1, y2y1, a,p);
    gcry_mpi_subm(y2y1, y2y1, p1.x,p);
    gcry_mpi_subm(p3.x, y2y1, p2.x,p);
    //x3 done
    gcry_mpi_subm(x2x1, p2.x, p1.x,p);
    gcry_mpi_mulm(buf, x2x1, x2x1,p);
    gcry_mpi_mulm(x2x1, buf, x2x1,p);
    if (!gcry_mpi_invm(x2x1,x2x1,p)) std::cout << "ATTENTION NO INVERSE ELEMENT\n";//TODO

    gcry_mpi_subm(y2y1, p2.y, p1.y,p);
    gcry_mpi_mulm(buf, y2y1, y2y1,p);
    gcry_mpi_mulm(y2y1, buf, y2y1,p);
    gcry_mpi_mulm(y2y1, b, y2y1,p);
    gcry_mpi_mulm(buf, x2x1, y2y1,p); //b*
    //gcry_mpi_subm(, y2y1, p1.y,p);//buf -

    gcry_mpi_subm(x2x1, p2.x, p1.x,p);
    if (!gcry_mpi_invm(x2x1,x2x1,p)) std::cout << "ATTENTION NO INVERSE ELEMENT\n";//TODO
    gcry_mpi_subm(y2y1, p2.y, p1.y,p);
    gcry_mpi_mulm(y2y1, x2x1, y2y1,p);

    gcry_mpi_t two;
    const char *qq = "2";
    gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, qq, 0, NULL);
    gcry_mpi_mulm(x2x1, p1.x, two,p);
    gcry_mpi_mulm(x2x1,x2x1,p2.x,p);
    gcry_mpi_mulm(x2x1,x2x1,a,p);
    gcry_mpi_mulm(y2y1, x2x1, y2y1,p);
    gcry_mpi_subm(y2y1, y2y1, buf,p);
    gcry_mpi_subm(p3.y, y2y1, p1.y,p);

}*/

/*
void EllipticCurve::comp_mult_point(Point &kp, const Point &p, const gcry_mpi_t m)
{
    //calculate P.^-1

    gcry_mpi_t invers = gcry_mpi_new(0);
    gcry_mpi_invm(invers , p.x ,this->p);
    gcry_mpi_set(kp.x ,invers);
    gcry_mpi_set(kp.y ,compY(kp.x));
    gcry_mpi_set_ui(kp.z, 1);

    Point invP;
    gcry_mpi_set(invP.x, kp.x);
    gcry_mpi_set(invP.y, kp.y);
    gcry_mpi_set(invP.z, kp.z);

    Point tempP;
    gcry_mpi_set(tempP.x, p.x);
    gcry_mpi_set(tempP.y, p.y);
    gcry_mpi_set(tempP.z, p.z);

    Point tempP2;
    gcry_mpi_set(tempP2.x, tempP.x);
    gcry_mpi_set(tempP2.y, tempP.y);
    gcry_mpi_set(tempP2.z, tempP.z);

    // 13P = 8P + 4P + P = (2^3)P + (2^2)P + (2^0)P
    int last = 0;
    for(int i = 0; i < (int)gcry_mpi_get_nbits(m); i++)
    {
        //std::cout << afBelong(kp)<<"****\n";
        //if (afBelong(kp)!=0) break;
        if (gcry_mpi_test_bit(m, i))
        {
            for(int j = last; j<i; ++j)
            {
                std::cout << afBelong(tempP)<<" DOUBLE ****\n";

                gcry_mpi_set(tempP2.x, tempP.x);
                gcry_mpi_set(tempP2.y, tempP.y);
                gcry_mpi_set(tempP2.z, tempP.z);

                doubling_point(tempP, tempP2);

                std::cout << afBelong(tempP)<<" DOUBLE2 ****\n";
            }
            last = i;
            if (!(gcry_mpi_cmp(kp.x, invP.x) && gcry_mpi_cmp(kp.y, invP.y)))
            {
                kp.x = gcry_mpi_copy(tempP.x);
                kp.y = gcry_mpi_copy(tempP.y);
                kp.z = gcry_mpi_copy(tempP.z);
            }
            else
            {
                std::cout << afBelong(tempP)<<" ELSE1 ****\n";
                addPointsAff(kp, kp, tempP);
                std::cout << afBelong(tempP)<<" ELSE2 ****\n";
            }

        }
    }
    std::cout << "end\n";
}
*/
/*
bool EllipticCurve::extra_check()
{
    Point Q, Q1, Q2;
    calculateMultiplePoint(Q1, P, k);
    calculateMultiplePoint(Q2, P, l);
    addPointsAff(Q, Q1, Q2);

    Point R;
    gcry_mpi_t m = gcry_mpi_new(0);
    gcry_mpi_addm(m, k, l, p);
    calculateMultiplePoint(R, P, m);
    gcry_mpi_release(m);

    if (Q.x == R.x && Q.y == R.y && Q.z == R.z)
        return true;
    else
        return false;
}
*/
