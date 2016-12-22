#ifndef EC_HPP
#define EC_HPP

#include <gcrypt.h>
#include <iostream>

// Набор параметров точки и кривой (HEX)
struct
{
    const char *x = "0091E38443A5E82C0D880923425712B2BB658B9196932E02C78B2582FE742DAA28";
    const char *y = "0032879423AB1A0375895786C4BB46E9565FDE0B5344766740AF268ADB32322E5C";
} point_param;


const unsigned char pstr[]={0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFD,0x97};
const unsigned char astr[]={0xFF,0xFF,0xFA,0x42,0x6C};
const unsigned char bstr[]={0x76,0x43,0xF2,0x6A};

struct
{
    const char *p = "00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97";
    const char *a = "fffffa426c";
    const char *b  = "7643f26a";
const char *q = "00400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67";


    /*
    const char *p = "00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97";
    const char *a = "00C2173F1513981673AF4892C23035A27CE25E2013BF95AA33B22C656F277E7335";
    const char *b = "00295F9BAE7428ED9CCC20E7C359A9D41A22FCCD9108E17BF7BA9337A6F8AE9513";
    const char *q = "00400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67";
    /*
    const char *p = "7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED";
    const char *a = "000000000000000000000000000000000000000000000000000000000000076d06";
    const char *b = "000000000000000000000000000000000000000000000000000000000000000001";
    const char *q = "00400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67";
    /*
    const char *p = "00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97";
    const char *a = "00C2173F1513981673AF4892C23035A27CE25E2013BF95AA33B22C656F277E7335";
    const char *b = "00295F9BAE7428ED9CCC20E7C359A9D41A22FCCD9108E17BF7BA9337A6F8AE9513";
    const char *q = "00400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67";*/
} ec_param;


class Point
{
public:
    gcry_mpi_t x;
    gcry_mpi_t y;
    gcry_mpi_t z;

    Point();
    ~Point();

    // Вывести координаты точки
    void print();
};

// Вейерштрасс: y^2 = x^3 + ax + b (p)
class EllipticCurve
{
public:
    gcry_mpi_t p;
    gcry_mpi_t a;
    gcry_mpi_t b;
    gcry_mpi_t q;
    gcry_mpi_t k;
    gcry_mpi_t l;

    Point P;
    Point Q;

    EllipticCurve();
    ~EllipticCurve();

    //Нахождение y^2 подстановкой x в уравнение кривой by^2=x^3+a*x^2+x (p)
    gcry_mpi_t GetF(gcry_mpi_t x);

    // Критерий Эйлера на то, f(x) - квадратичный вычет в p или нет?
    int euler_criteria(gcry_mpi_t fx0);

    // Вычисление y = f(x)^((p+1)/4) (p) при условии, что p = 3 (4)
    gcry_mpi_t compY(gcry_mpi_t x0);

    // Проверка принадлежности аффинной точки к кривой
    int afBelong( Point &point);

    // Проверка принадлежности проективной точки к кривой
    int prBelong(const Point &point);

    // Построение точки на кривой
    int build_point(Point &pk, int mode);

    // Операция удвоения точки в проективных координатах
    void doubling_point(Point &dp, const Point &p1);

    void doubling_pointP(Point &dp, const Point &p1);

    // Операция сложения двух точек в проективных координатах
    void addPointsAff(Point &p3, const Point &p1, const Point &p2);

    // Вычисление кратной точки Q = k*P
    void comp_mult_point(Point &k_point, const Point &point, const gcry_mpi_t m);

    // Проверка правильности работы программы
    bool extra_check();
};

void show_mpi(gcry_mpi_t mpi);

#endif // EC_HPP
