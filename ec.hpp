#ifndef EC_HPP
#define EC_HPP

#include <gcrypt.h>
#include <iostream>

// Набор параметров точки и кривой (HEX)
struct
{
    const char *x = "0034470a1427fb7946dc1404e81ab9be8cebb35384fbb88e84f15721dfaba59374";
    const char *y = "00c8f1c5b2c7ffa6de10edd8fe2970bfe37497a3f6898459bb7506b3d620e55842";
} point_param;

struct MontgometyCurveParameters
{
    const char *p = "00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97";
    const char *a = "00C2173F1513981673AF4892C23035A27CE25E2013BF95AA33B22C656F277E7335";
    const char *b = "00295F9BAE7428ED9CCC20E7C359A9D41A22FCCD9108E17BF7BA9337A6F8AE9513";
    const char *q = "00400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67";
};// ec_param2;

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

    EllipticCurve(MontgometyCurveParameters ec_param);
    ~EllipticCurve();

    //Нахождение y^2 подстановкой x в уравнение кривой by^2=x^3+a*x^2+x (p)
    gcry_mpi_t calculateRightPart(gcry_mpi_t x);

    // Критерий Эйлера на то, f(x) - квадратичный вычет в p или нет?
    int eulerCriteria(gcry_mpi_t fx0);

    // Вычисление y = f(x)^((p+1)/4) (p) при условии, что p = 3 (4)
    gcry_mpi_t calculateY(gcry_mpi_t x0);

    // Проверка принадлежности аффинной точки к кривой
    int belongingAffine( Point &point);

    // Проверка принадлежности проективной точки к кривой
    int prBelong(const Point &point);

    // Построение точки на кривой
    int buildPoint(Point &pk, int mode);

    // Операция удвоения точки в проективных координатах
    void doublePoint(Point &dp, const Point &p1);

    void doubling_pointP(Point &dp, const Point p1);

    // Операция сложения двух точек в проективных координатах
    void addPointsAff(Point &p3, const Point &p1, const Point &p2);

    // Вычисление кратной точки Q = k*P
    //void comp_mult_point(Point &k_point, const Point &point, const gcry_mpi_t m);
    void calculateMultiplePoint(Point &kp, const Point &p, const gcry_mpi_t m);
};

void show_mpi(gcry_mpi_t mpi);

#endif // EC_HPP
