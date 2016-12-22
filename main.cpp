#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

#define MODE_GOST 0
#define MODE_RANDOM 1

using namespace std;

int main()
{
    gcry_control (GCRYCTL_DISABLE_SECMEM, 0);
    gcry_control (GCRYCTL_INITIALIZATION_FINISHED, 0);
    EllipticCurve curve;
    Point Q1, Q2;

    curve.build_point(Q1,1);
/*
    curve.doubling_point(Q1, Q1);
    if(curve.afBelong(Q1) == 0) cout << "++++++++++";
    Q1.print();

    cout << "Тестирование операции удвоения:\n";
    Point DP;
    curve.doubling_point(DP,Q1);// curve.P);
    if (curve.afBelong(DP) == 0) //check_affine_point_belongs
    {
        cout << "Удвоенная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Удвоенная точка не принадлежит кривой.\n";
        return 1;
    }
    DP.print();
    cout << endl;
    //cout << "Тестирование вычисления у:\n";
    //curve.comp_y0(curve.P.x);
*/
    cout << endl;


    cout << "Тестирование операции сложения точек:\n";
    Point SumP;
    curve.addPointsAff(SumP, Q1, Q1);
    if (curve.afBelong(SumP) == 0)
    {
        cout << "Полученная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Полученная точка не принадлежит кривой.\n";
        return 1;
    }
    SumP.print();
    cout << endl;
/*
    cout << "Тестирование операции нахождения кратной точки 1:\n";
    Point KP1;
    curve.comp_mult_point(KP1, curve.P, curve.k);
    if (curve.check_projective_point_belongs(KP1) == 0)
    {
        cout << "Кратная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Кратная точка не принадлежит кривой.\n";
    }
    KP1.print();
    cout << endl;

    cout << "Тестирование операции нахождения кратной точки 2:\n";
    Point KP2;
    curve.comp_mult_point(KP2, curve.P, curve.l);
    if (curve.check_projective_point_belongs(KP1) == 0)
    {
        cout << "Кратная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Кратная точка не принадлежит кривой.\n";
    }
    KP2.print();
    cout << endl;

    cout << "Extra проверка правильности всего происходящего: ";
    if (curve.extra_check())
        cout << "[OK]\n";
    else
        cout << "[NOT OK]\n";*/
    return 0;
}

