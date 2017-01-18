#ifndef TEST
#define TEST
#include <ec.hpp>
#include <iostream>

using namespace std;

int testBuilding( MontgometyCurveParameters ec_param){
    EllipticCurve curve(ec_param);
    Point Q1, Q2;
    curve.buildPoint(Q1,1);
    std::cout << "Q1\n";
    if(curve.belongingAffine(Q1) == 0) cout << "belong\n";
    Q1.print();
    curve.buildPoint(Q2,0);
    std::cout << "Q2\n";
    if(curve.belongingAffine(Q2) == 0) cout << "belong\n";
    Q2.print();
    return 0;
}
int testDoubling(MontgometyCurveParameters ec_param){
    EllipticCurve curve(ec_param);
    Point Q1, Q2;

    curve.buildPoint(Q1,0);

    curve.buildPoint(Q2,0);

    curve.doublePoint(Q2, Q1);
    std::cout << "go\n";
    if(curve.belongingAffine(Q2) == 0){
        cout << "belong\n";
    }
    else{
        cout << "non belong\n";
    }

    curve.doublePoint(Q1, Q2);
    if(curve.belongingAffine(Q1) == 0){
        cout << "belong\n";
    }
    else{
        cout << "non belong\n";
    }
    curve.doublePoint(Q2, Q1);
    if(curve.belongingAffine(Q2) == 0){
        cout << "belong\n";
    }
    else{
        cout << "non belong\n";
    }
    curve.doublePoint(Q2, Q2);
    if(curve.belongingAffine(Q2) == 0){
        cout << "belong\n";
    }
    else{
        cout << "non belong\n";
    }
    curve.doublePoint(Q2, Q2);
    if(curve.belongingAffine(Q2) == 0){
        cout << "belong\n";
    }
    else{
        cout << "non belong\n";
    }
    return 0;
}
int testSum(MontgometyCurveParameters ec_param){
    EllipticCurve curve(ec_param);


    Point Q1, Q2;

    curve.buildPoint(Q1,0);
    curve.buildPoint(Q2,0);
    Point DP;
    curve.doublePoint(DP,Q1);
    if (curve.belongingAffine(DP) != 0) return 1;
    cout << endl;

    cout << "Тестирование операции сложения точек:\n";
    Point SumP;

    Q1.print();
    std::cout << "\n";
    DP.print();
    std::cout << "\n";
    curve.addPointsAff(SumP, Q1, DP);
    if (curve.belongingAffine(SumP) == 0)
    {
        cout << "Сложение - выполнено.\n";
    }
    else
    {
        cout << "Сложение - НЕ выполнено.\n";
        return 1;
    }
    SumP.print();
    cout << endl;
    return 0;
}


int test(MontgometyCurveParameters ec_param){
    EllipticCurve curve(ec_param);
    Point Q1, Q2;

    std::cout <<"Q1 \n";
    curve.buildPoint(Q1,1);
    Q1.print();
    std::cout <<"Q2 \n";
    curve.buildPoint(Q2,0);
    cout << "Тестирование операции удвоения:\n";
    Point DP;
    curve.doublePoint(DP,Q1);// curve.P);
    if (curve.belongingAffine(DP) == 0)
    {
        cout << "Удвоенная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Удвоенная точка не принадлежит кривой.\n";
    }
    DP.print();
    cout << endl;
    cout << endl;

    cout << "Тестирование операции сложения точек:\n";
    Point SumP;
    curve.addPointsAff(SumP, Q1, Q2);
    if (curve.belongingAffine(SumP) == 0)
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

    cout << "Тестирование операции нахождения кратной точки:\n";
    Point KP1;
    curve.calculateMultiplePoint(KP1, Q2, curve.l);
    if (curve.belongingAffine(KP1) == 0)
    {
        cout << "Кратная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Кратная точка не принадлежит кривой.\n";
    }
    KP1.print();
    cout << endl;
    return 0;
}
#endif // TEST

