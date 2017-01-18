#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>
#include <test.hpp>

using namespace std;

int main()
{
    gcry_control (GCRYCTL_DISABLE_SECMEM, 0);
    gcry_control (GCRYCTL_INITIALIZATION_FINISHED, 0);
    MontgometyCurveParameters param;

    const char *p = "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffd97";
    const char *a = "5432954bd3087427c0da951d18f761cf0ce9c8f792b1be02b1306ad0d4938f11";
    const char *b = "5432954bd3087427c0da951d18f761cf0ce9c8f792b1be02b1306ad0d4938f0f";
    const char *q = "400000000000000000000000000000000fd8cddfc87b6635c115af556c360c66";
    param.a = a;
    param.b = b;
    param.p = p;
    param.q = q;

    test(param);
    return 0;
}

