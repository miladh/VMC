#include <unittest++/UnitTest++.h>
#include <src/Jastrow/nojastrow.h>

TEST(JastrowTest) {
    NoJastrow jas(2);
    mat A = zeros(2,3);
    double result = jas.evaluateJastrow(A);
    CHECK(result == 1);
}
