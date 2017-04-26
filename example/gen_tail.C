#include "CElasTails.h"
#include "CExpData.h"
#include "canalib.h"
#include <string>
#include <iostream>

using namespace std;

void gen_tail()
{
    CExpData data("configs/data_sets_9deg.conf");
    CElasTails eltail("configs/elas_tail.conf");

    eltail.Initialize(data, 0);
    eltail.Generate();
    eltail.Output("output/tail_1147.out");
}

