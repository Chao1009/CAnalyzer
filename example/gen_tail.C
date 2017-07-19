#include "CElasTails.h"
#include "CExpData.h"
#include "canalib.h"
#include <string>
#include <iostream>

using namespace std;

void gen_tail(const string &data_conf = "configs/data_sets_9deg.conf")
{
    CExpData data(data_conf);
    CElasTails eltail("configs/elas_tail.conf");

    for(size_t i = 0; i < data.Size(); ++i)
    {
        eltail.Initialize(data, i);
        eltail.Generate();
        string output = "output/tail_" + to_string((int)data.GetSet(i).energy) + "_perp.out";
        eltail.Output(output);
    }
}

