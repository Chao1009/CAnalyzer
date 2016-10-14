// loading the library while using root
{
    // load library
    gSystem->Load("/home/chao/saGDH/analysis/tools/libCAna.so");
    // load include folder
    gSystem->AddIncludePath("/home/chao/saGDH/analysis/tools/include");
    // add include folder for the interpreter
    gInterpreter->AddIncludePath("/home/chao/saGDH/analysis/tools/include");
}
