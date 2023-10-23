// loading the library while using root
{
    string cana_dir = string(getenv("PWD")) + "/../";

    // load library
    gSystem->Load((cana_dir + "/libCAna.so").c_str());
    // gSystem->Load("libgfortran.so.3");
    // load include folder
    gSystem->AddIncludePath((cana_dir + "/include").c_str());
    // add include folder for the interpreter
    gInterpreter->AddIncludePath((cana_dir + "/include").c_str());
}
