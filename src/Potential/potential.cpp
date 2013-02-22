#include "potential.h"

Potential::Potential()
{
}


/************************************************************
Name:               loadConfiguration
Description:        loads different variables
*/
void Potential::loadConfiguration(Config *cfg){

    charge=cfg->lookup("PotentialSettings.charge");
}
