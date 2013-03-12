#include "jastrow.h"

Jastrow::Jastrow(const uint &nParticles):
    nParticles(nParticles),
    a(zeros(nParticles,nParticles)),
    rOld(zeros(nParticles,3)),
    rNew(zeros(nParticles,3))
{
}
