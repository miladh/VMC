#include "jastrow.h"

Jastrow::Jastrow(const uint nParticles):
    nParticles(nParticles),
    a(zeros(nParticles,nParticles))
{
}
