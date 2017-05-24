#ifndef SPARPY_H_
#define SPARPY_H_

#include "Aboria.h"
#include <boost/python.hpp>

namespace sparpy {

using namespace Aboria;
using namespace boost::python;

ABORIA_VARIABLE(scalar, double, "an example scalar")
ABORIA_VARIABLE_VECTOR(force_d,double,"force")
template <unsigned int D>
using ParticlesType = Particles<std::tuple<scalar,force_d<D>>,D>;


}

#include "interactions.hpp"
#include "simulation.hpp"

#endif
