#ifndef INTERACTIONS_H_
#define INTERACTIONS_H_

#include "sparpy.h"

namespace sparpy {


template <unsigned int D>
struct exponential_force {
    typedef ParticlesType<D> particles_type;
    typedef Vector<double,D> double_d;
    typedef typename particles_type::position position;
    typedef force_d<D> force;
    typedef std::shared_ptr<ParticlesType<D>> particles_pointer;

    double m_cutoff;
    double m_epsilon;
    exponential_force(const double cutoff, const double epsilon):
        m_cutoff(cutoff),m_epsilon(epsilon)
    {}

    void operator()(particles_pointer particles1, particles_pointer particles2) {
        Symbol<position> p;
        Symbol<force> f;
        Symbol<id> id_;
        Label<0,particles_type> a(*particles1);
        Label<1,particles_type> b(*particles2);
        auto dx = create_dx(a,b);
        AccumulateWithinDistance<std::plus<double_d> > sum(m_cutoff);

        f[a] += sum(b,if_else(norm(dx)!=0,-(1.0/m_epsilon)*exp(-norm(dx)/m_epsilon)/norm(dx),0)*dx);
    }
};

template <unsigned int D>
struct lennard_jones_force {
    typedef ParticlesType<D> particles_type;
    typedef Vector<double,D> double_d;
    typedef typename particles_type::position position;
    typedef force_d<D> force;
    typedef std::shared_ptr<ParticlesType<D>> particles_pointer;

    double m_cutoff;
    double m_epsilon;
    lennard_jones_force(const double cutoff, const double epsilon):
        m_cutoff(cutoff),m_epsilon(epsilon)
    {}

    void operator()(particles_pointer particles1, particles_pointer particles2) {
        Symbol<position> p;
        Symbol<force> f;
        Symbol<id> id_;
        Label<0,particles_type> a(*particles1);
        Label<1,particles_type> b(*particles2);
        auto dx = create_dx(a,b);
        AccumulateWithinDistance<std::plus<double_d> > sum(m_cutoff);

        f[a] += sum(b,if_else(norm(dx)!=0, (1.0/dot(dx,dx))*(
                                    12*pow(m_epsilon/norm(dx),12) 
                                  - 6 *pow(m_epsilon/norm(dx), 6)
                                                           ),0)*dx);
    }
};

template <unsigned int D>
struct hard_sphere {
    typedef ParticlesType<D> particles_type;
    typedef Vector<double,D> double_d;
    typedef typename particles_type::position position;
    typedef std::shared_ptr<ParticlesType<D>> particles_pointer;
    typedef force_d<D> force;

    double m_diameter;
    hard_sphere(const double radius):
        m_diameter(2*radius)
    {}

    void operator()(particles_pointer particles1, particles_pointer particles2) {
        Symbol<position> p;
        Symbol<force> f;
        Symbol<id> id_;
        Label<0,particles_type> a(*particles1);
        Label<1,particles_type> b(*particles2);
        auto dx = create_dx(a,b);
        AccumulateWithinDistance<std::plus<double_d> > sum(m_diameter);

        p[a] += sum(b, if_else(id_[a]!=id_[b],(m_diameter/norm(dx)-1),0)*dx);
    }
};


}

#endif
