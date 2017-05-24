#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "sparpy.h"
#include "interactions.hpp"
#include "timestepping.hpp"

namespace sparpy {

template <unsigned int D>
class Simulation {
    typedef ParticlesType<D> particles_type;
    typedef typename ParticlesType<D>::position position;
    typedef force_d<D> force;
    typedef Vector<double,D> double_d;
    typedef Vector<bool,D> bool_d;
    typedef std::shared_ptr<ParticlesType<D>> particles_pointer;
    typedef std::map<particles_pointer,double> particles_storage_type;
    typedef std::vector<std::function<void()>> actions_storage_type;
    typedef std::vector<std::function<void()>> forces_storage_type;

    particles_storage_type particle_sets;
    actions_storage_type actions;
    forces_storage_type forces;
    bool m_domain_has_been_set;
    double_d m_min;
    double_d m_max;
    double_d m_min_reflect;
    double_d m_max_reflect;
    bool_d m_periodic;


public:

    template <typename F>
    void add_force(particles_pointer particles1, particles_pointer particles2, 
            const F& calc_force) {
        forces.push_back(std::bind(calc_force,particles1,particles2));
    }

    void add_action(particles_pointer particles1, particles_pointer particles2, 
            const std::function<void(particles_pointer,particles_pointer)>& calc_action) {
        actions.push_back(std::bind(calc_action,particles1,particles2));
    }


    void set_domain(const double_d& min, const double_d& max, const bool_d periodic) {
        m_min = min;
        m_max = max;
        m_min_reflect = min;
        m_max_reflect = max;
        m_periodic = periodic;
        for (int i = 0; i < D; ++i) {
            if (!m_periodic[i]) {
                m_min[i] -= m_max[i]-m_min[i];
                m_max[i] += m_max[i]-m_min[i];
            }
        }
        m_domain_has_been_set = true;
        for (auto& particle_set: particle_sets) {
            auto& particles = particle_set.first;
            particles->init_neighbour_search(min,max,periodic);
        }
    }

    void add_particles(particles_pointer particles, const double diffusion_constant) {
        if (m_domain_has_been_set) {
            particles->init_neighbour_search(m_min,m_max,m_periodic);
        }
        typename particles_storage_type::iterator search = particle_sets.find(particles);
        if (search != particle_sets.end()) {
            search->second = diffusion_constant;
        } else {
            particle_sets.insert(std::make_pair(particles,diffusion_constant));
        }
    }

    void euler_integration(const double dt, 
                           particles_pointer particles, 
                           const double diffusion_constant) {
        /*
        Symbol<position> p;
        Symbol<force> f;
        Label<0,particles_type> a(*particles);
        VectorSymbolic<double,D> vector;
        Normal N;

        p[a] += std::sqrt(2*diffusion_constant*dt)*vector(N[a],N[a],N[a]) + dt*f[a]/(1+norm(f[a])*dt);
        */

        std::normal_distribution<double> N;

        for (typename particles_type::reference i: *particles) {
            double_d& p = get<position>(i);
            double_d& f = get<force>(i);
            auto& g = get<Aboria::random>(i);
            const double scale = 1.0/(1+f.norm()*dt);
            const double diffusion = std::sqrt(2*diffusion_constant*dt);
            for (int i = 0; i < D; ++i) {
                p[i] += diffusion*N(g) + dt*f[i]*scale;
            }
        }
        particles->update_positions();
    }

    void reflective_boundaries(particles_pointer particles) {
        for (typename particles_type::reference i: *particles) {
            double_d& p = get<position>(i);
            for (int i = 0; i < D; ++i) {
                if (!m_periodic[i]) {
                    if (p[i] > m_max_reflect[i]) {
                        p[i] = 2*m_max_reflect[i]-p[i];
                    }
                    if (p[i] < m_min_reflect[i]) {
                        p[i] = 2*m_min_reflect[i]-p[i];
                    }
                }
                
            }
        }
        particles->update_positions();
    }


    void time_step(const double dt) {
        // zero forces
        for (auto& particle_set: particle_sets) {
            Symbol<force> f;
            Label<0,particles_type> a(*particle_set.first);
            f[a] = 0.0;
        }

        // calculate forces
        for (auto& calc_force: forces) {
            calc_force();
        }

        // integrate
        for (auto& particle_set: particle_sets) {
            euler_integration(dt,particle_set.first,particle_set.second);
        }

        // calculate actions
        for (auto& calc_action: actions) {
            calc_action();
        }

        // deal with reflective boundaries if needed
        for (auto& particle_set: particle_sets) {
            reflective_boundaries(particle_set.first);
        }
    }

    void integrate(const double for_time, const double dt) {
        const int timesteps = std::floor(for_time/dt);
        const double remainder_dt = for_time-timesteps*dt;

        for (int i = 0; i < timesteps; ++i) {
            time_step(dt);
        }
        time_step(remainder_dt);
    }

};

}

#endif
