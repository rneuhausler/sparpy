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
        Symbol<species> w;
        Symbol<id> id_;
        Label<0,particles_type> a(*particles1);
        Label<1,particles_type> b(*particles2);
        auto dx = create_dx(a,b);
        AccumulateWithinDistance<std::plus<double_d> > sum(m_cutoff);
        //std::cout << "calc force" << m_epsilon << m_cutoff<< std::endl;
        //f[a] += sum(b,if_else(norm(dx)!=0 && w[b]==0,(1.0/m_epsilon)*exp(-norm(dx)/m_epsilon)/norm(dx),0)*dx);
        //f[a] += sum(b,if_else(norm(dx)!=0,(1.0/m_epsilon)*exp(-norm(dx)/m_epsilon)/norm(dx),0)*dx);
	for (int i=0; i<particles1->size(); ++i) {
		for (int j=0; j<particles2->size(); ++j) {
			if (get<species>(*particles2)[j] != 0) continue;
			const double_d dx = get<position>(*particles2)[j]-get<position>(*particles1)[i];
			const double r2 = dx.squaredNorm();
			if (r2 != 0) {
				const double r = std::sqrt(r2);
				get<force>(*particles1)[i] += (1.0/m_epsilon)*std::exp(-r/m_epsilon)*dx/r; 
			}
		}
	}
    }
};



template <unsigned int D>
struct morse_force {
    typedef ParticlesType<D> particles_type;
    typedef Vector<double,D> double_d;
    typedef typename particles_type::position position;
    typedef force_d<D> force;
    typedef std::shared_ptr<ParticlesType<D>> particles_pointer;
  
    double m_cutoff;
    double m_Ca;
    double m_la;
    double m_Cr;
    double m_lr;
    double m_type;
    morse_force(const double cutoff, const double Ca, const double la, const double Cr, const double lr, const double type):
        m_cutoff(cutoff),m_Ca(Ca),m_la(la),m_Cr(Cr),m_lr(lr),m_type(type)
    {}
  
    void operator()(particles_pointer particles1, particles_pointer particles2) {
        Symbol<position> p;
        Symbol<force> f;
        Symbol<id> id_;
        Label<0,particles_type> a(*particles1);
        Label<1,particles_type> b(*particles2);
        auto dx = create_dx(a,b);
        AccumulateWithinDistance<std::plus<double_d> > sum(m_cutoff);
        //f[a] += sum(b, if_else(norm(dx)!=0, (m_Ca/m_la*exp(-norm(dx)/m_la) - m_Cr/m_lr*exp(-norm(dx)/m_lr))/norm(dx),0)*dx);
  for (int i=0; i<particles1->size(); ++i) {
    for (int j=0; j<particles2->size(); ++j) {
      if (get<species>(*particles2)[j] != m_type) continue;
      const double_d dx = get<position>(*particles2)[j]-get<position>(*particles1)[i];
      const double r2 = dx.squaredNorm();
      if (r2 != 0) {
        const double r = std::sqrt(r2);
        get<force>(*particles1)[i] += (m_Ca/m_la*std::exp(-r/m_la) - m_Cr/m_lr*std::exp(-r/m_lr))*dx/r;
      }                               
    }
  }
}
};
  
  
  
  
  

template <unsigned int D>
struct yukawa_force {
    typedef ParticlesType<D> particles_type;
    typedef Vector<double,D> double_d;
    typedef typename particles_type::position position;
    typedef force_d<D> force;
    typedef std::shared_ptr<ParticlesType<D>> particles_pointer;

    double m_cutoff;
    double m_epsilon;
    yukawa_force(const double cutoff, const double epsilon):
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

        f[a] += sum(b, if_else(norm(dx)!=0, (exp(-norm(dx)/m_epsilon)*(m_epsilon+norm(dx))/pow(norm(dx),2)) /norm(dx),0)* dx);
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

template <unsigned int D>
struct calculate_density {
  typedef ParticlesType<D> particles_type;
  typedef Vector<double,D> double_d;
  typedef typename particles_type::position position;
  typedef std::shared_ptr<ParticlesType<D>> particles_pointer;
  typedef force_d<D> force;
  typedef typename particles_type::reference reference;
  typedef typename particles_type::const_reference const_reference;
  
  double m_radius;
  double m_dt;
  calculate_density(const double radius, const double dt):
    m_radius(radius),m_dt(dt)
  {}
  
  void operator()(particles_pointer particles1, particles_pointer particles2) {
    for (reference i: *particles1) {
      for (const auto& tpl: euclidean_search(particles2->get_query(),get<position>(i),m_radius)) {
        const_reference j = std::get<0>(tpl);
        const double_d& dx = std::get<1>(tpl);
        const double r2 = dx.squaredNorm();
        get<density>(i)[get<species>(j)] += m_dt;
      }
    }
  }
};


}

#endif
