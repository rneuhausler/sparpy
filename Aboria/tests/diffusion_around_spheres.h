/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#ifndef DIFFUSION_AROUND_SPHERES_H_
#define DIFFUSION_AROUND_SPHERES_H_


#include <cxxtest/TestSuite.h>

#include <random>
#include "Aboria.h"

using namespace Aboria;


class DiffusionAroundSpheres : public CxxTest::TestSuite {
public:
    typedef std::mt19937 generator_type;
    generator_type generator;


    template<template <typename> class SearchMethod>
	void helper_diffusion_around_spheres(void) {
		//const double tol = GEOMETRY_TOLERANCE;

		ABORIA_VARIABLE(radius,double,"radius")

		typedef Particles<std::tuple<radius>,3,std::vector,SearchMethod> spheres_type;
		typedef Particles<> points_type;
        typedef position_d<3> position;
		spheres_type spheres;

		const double L = 10.0;
		const double D = 1.0;
		const double dt = 0.1;
		const double timesteps = 100;

		spheres.push_back(double3(0,0,0));
		get<radius>(spheres[0]) = 1.0;
		spheres.push_back(double3(5,0,0));
		get<radius>(spheres[1]) = 2.0;
		spheres.push_back(double3(0,-5,0));
		get<radius>(spheres[2]) = 1.5;
		spheres.push_back(double3(0,0,5));
		get<radius>(spheres[3]) = 1.0;

    	spheres.init_neighbour_search(double3(-L,-L,-L),double3(L,L,L),bool3(true,true,true));

		points_type points;
		std::uniform_real_distribution<double> uni(-L,L);
		for (int i = 0; i < 1000; ++i) {
			points.push_back(double3(uni(generator),uni(generator),uni(generator)));
		}

    	points.init_neighbour_search(double3(-L,-L,-L),double3(L,L,L),bool3(true,true,true));

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;
        Symbol<radius> r;
        auto a_s = create_label<0>(spheres);
        auto b_s = create_label<1>(spheres);
        auto a_p = create_label<0>(points);
        auto b_p = create_label<1>(points);
		auto dx = create_dx(a_p,b_s);

		Normal N;
		VectorSymbolic<double,3> vector;		
        AccumulateWithinDistance<std::bit_or<bool> > any(2);
        AccumulateWithinDistance<std::plus<double3> > sum(2);

		/*
		 * Kill any points within spheres
		 */
		alive_[a_p] = !any(b_s, if_else(norm(dx) <= r[b_s],true,false));

		/*
		 * Check no points within spheres
		 */
		for(auto i: points) {
			TS_ASSERT_EQUALS(get<alive>(i), true);
			TS_ASSERT_RELATION(std::greater<double>, (get<position>(i) - double3(0,0,0)).norm(), 1.0);
			TS_ASSERT_RELATION(std::greater<double>, (get<position>(i) - double3(5,0,0)).norm(), 2.0);
			TS_ASSERT_RELATION(std::greater<double>, (get<position>(i) - double3(0,-5,0)).norm(), 1.5);
			TS_ASSERT_RELATION(std::greater<double>, (get<position>(i) - double3(0,0,5)).norm(), 1.0);
		}

		/*
		 * Diffusion step for points and "reflect" off spheres
		 */
		for (int i = 0; i < timesteps; ++i) {
			p[a_p] += std::sqrt(2*D*dt)*vector(N[a_p],N[a_p],N[a_p]);
            p[a_p] += sum(b_s, if_else(norm(dx) < r[b_s]
                                    ,-2*(r[b_s]/norm(dx)-1)*dx
                                    ,0));
            /*
            const double3 pos = get<position>(points[0]);
            if (((pos - double3(0,0,0)).norm() < 1.0) || 
                ((pos - double3(5,0,0)).norm() < 2.0) ||
                ((pos - double3(0,-5,0)).norm() < 1.5)||
                ((pos - double3(0,0,5)).norm() < 1.0)) {
                std::cout << "BOUNCE" << std::endl;
                std::cout << "after step "<<get<position>(points[0])<<std::endl;
                if ((pos - double3(0,0,0)).norm() < 1.0) {
                    std::cout << "dx "<<pos-double3(0,0,0)<<std::endl;
                    std::cout << "|dx| "<<(pos-double3(0,0,0)).norm()<<std::endl;
                }
                if ((pos - double3(5,0,0)).norm() < 2.0) {
                    std::cout << "dx "<<pos-double3(5,0,0)<<std::endl;
                    std::cout << "|dx| "<<(pos-double3(5,0,0)).norm()<<std::endl;
                }
                if ((pos - double3(0,-5,0)).norm() < 1.5) {
                    std::cout << "dx "<<pos-double3(0,-5,0)<<std::endl;
                    std::cout << "|dx| "<<(pos-double3(0,-5,0)).norm()<<std::endl;
                }
                if ((pos - double3(0,0,5)).norm() < 1.0) {
                    std::cout << "dx "<<pos-double3(0,0,5)<<std::endl;
                    std::cout << "|dx| "<<(pos-double3(0,0,5)).norm()<<std::endl;
                }
                p[a_p] += sum(b_s, norm(dx) < r[b_s],
                    -2*(r[b_s]/norm(dx)-1)*dx );
                std::cout << "after bounce "<<get<position>(points[0])<<std::endl;
                const double3 pos2 = get<position>(points[0]);
                if (((pos2 - double3(0,0,0)).norm() < 1.0) || 
                ((pos2 - double3(5,0,0)).norm() < 2.0) ||
                ((pos2 - double3(0,-5,0)).norm() < 1.5)||
                ((pos2 - double3(0,0,5)).norm() < 1.0)) {
                    std::cout << "STILL IN" << std::endl;
                    if ((pos2 - double3(0,0,0)).norm() < 1.0) {
                        std::cout << "dx "<<pos2-double3(0,0,0)<<std::endl;
                        std::cout << "|dx| "<<(pos2-double3(0,0,0)).norm()<<std::endl;
                    }
                    if ((pos2 - double3(5,0,0)).norm() < 2.0) {
                        std::cout << "dx "<<pos2-double3(5,0,0)<<std::endl;
                        std::cout << "|dx| "<<(pos2-double3(5,0,0)).norm()<<std::endl;
                    }
                    if ((pos2 - double3(0,-5,0)).norm() < 1.5) {
                        std::cout << "dx "<<pos2-double3(0,-5,0)<<std::endl;
                        std::cout << "|dx| "<<(pos2-double3(0,-5,0)).norm()<<std::endl;
                    }
                    if ((pos2 - double3(0,0,5)).norm() < 1.0) {
                        std::cout << "dx "<<pos2-double3(0,0,5)<<std::endl;
                        std::cout << "|dx| "<<(pos2-double3(0,0,5)).norm()<<std::endl;
                    }
                }

            }
            */

		}

		/*
		 * Check still no points within spheres
		 */
		for(auto i: points) {
			TS_ASSERT_RELATION(std::greater<double>, (get<position>(i) - double3(0,0,0)).norm(), 1.0);
			TS_ASSERT_RELATION(std::greater<double>, (get<position>(i) - double3(5,0,0)).norm(), 2.0);
			TS_ASSERT_RELATION(std::greater<double>, (get<position>(i) - double3(0,-5,0)).norm(), 1.5);
			TS_ASSERT_RELATION(std::greater<double>, (get<position>(i) - double3(0,0,5)).norm(), 1.0);
		}
	}


    void test_bucket_search_parallel() {
        helper_diffusion_around_spheres<bucket_search_parallel>();
    }

    void test_bucket_search_serial() {
        helper_diffusion_around_spheres<bucket_search_serial>();
    }

    void test_nanoflann_adaptor() {
        helper_diffusion_around_spheres<nanoflann_adaptor>();
    }

};

#endif /* DIFFUSION_AROUND_SPHERES_H_ */
