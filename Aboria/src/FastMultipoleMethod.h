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

#ifndef FAST_MULTIPOLE_METHOD_H_
#define FAST_MULTIPOLE_METHOD_H_

#include "detail/FastMultipoleMethod.h"

namespace Aboria {

   
template <typename Expansions, typename NeighbourQuery>
class FastMultipoleMethod {
    typedef typename NeighbourQuery::traits_type traits_type;
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::pointer pointer;
    typedef typename NeighbourQuery::child_iterator child_iterator;
    typedef typename Expansions::expansion_type expansion_type;
    typedef typename traits_type::template vector_type<expansion_type>::type storage_type;
    typedef typename traits_type::template vector_type<
        child_iterator
        //typename std::remove_const<child_iterator>::type
        >::type child_iterator_vector_type;
    typedef typename traits_type::template vector_type<child_iterator_vector_type>::type connectivity_type;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef typename particle_iterator::reference particle_reference;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;
    storage_type m_W;
    storage_type m_g;
    connectivity_type m_connectivity; 
    const NeighbourQuery *m_query;
    Expansions m_expansions;

public:

    FastMultipoleMethod(const NeighbourQuery &query, const Expansions& expansions):
        m_query(&query),m_expansions(expansions)
    {}

    template <typename VectorType>
    expansion_type& calculate_dive_P2M_and_M2M(const child_iterator& ci, 
                                               const VectorType& source_vector) {
        const size_t my_index = m_query->get_bucket_index(*ci);
        const box_type& my_box = m_query->get_bounds(ci);
        LOG(3,"calculate_dive_P2M_and_M2M with bucket "<<my_box);
        expansion_type& W = m_W[my_index];
        std::fill(std::begin(W),std::end(W),0.0);
        if (m_query->is_leaf_node(*ci)) { // leaf node
            detail::calculate_P2M(W, my_box, 
                    m_query->get_bucket_particles(*ci),
                    source_vector,m_query->get_particles_begin(),m_expansions);
        } else { 
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                expansion_type& child_W = calculate_dive_P2M_and_M2M(cj,source_vector);
                const box_type& child_box = m_query->get_bounds(cj);
                m_expansions.M2M(W,my_box,child_box,child_W);
            }
        }
        return W;
    }

    template <typename VectorType>
    void calculate_dive_M2L_and_L2L(
            VectorType& target_vector,
            const child_iterator_vector_type& connected_buckets_parent,
            const expansion_type& g_parent, 
            const box_type& box_parent, 
            const child_iterator& ci,
            const VectorType& source_vector) {
        const box_type& target_box = m_query->get_bounds(ci);
        LOG(3,"calculate_dive_M2L_and_L2L with bucket "<<target_box);
        size_t target_index = m_query->get_bucket_index(*ci);
        expansion_type& g = m_g[target_index];
        std::fill(std::begin(g),std::end(g),0.0);
        typename connectivity_type::reference 
            connected_buckets = m_connectivity[target_index];
        connected_buckets.clear();
        
        detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);

        if (connected_buckets_parent.empty()) {
            for (child_iterator cj = m_query->get_children(); cj != false; ++cj) {
                const box_type& source_box = m_query->get_bounds(cj);
                if (theta.check(source_box.bmin,source_box.bmax)) {
                    connected_buckets.push_back(cj);
                } else {
                    size_t source_index = m_query->get_bucket_index(*cj);
                    m_expansions.M2L(g,target_box,source_box,m_W[source_index]);
                }
            }
        } else {
            // expansion from parent
            m_expansions.L2L(g,target_box,box_parent,g_parent);

            // expansions from weakly connected buckets on this level
            // and store strongly connected buckets to connectivity list
            for (const child_iterator& source: connected_buckets_parent) {
                if (m_query->is_leaf_node(*source)) {
                    connected_buckets.push_back(source);
                } else {
                    for (child_iterator cj = m_query->get_children(source); cj != false; ++cj) {
                        const box_type& source_box = m_query->get_bounds(cj);
                        if (theta.check(source_box.bmin,source_box.bmax)) {
                            connected_buckets.push_back(cj);
                        } else {
                            size_t source_index = m_query->get_bucket_index(*cj);
                            m_expansions.M2L(g,target_box,source_box,m_W[source_index]);
                        }
                    }
                }
            }
        }
        if (!m_query->is_leaf_node(*ci)) { // leaf node
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                calculate_dive_M2L_and_L2L(target_vector,connected_buckets,g
                                            ,target_box,cj,source_vector);
            }
        } else if (target_vector.size() > 0) {
            detail::calculate_L2P(target_vector,g,target_box,
                    m_query->get_bucket_particles(*ci),
                    m_query->get_particles_begin(),m_expansions);

            for (child_iterator& cj: connected_buckets) { 
                if (m_query->is_leaf_node(*cj)) {
                    LOG(3,"calculate_P2P: target = "<<target_box<<" source = "<<m_query->get_bounds(cj));
                    detail::calculate_P2P(target_vector,source_vector,
                        m_query->get_bucket_particles(*ci),m_query->get_bucket_particles(*cj),
                        m_query->get_particles_begin(),m_query->get_particles_begin(),
                        m_expansions);
                } else {
                    for (reference ref_j: m_query->get_subtree(cj)) {
                        if (m_query->is_leaf_node(ref_j)) {
                            detail::calculate_P2P(target_vector,source_vector,
                                m_query->get_bucket_particles(*ci),m_query->get_bucket_particles(ref_j),
                                m_query->get_particles_begin(),m_query->get_particles_begin(),
                                m_expansions);
                        }
                    }
                }
            }
        }
    }


    // note: assumes identical source and target particle sets
    // target_vector += A*source_vector
    template <typename VectorType>
    void gemv(VectorType& target_vector, const VectorType& source_vector) {
        CHECK(target_vector.size() == source_vector.size(), "source and target vector not same length")

        const size_t n = m_query->number_of_buckets();
        m_W.resize(n);
        m_g.resize(n);
        m_connectivity.resize(n);

        // upward sweep of tree
        //
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            calculate_dive_P2M_and_M2M(ci,source_vector);
        }

        // downward sweep of tree. 
        //
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            child_iterator_vector_type dummy;
            expansion_type g = {};
            calculate_dive_M2L_and_L2L(target_vector,dummy,g,box_type(),ci,source_vector);
        }
    }


    template <typename VectorType>
    void calculate_expansions(const VectorType& source_vector) {

        const size_t n = m_query->number_of_buckets();
        m_W.resize(n);
        m_g.resize(n);
        m_connectivity.resize(n);

        // upward sweep of tree
        //
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            calculate_dive_P2M_and_M2M(ci,source_vector);
        }

        // downward sweep of tree.
        //
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            child_iterator_vector_type dummy;
            VectorType dummy2;
            expansion_type g = {};
            calculate_dive_M2L_and_L2L(dummy2,dummy,g,box_type(),ci,source_vector);
        }
    }


    // evaluate expansions for given point
    template <typename VectorType>
    double evaluate_expansion(const Vector<double,dimension>& p, const VectorType& source_vector) {
        pointer bucket;
        box_type box;
        m_query->get_bucket(p,bucket,box);
        LOG(3,"evaluating expansion at point "<<p<<" with box "<<box);
        const size_t index = m_query->get_bucket_index(*bucket); 

        double sum = Expansions::L2P(p,box,m_g[index]);
        for (child_iterator& ci: m_connectivity[index]) { 
            if (m_query->is_leaf_node(*ci)) {
                sum += detail::calculate_P2P_position(p
                    ,m_query->get_bucket_particles(*ci)
                    ,m_expansions,source_vector,m_query->get_particles_begin());
            } else {
                for (reference subtree_reference: m_query->get_subtree(ci)) {
                    if (m_query->is_leaf_node(subtree_reference)) {
                        sum += detail::calculate_P2P_position(p
                                ,m_query->get_bucket_particles(subtree_reference)
                                ,m_expansions,source_vector,m_query->get_particles_begin());
                    }
                }
            }
        }
        return sum;
    }

};


template <unsigned int D, unsigned int N, typename Function> 
detail::BlackBoxExpansions<D,N,Function> make_black_box_expansion(const Function& function) {
    return detail::BlackBoxExpansions<D,N,Function>(function);
}


template <typename Expansions, typename NeighbourQuery>
FastMultipoleMethod<Expansions,NeighbourQuery>
make_fmm_query(const NeighbourQuery &query, const Expansions& expansions) {
    return FastMultipoleMethod<Expansions,NeighbourQuery>(query,expansions);
}

}

#endif

