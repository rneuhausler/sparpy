#include "python.hpp"

using namespace sparpy;

#define VTK_PYTHON_CONVERSION(type) \
    /* register the to-python converter */ \
    to_python_converter<vtkSmartPointer<type>,vtkSmartPointer_to_python<type> >(); \
    /* register the from-python converter */ \
    converter::registry::insert(&extract_vtk_wrapped_pointer, type_id<type>());


BOOST_PYTHON_MODULE(sparpy) {

	VTK_PYTHON_CONVERSION(vtkUnstructuredGrid);

    #define ADD_PROPERTY(name_string, name, D) \
		.add_property(name_string, \
					make_function(&get_non_const<name,ParticlesType<D>::value_type>, \
									return_value_policy<copy_non_const_reference>()), \
					&set_data<name,ParticlesType<D>::value_type>)

    #define ADD_PROPERTY_REF(name_string, name, D) \
		.add_property(name_string, \
					make_function(&get_non_const<name,ParticlesType<D>::reference>, \
									return_value_policy<copy_non_const_reference>()), \
					&set_data<name,ParticlesType<D>::reference>)




    #define ADD_DIMENSION(D) \
        VectFromPythonList<double,D>(); \
        VectFromPythonList<bool,D>();   \
                                        \
        to_python_converter<            \
            Vector<double,D>,           \
            VectToPython<double,D> >(); \
                                        \
        class_<ParticlesType<D>,std::shared_ptr<ParticlesType<D>>>("Particles"#D) \
            .def(init<size_t>())                                                \
            .def("__getitem__", &getitem_particles<ParticlesType<D>>)           \
            .def("__setitem__", &setitem_particles_from_reference<ParticlesType<D>>) \
            .def("__setitem__", &setitem_particles_from_value<ParticlesType<D>>)     \
            .def("__len__", &ParticlesType<D>::size)     \
            .def("init_neighbour_search",&ParticlesType<D>::init_neighbour_search)\
            .def("get_grid",&ParticlesType<D>::get_grid,             \
                    return_value_policy<return_by_value>())                          \
            .def("append",&particles_push_back<ParticlesType<D>>)                          \
            ;                                                                   \
                                                                                    \
        class_<ParticlesType<D>::reference >("ParticleRef"#D,no_init) \
            ADD_PROPERTY_REF("id",id,D)                               \
            ADD_PROPERTY_REF("position",position_d<D>,D)              \
            ADD_PROPERTY_REF("alive",alive,D)                         \
            ADD_PROPERTY_REF("scalar",scalar,D)          \
            ADD_PROPERTY_REF("force",force_d<D>,D)          \
            ;                                                   \
        class_<ParticlesType<D>::value_type >("Particle"#D,init<>()) \
            ADD_PROPERTY("id",id,D)         \
            ADD_PROPERTY("position",position_d<D>,D)     \
            ADD_PROPERTY("alive",alive,D)                \
            ADD_PROPERTY("scalar",scalar,D)              \
            ADD_PROPERTY("force",force_d<D>,D)          \
            ;                                            \
                                                        \
        class_<Simulation<D>>("Simulation"#D,init<>()) \
            .def("add_force", &Simulation<D>::add_force<exponential_force<D>>)   \
            .def("add_action", &Simulation<D>::add_action)   \
            .def("set_domain", &Simulation<D>::set_domain)   \
            .def("add_particles", &Simulation<D>::add_particles)   \
            .def("integrate", &Simulation<D>::integrate)   \
            ;                                            \
                                                        \
        class_<exponential_force<D>>("exponential_force"#D,init<double,double>()) \
            ;                                            \
                                                        \
        class_<hard_sphere<D>>("hard_sphere"#D,init<double>()) \
            ;                                            \


    ADD_DIMENSION(1)
    ADD_DIMENSION(2)
    ADD_DIMENSION(3)
    //.def("copy_from_vtk_grid",&ParticlesType<D>::copy_from_vtk_grid)      \

}
