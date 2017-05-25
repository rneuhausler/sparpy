/*
 * Python.h
 *
 *  Created on: 28 Nov 2014
 *      Author: mrobins
 */

#ifndef SPARPY_PYTHON_H_
#define SPARPY_PYTHON_H_

#include "sparpy.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

namespace sparpy {

template<typename T, unsigned int D>
struct VectFromPythonList
{

	VectFromPythonList() {
		boost::python::converter::registry::push_back(
				&convertible,
				&construct,
				boost::python::type_id<Vector<T,D> >());
	}

	static void* convertible(PyObject* obj_ptr) {
		if (!PyList_Check(obj_ptr)) return 0;
		if (PyList_Size(obj_ptr) != D) return 0;
		return obj_ptr;
	}

	static void construct(
			PyObject* obj_ptr,
			boost::python::converter::rvalue_from_python_stage1_data* data) {
		const int size = PyList_Size(obj_ptr);

		// Grab pointer to memory into which to construct the new QString
		void* storage = (
				(boost::python::converter::rvalue_from_python_storage<Vector<T,D> >*)
				data)->storage.bytes;

		// construct the new vector using the character data
		// extraced from the python object
        for (int i = 0; i < D; ++i) {
            (*static_cast<Vector<T,D>*>(storage))[i] = extract<T>(PyList_GetItem(obj_ptr,i));
        }

		// Stash the memory chunk pointer for later use by boost.python
		data->convertible = storage;
	}

};

template<typename T, unsigned int D>
struct VectToPython
{
    static PyObject* convert(Vector<T,D> const& s)
      {
    	boost::python::list ret;
        for (int i = 0; i < D; ++i) {
    	    ret.append(s[i]);
        }
        return boost::python::incref(ret.ptr());
      }
};

template <typename T, unsigned int D>
double getitem_vector(Vector<T,D>& v, int index) {
  if(index < 0 || index >=D)
  {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    throw boost::python::error_already_set();;
  }
  return v[index];
}

template <typename T, unsigned int D>
void setitem_vector(Vector<T,D>& v, int index, float val)
{
  if(index < 0 || index >=D)
  {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    throw boost::python::error_already_set();;
  }
  v[index] = val;
}

template <typename ParticlesType>
typename ParticlesType::reference getitem_particles(ParticlesType& p, int index) {
  if(index < 0 || index >= p.size())
  {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    throw boost::python::error_already_set();;
  }
  return p[index];
}

template <typename ParticlesType>
void setitem_particles_from_value(ParticlesType& p, int index, const typename ParticlesType::value_type& val)
{
  if(index < 0 || index >= p.size())
  {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    throw boost::python::error_already_set();;
  }
  p[index] = val;
}

template <typename ParticlesType>
void setitem_particles_from_reference(ParticlesType& p, int index, const typename ParticlesType::reference& ref)
{
  if(index < 0 || index >= p.size())
  {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    throw boost::python::error_already_set();;
  }
  p[index] = ref;
}


template <typename ParticlesType>
void particles_push_back(ParticlesType& p, const typename ParticlesType::value_type& val) {
    p.push_back(val);
}


template<class T>
struct vtkSmartPointer_to_python {
	static PyObject *convert(const vtkSmartPointer<T> &p) {
		std::ostringstream oss;
		oss << (vtkObjectBase*) p.GetPointer();
		std::string address_str = oss.str();

		using namespace boost::python;
		object obj = import("vtk").attr("vtkObjectBase")(address_str);
		//object obj2 = import("tvtk.api").attr("tvtk").attr("to_tvtk")(obj);
		return incref(obj.ptr());
	}
};

//
// This python to C++ converter uses the fact that VTK Python objects have an
// attribute called __this__, which is a string containing the memory address
// of the VTK C++ object and its class name.
// E.g. for a vtkPoints object __this__ might be "_0000000105a64420_p_vtkPoints"
//
void* extract_vtk_wrapped_pointer(PyObject* obj)
{
    char thisStr[] = "__this__";
    //first we need to get the __this__ attribute from the Python Object
    if (!PyObject_HasAttrString(obj, thisStr))
        return NULL;

    PyObject* thisAttr = PyObject_GetAttrString(obj, thisStr);
    if (thisAttr == NULL)
        return NULL;

    char* str = PyString_AsString(thisAttr);
    if(str == 0 || strlen(str) < 1)
        return NULL;

    char hex_address[32], *pEnd;
    char *_p_ = strstr(str, "_p_vtk");
    if(_p_ == NULL) return NULL;
    char *class_name = strstr(_p_, "vtk");
    if(class_name == NULL) return NULL;
    strcpy(hex_address, str+1);
    hex_address[_p_-str-1] = '\0';

    long address = strtol(hex_address, &pEnd, 16);

    vtkObjectBase* vtk_object = (vtkObjectBase*)((void*)address);
    if(vtk_object->IsA(class_name))
    {
        return vtk_object;
    }

    return NULL;
}


template<typename T, typename VT>
const typename T::value_type & get_const(const VT & arg) {
    return get<T>(arg);
}

template<typename T, typename VT>
typename T::value_type & get_non_const(VT & arg) {
    return get<T>(arg);
}

template<typename T, typename VT>
void set_data(VT & arg, const typename T::value_type & data) {
    get<T>(arg) = data;
}


}


#endif /* PYTHON_H_ */
