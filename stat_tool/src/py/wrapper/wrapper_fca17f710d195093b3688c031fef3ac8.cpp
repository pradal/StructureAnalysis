#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::bidirectional_iterator_tag const volatile * get_pointer<struct ::std::bidirectional_iterator_tag const volatile >(struct ::std::bidirectional_iterator_tag const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_fca17f710d195093b3688c031fef3ac8()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< struct ::std::bidirectional_iterator_tag, autowig::Held< struct ::std::bidirectional_iterator_tag >::Type, boost::python::bases< struct ::std::forward_iterator_tag > > class_fca17f710d195093b3688c031fef3ac8("BidirectionalIteratorTag", "", boost::python::no_init);
    class_fca17f710d195093b3688c031fef3ac8.def(boost::python::init<  >(""));
    class_fca17f710d195093b3688c031fef3ac8.def(boost::python::init< struct ::std::bidirectional_iterator_tag const & >(""));

    if(autowig::Held< struct ::std::bidirectional_iterator_tag >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< struct ::std::bidirectional_iterator_tag >::Type, autowig::Held< struct ::std::forward_iterator_tag >::Type >();
    }

}