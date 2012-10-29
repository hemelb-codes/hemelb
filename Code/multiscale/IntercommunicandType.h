// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_MULTISCALE_INTERCOMMUNICANDTYPE_H
#define HEMELB_MULTISCALE_INTERCOMMUNICANDTYPE_H
#include <string>
#include <vector>

namespace hemelb
{
  namespace multiscale
  {
    // Different implementations of intercommunicator will use different systems to specify type at runtime
    // template on this choice.
    // Example: MPI Datatype
    /***
     * This class represents run-time type information for the multiscale system to know which shared values in a given intercommunicand
     * are of which types.
     * It uses a vector of pairs of string labels and runtime type values, using some runtime type implementation.
     * The string labels are used to identify the shared values to the intercommunication implementation.
     * A vector is used (rather than a map) as the order of values corresponds to the order with which the shared values are registered in an intercommunicand.
     * This could, for example, be MPI_DATATYPE.
     * @tparam RuntimeTypeImplementation A choice of how to represent types at runtime, such as MPI datatype.
     */
    template<class RuntimeTypeImplementation> class IntercommunicandType
    {
      public:
       
        IntercommunicandType(const std::string & alabel): fields(), label(alabel){
 hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Icand created with label %s", alabel.c_str());
}
        
        /***
         * The order of the fields in the intercommunicand type corresponds to the order they are registered in the intercommunicand.
         * @return vector of pairs of string labels and runtime type values.
         */
        std::vector<std::pair<std::string, typename RuntimeTypeImplementation::RuntimeType> > &Fields(){return fields;}
        
        void RegisterSharedValue(const std::string &label,typename RuntimeTypeImplementation::RuntimeType type)
         {
            fields.push_back(std::make_pair(label,type));
             hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("shared value created with label %s", label.c_str());
         }
         template<class T> void RegisterSharedValue(const std::string &label){
            fields.push_back(std::make_pair(label,RuntimeTypeImplementation::template GetType<T>()));
         } 
      private:
        std::vector<std::pair<std::string, typename RuntimeTypeImplementation::RuntimeType> >fields;
        std::string label;
    };
  }
}

#endif // HEMELB_MULTISCALE_INTERCOMMUNICAND_H
