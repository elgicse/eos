/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <map>
#include <random>
#include <vector>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <yaml-cpp/yaml.h>

#include <iostream>

namespace fs = boost::filesystem;

namespace eos
{
    bool
    operator== (const ParameterDescription & lhs, const ParameterDescription & rhs)
    {
    	if (lhs.min != rhs.min)
    		return false;
    	if (lhs.max != rhs.max)
    		return false;
    	if (lhs.nuisance!= rhs.nuisance)
    		return false;
    	if (lhs.parameter->name() != rhs.parameter->name())
    		return false;

    	return true;
    }

    struct Parameter::Template
    {
        std::string name;

        double min, central, max;
    };

    struct Parameter::Data :
        Parameter::Template
    {
        double value;

        Parameter::Id id;

        Data(const Parameter::Template & t, const Parameter::Id & i) :
            Parameter::Template(t),
            value(t.central),
            id(i)
        {
        }
    };

    struct Parameters::Data
    {
        std::vector<Parameter::Data> data;
    };

    template <>
    struct WrappedForwardIteratorTraits<Parameters::IteratorTag>
    {
        typedef std::vector<Parameter>::iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<Parameters::IteratorTag, Parameter>;

    template <>
    struct Implementation<Parameters>
    {
        std::shared_ptr<Parameters::Data> parameters_data;

        std::map<std::string, unsigned> parameters_map;

        std::vector<Parameter> parameters;

        Implementation(const std::initializer_list<Parameter::Template> & list) :
            parameters_data(new Parameters::Data)
        {
            unsigned idx(0);
            for (auto i(list.begin()), i_end(list.end()) ; i != i_end ; ++i, ++idx)
            {
                parameters_data->data.push_back(Parameter::Data(*i, idx));
                parameters_map[i->name] = idx;
                parameters.push_back(Parameter(parameters_data, idx));
            }
        }

        Implementation(const Implementation & other) :
            parameters_data(new Parameters::Data(*other.parameters_data)),
            parameters_map(other.parameters_map)
        {
            parameters.reserve(other.parameters.size());
            for (unsigned i = 0 ; i != parameters.size() ; ++i)
            {
                parameters.push_back(Parameter(parameters_data, i));
            }
        }

        void
        load_defaults()
        {
            fs::path base;
            if (std::getenv("EOS_TESTS_PARAMETERS"))
            {
                std::string envvar = std::string(std::getenv("EOS_TESTS_PARAMETERS"));
                base = fs::system_complete(envvar);
            }
            else
            {
                base = fs::system_complete(EOS_DATADIR "/eos/parameters/");
            }

            if (! fs::exists(base))
            {
                throw InternalError("Could not find the parameter input files");
            }

            if (! fs::is_directory(base))
            {
                throw InternalError("Expect '" + base.string() + " to be a directory");
            }

            unsigned idx = parameters.size();
            for (fs::directory_iterator f(base), f_end ; f != f_end ; ++f)
            {
                if (! fs::is_regular_file(f->status()))
                    continue;

                auto file_path = f->path();

                if (".yaml" != file_path.extension().string())
                    continue;

                std::string file = file_path.string();
                try
                {
                    YAML::Node node = YAML::LoadFile(file);

                    for (auto && p : node)
                    {
                        std::string name = p.first.Scalar();
                        double central, min, max;

                        if (! p.second["central"])
                        {
                            throw ParameterInputFileNodeError(file, name, "has no entry named 'central'");
                        }
                        else if (YAML::NodeType::Scalar != p.second["central"].Type())
                        {
                            throw ParameterInputFileNodeError(file, name + ".central", "is not a scalar");
                        }
                        central = p.second["central"].as<double>();

                        if (! p.second["min"])
                        {
                            throw ParameterInputFileNodeError(file, name, "has no entry named 'min'");
                        }
                        else if (YAML::NodeType::Scalar != p.second["min"].Type())
                        {
                            throw ParameterInputFileNodeError(file, name, "is not a scalar");
                        }
                        min = p.second["min"].as<double>();

                        if (! p.second["max"])
                        {
                            throw ParameterInputFileNodeError(file, name, "has no entry named 'max'");
                        }
                        else if (YAML::NodeType::Scalar != p.second["max"].Type())
                        {
                            throw ParameterInputFileNodeError(file, name, "is not a scalar");
                        }
                        max = p.second["max"].as<double>();

                        parameters_data->data.push_back(Parameter::Data(Parameter::Template { name, min, central, max }, idx));
                        parameters_map[name] = idx;
                        parameters.push_back(Parameter(parameters_data, idx));

                        ++idx;
                    }
                }
                catch (std::exception & e)
                {
                    throw ParameterInputFileParseError(file, e.what());
                }
            }
        }
    };

    Parameters::Parameters(Implementation<Parameters> * imp) :
        PrivateImplementationPattern<Parameters>(imp)
    {
    }

    Parameters::~Parameters()
    {
    }

    Parameters
    Parameters::clone() const
    {
        return Parameters(new Implementation<Parameters>(*_imp));
    }

    Parameter
    Parameters::operator[] (const std::string & name) const
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        return Parameter(_imp->parameters_data, i->second);
    }

    Parameter
    Parameters::operator[] (const Parameter::Id & id) const
    {
        if (id >= _imp->parameters.size())
            throw InternalError("Parameters::operator[] (Parameter::Id): invalid id '" + stringify(id) + "'");

        return _imp->parameters[id];
    }

    Parameter
    Parameters::declare(const std::string & name, double value)
    {
        // return existing parameter
        auto i(_imp->parameters_map.find(name));
        if (_imp->parameters_map.end() != i)
            return Parameter(_imp->parameters_data, i->second);

        // create new parameter
        unsigned idx = _imp->parameters.size();
        _imp->parameters_data->data.push_back(Parameter::Data(Parameter::Template { name, value, value, value }, idx));
        _imp->parameters_map[name] = idx;
        _imp->parameters.push_back(Parameter(_imp->parameters_data, idx));

        return _imp->parameters.back();
    }

    void
    Parameters::set(const std::string & name, const double & value)
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        _imp->parameters_data->data[i->second].value = value;
    }

    Parameters::Iterator
    Parameters::begin() const
    {
        return Parameters::Iterator(_imp->parameters.begin());
    }

    Parameters::Iterator
    Parameters::end() const
    {
        return Parameters::Iterator(_imp->parameters.end());
    }

    bool
    Parameters::operator!= (const Parameters & rhs) const
    {
        return rhs._imp.get() != this->_imp.get();
    }

    Parameters
    Parameters::Defaults()
    {
        auto imp = new Implementation<Parameters>{};
        imp->load_defaults();

        return Parameters(imp);
    }

    Parameter::Parameter(const std::shared_ptr<Parameters::Data> & parameters_data, unsigned index) :
        _parameters_data(parameters_data),
        _index(index)
    {
    }

    Parameter::Parameter(const Parameter & other) :
        _parameters_data(other._parameters_data),
        _index(other._index)
    {
    }

    Parameter::~Parameter()
    {
    }

    MutablePtr
    Parameter::clone() const
    {
        return MutablePtr(new Parameter(_parameters_data, _index));
    }

    Parameter::operator double () const
    {
        return _parameters_data->data[_index].value;
    }

    double
    Parameter::operator() () const
    {
        return _parameters_data->data[_index].value;
    }

    double
    Parameter::evaluate() const
    {
        return _parameters_data->data[_index].value;
    }

    const Parameter &
    Parameter::operator= (const double & value)
    {
        _parameters_data->data[_index].value = value;

        return *this;
    }

    void
    Parameter::set(const double & value)
    {
        _parameters_data->data[_index].value = value;
    }

    const double &
    Parameter::central() const
    {
        return _parameters_data->data[_index].central;
    }

    const double &
    Parameter::max() const
    {
        return _parameters_data->data[_index].max;
    }

    const double &
    Parameter::min() const
    {
        return _parameters_data->data[_index].min;
    }

    const std::string &
    Parameter::name() const
    {
        return _parameters_data->data[_index].name;
    }

    Parameter::Id
    Parameter::id() const
    {
        return _parameters_data->data[_index].id;
    }

    template <>
    struct WrappedForwardIteratorTraits<ParameterUser::ConstIteratorTag>
    {
        typedef std::set<Parameter::Id>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<ParameterUser::ConstIteratorTag, const Parameter::Id>;

    ParameterUser::ConstIterator
    ParameterUser::begin() const
    {
        return ConstIterator(_ids.cbegin());
    }

    ParameterUser::ConstIterator
    ParameterUser::end() const
    {
        return ConstIterator(_ids.cend());
    }

    void
    ParameterUser::drop(const Parameter::Id & id)
    {
        _ids.erase(id);
    }

    void
    ParameterUser::uses(const Parameter::Id & id)
    {
        _ids.insert(id);
    }

    void
    ParameterUser::uses(const ParameterUser & other)
    {
        _ids.insert(other._ids.cbegin(), other._ids.cend());
    }

    UsedParameter::UsedParameter(const Parameter & parameter, ParameterUser & user) :
        Parameter(parameter)
    {
        user.uses(parameter.id());
    }

    UnknownParameterError::UnknownParameterError(const std::string & name) throw () :
        Exception("Unknown parameter: '" + name + "'")
    {
    }

    ParameterInputFileParseError::ParameterInputFileParseError(const std::string & file, const std::string & msg) throw () :
        Exception("Malformed parameter input file '" + file + "': " + msg)
    {
    }

    ParameterInputFileNodeError::ParameterInputFileNodeError(const std::string & file, const std::string & node, const std::string & msg) throw () :
        Exception("Malformed parameter input file '" + file + "': Node '" + node + "' " + msg)
    {
    }
}
