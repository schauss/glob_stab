/*
    Copyright (C) 2017 Thomas Schauss

    This file is part of glob_stab.

    glob_stab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    glob_stab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with glob_stab. If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#include <boost/program_options.hpp>
#include <iostream>

namespace po = boost::program_options;

Config* Config::config = 0;

Config* Config::instance()
{
    if (!config) {
        config = new Config();
        boost::program_options::options_description opt("Help");
        opt.add_options()
            ("help,h", "help message")
        ;
        config->addOptions(opt);
    }

    return config;
}

void Config::clear()
{
    if (config) {
        delete config;
        config = 0;
    }
}

void Config::addOptions(const boost::program_options::options_description &opt)
{
    options.add(opt);
}

void Config::parseOptions(int argc, char *argv[])
{
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        printDescription();
        exit(1);
    }
    /*
    if (vm.count("verbosity")) {
        std::cout << "Verbosity was set to "
     << vm["verbosity"].as<int>() << ".\n";
    } else {
        std::cout << "Verbosity was not set.\n";
    }
    */
}

void Config::printDescription() const
{
    std::cout << options << "\n";
}

bool Config::varIsSet(const std::string &name) const
{
    return (vm.count(name)>0);
}

const po::variable_value* Config::var(const std::string &name) const
{
    if (vm.count(name))
        return (&(vm[name]));
    else {
        std::cerr << "Trying to read variable \"" << name << "\"!\n";
        throw ConfigException("Required configuration option not present!");
    }
}

void Config::printOptions() const
{
    if (vm.empty()) {
        std::cout << "Configuration empty!\n";
        return;
    }

    std::cout << "Configuration:\n";

    for (boost::program_options::variables_map::const_iterator var = vm.begin(); var != vm.end(); var++) {
        std::cout << var->first << ": ";
        boost::any val = var->second.value();
        if (val.type() == typeid(int))
            std::cout << var->second.as<int>() << std::endl;
        else if (val.type() == typeid(double))
            std::cout << var->second.as<double>() << std::endl;
        else if (val.type() == typeid(std::string))
            std::cout << "\"" << var->second.as<std::string>() << "\"" << std::endl;
        else
            std::cout << "Unknown type!\n";
    }
    std::cout << std::endl;
}
