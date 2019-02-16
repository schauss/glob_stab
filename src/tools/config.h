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

#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <boost/program_options.hpp>

//! Return instance
#define CONFIG (Config::instance())
//! Return variable with given name as given type
#define CONFIG_VAR(name,type)    (CONFIG->var(#name)->as<type>())
//! Check whether config with given name is set
#define CONFIG_ISSET(name)       (CONFIG->varIsSet(#name))

/*!
 * \brief The Config class is a simple wrapper class for boost::program_options which provides a global configuration.
 *
 * It is implemented as singleton.
 * A pointer to the instance can be aquired using Config::instance() or the macro CONFIG.
 * Convenience macros exist to easily read configuration values:
 *  * CONFIG_VAR(name, type) => Returns configuration value with "name" as "type"
 *  * CONFIG_ISSET(name)     => Check whether config with given name is set
 */
class Config
{
public:
    //! Return instance of Config
    static Config* instance();
    //! Call destructor
    static void clear();
    const boost::program_options::variable_value* var(const std::string &name) const;
    bool varIsSet(const std::string &name) const;

    /*!
     * \brief Add program_options to Config
     * \param opt Options description, may describe several program options
     */
    void addOptions(const boost::program_options::options_description &opt);
    /*!
     * \brief Parse command line options
     * \param argc Number of command-line arguments
     * \param argv Array of pointers to char (i.e., array of argc command-line strings)
     */
    void parseOptions(int argc, char *argv[]);
    //! Print all configuration values to std::cout
    void printOptions() const;
    //! Print program options (i.e., usage information)
    void printDescription() const;

private:
    //! Private constructor (singleton)
    Config() {}
    //! Private destructor (singleton)
    Config(Config const&) {}
    //! Assignment operator => Returns instance
    Config& operator=(Config const &) {return *config;}
    //! The instance.
    static Config* config;

    //! Map from program_option-name to value of program_option
    boost::program_options::variables_map vm;
    //! Options description
    boost::program_options::options_description options;
};

#include <stdexcept>
//! Exception which is raised if there is a configuration error
class ConfigException : public std::runtime_error
{
public:
    ConfigException() : std::runtime_error("Configuration Error") { }
    ConfigException(const std::string& what_arg) : std::runtime_error(what_arg) { }
};

#endif // CONFIG_H
