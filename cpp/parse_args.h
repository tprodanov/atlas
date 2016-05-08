#pragma once

#include <string>
#include <istream>
#include <ostream>
#include <map>
#include <vector>
#include <sstream>

namespace atlas {

std::string base_name(std::string const& name);

struct parse_args {

    parse_args(std::string const& usage);
    void mark_synonyms(std::string const& alias_name, std::string const& name);
    void specify_min_nargs(std::string const& name, size_t nargs);

    void parse(int argc, char const** argv);
    bool flag_up(std::string const& name) const;
    size_t nargs(std::string const& name) const;

    template<class T>
    T get_value(std::string const& name, std::string const& var_name, size_t arg_num) const;
    template<class T>
    T get_value(std::string const& name, std::string const& var_name, size_t arg_num, T const& def) const;

    std::string get_string_value(std::string const& name, std::string const& var_name, size_t arg_num) const;
    std::string get_string_value(std::string const& name, std::string const& var_name,
                                 size_t arg_num, std::string const& def) const;

    std::istream* get_istream(std::string const& name, std::string const& var_name, size_t arg_num,
                              std::istream* def = nullptr) const;
    std::ostream* get_ostream(std::string const& name, std::string const& var_name, size_t arg_num,
                              std::ostream* def = nullptr) const;
    std::vector<std::string> const& leading_words() const;

    void print() const;
    void exit_with_help(std::string const& error, bool stderr = true) const;

private:

    std::string const& true_name(std::string const& name) const;

    static int const error_code;
    std::string usage_;
    std::map<std::string, size_t> min_nargs_;
    std::map<std::string, std::string> synonyms_;
    std::map<std::string, std::vector<std::string>> args_;

};

template<class T>
T parse_args::get_value(std::string const& name, std::string const& var_name, size_t arg_num) const {
    if (nargs(name) < arg_num) {
        exit_with_help(var_name + " is required argument");
    }

    std::ostringstream oss(args_.at(true_name(name))[arg_num - 1]);
    T value;
    oss >> value;
    if (!oss) {
        exit_with_help("Cannot parse " + var_name);
    }
    return value;
}

template<class T>
T parse_args::get_value(std::string const& name, std::string const& var_name, size_t arg_num, T const& def) const {
    if (nargs(name) < arg_num) {
        return def;
    }

    std::istringstream iss(args_.at(true_name(name))[arg_num - 1]);
    T value;
    iss >> value;
    if (!iss) {
        exit_with_help("Cannot parse " + var_name + ": " + args_.at(true_name(name))[arg_num - 1]);
    }
    return value;
}

} // atlas
