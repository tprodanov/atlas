#include "parse_args.h"

#include <iostream>
#include <fstream>

namespace atlas {

std::string base_name(const std::string& name) {
    size_t i = name.find_last_of('/');

    // Fail
    if (i == std::string::npos) {
        if (name[0] == '/')
            return name.substr(1);
        return name;
    }

    return name.substr(i + 1);
}

int const parse_args::error_code = 64;

parse_args::parse_args(std::string const& usage)
    : usage_(usage)
{
    mark_synonyms("h", "help");
}

void parse_args::mark_synonyms(std::string const& alias_name, std::string const& name) {
    synonyms_[alias_name] = name;
}

void parse_args::specify_min_nargs(std::string const& name, size_t nargs) {
    min_nargs_[true_name(name)] = nargs;
}

bool parse_args::flag_up(std::string const& name) const {
    return args_.find(true_name(name)) != args_.end();
}

size_t parse_args::nargs(std::string const& name) const {
    auto it = args_.find(true_name(name));
    if (it == args_.end())
        return 0;
    return it->second.size();
}

std::string const& parse_args::true_name(std::string const& name) const {
    auto it = synonyms_.find(name);
    if (it == synonyms_.end())
        return name;
    return it->second;
}

void parse_args::exit_with_help(std::string const& error, bool stderr) const {
    (stderr ? std::cerr : std::cout) << error << "\n\n" << usage_;
    exit(stderr ? error_code : 0);
}

std::istream* parse_args::get_istream(std::string const& name, std::string const& var_name, size_t arg_num,
                                      std::istream* def) const {
    if (nargs(name) < arg_num) {
        if (def) {
            return def;
        } else {
            exit_with_help(var_name + " is required argument");
        }
    }

    std::string const& filename = args_.at(true_name(name))[arg_num - 1];
    std::istream *stream = new std::ifstream(filename);
    if (!*stream) {
        std::cerr << "Cannot open " << var_name << ": " << filename << std::endl;
        exit(error_code);
    }
    return stream;
}

std::ostream* parse_args::get_ostream(std::string const& name, std::string const& var_name, size_t arg_num,
                                      std::ostream* def) const {
    if (nargs(name) < arg_num) {
        if (def) {
            return def;
        } else {
            exit_with_help(var_name + " is required argument");
        }
    }

    std::string const& filename = args_.at(true_name(name))[arg_num - 1];
    std::ostream *stream = new std::ofstream(filename);
    if (!*stream) {
        std::cerr << "Cannot write to " << var_name << ": " << filename << std::endl;
        exit(error_code);
    }
    return stream;
}


std::string parse_args::get_string_value(std::string const& name, std::string const& var_name,
                                         size_t arg_num) const {
    if (nargs(name) < arg_num) {
        exit_with_help(var_name + " is required argument");
    }
    return args_.at(true_name(name))[arg_num - 1];
}

std::string parse_args::get_string_value(std::string const& name, std::string const&,
                                         size_t arg_num, std::string const& def) const {
    if (nargs(name) < arg_num) {
        return def;
    }
    return args_.at(true_name(name))[arg_num - 1];
}

std::vector<std::string> const& parse_args::leading_words() const {
    return args_.at("");
}

void parse_args::parse(int argc, char const** argv) {
    std::string current_flag = "";
    args_[current_flag];

    char short_flag[2] {'\0', '\0'};
    for (int i = 0; i < argc; ++i) {
        char const* current_arg = argv[i];

        if (current_arg[0] != '-') {
            args_[current_flag].push_back(current_arg);
            continue;
        } else if (current_arg[1] != '-') {
            ++current_arg;
            while (*current_arg) {
                short_flag[0] = *current_arg;
                current_flag = true_name(short_flag);
                args_[current_flag];
                ++current_arg;
            }
        } else {
            current_flag = current_arg + 2;
            args_[current_flag];
        }

        auto it = min_nargs_.find(current_flag);
        if (it != min_nargs_.end()) {
            size_t nargs = it->second;
            while (i + 1 < argc && nargs) {
                ++i;
                args_[current_flag].push_back(argv[i]);
                --nargs;
            }
        }
    }
    if (current_flag == "" || flag_up("help")) {
        std::cout << usage_;
        exit(0);
    }
}

void parse_args::print() const {
    for (auto const& entry : args_) {
        std::cout << entry.first << ": ";
        for (std::string const& arg : entry.second) {
            std::cout << arg << " ";
        }
        std::cout << std::endl;
    }
}

} // atlas
