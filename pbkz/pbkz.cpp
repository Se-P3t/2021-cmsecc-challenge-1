/*
Run progressive BKZ for a given lattice

Reference: bkztest.cpp - solve_lattice_challenge
 */
#include <cxxopts.hpp>

#include "lattice/pbkz.hpp"



int main(int argc, char** argv)
{
    cxxopts::Options options(argv[0], "Run progressive BKZ for a given lattice");

    options.add_options()
        ("matrix", "matrix file", cxxopts::value<std::string>())
        ("block-size", "an integer from 2 to `dim`", cxxopts::value<int>())
        ("verbose", "verbose level", cxxopts::value<int>()->default_value("0"))
        ("threads", "threads", cxxopts::value<int>()->default_value("4"))
        ("overwrite", "overwrite the input matrix file", cxxopts::value<bool>()
          ->implicit_value("true")->default_value("false"))

        ("delta", "LLL: delta", cxxopts::value<double>()->default_value("0.999"))

        ("beta-start", "PBKZ: beta_start", cxxopts::value<int>()->default_value("15"))
        ("beta-shift", "PBKZ: beta_shift", cxxopts::value<int>()->default_value("6"))
        //("beta-shift-update", "PBKZ: flag for update beta_shift",
        //  cxxopts::value<bool>()->implicit_value("true")->default_value("false"))
        //("ignore-flat", "PBKZ: ignore-flat", cxxopts::value<bool>()
        //  ->implicit_value("true")->default_value("false"))
        //("logfile", "PBKZ: logfile", cxxopts::value<std::string>()->default_value(""))

        ("debug", "print more info.", cxxopts::value<bool>()
          ->implicit_value("true")->default_value("false"))

        ("h,help", "Print usage")
    ;
    options.parse_positional({"matrix", "block-size"});
    options.positional_help("<matrix> <block-size>");
    //options.show_positional_help();

    auto args = options.parse(argc, argv);

    if (args.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    bool debug = args["debug"].as<bool>();

    std::string input_file = args["matrix"].as<std::string>();
    std::string output_file = input_file;
    int block_size = args["block-size"].as<int>();
    int verbose = args["verbose"].as<int>();
    int vl = 0;
    int threads = args["threads"].as<int>();
    if (threads > omp_get_num_procs()) {
        threads = omp_get_num_procs();
    }

    bool overwrite = args["overwrite"].as<bool>();
    if (!overwrite) {
        output_file = output_file + ".reduced";
    }

    if (debug && verbose >= 1) {
        std::cout << "Input Param:" << std::endl;
        std::cout << "  matrix: " << input_file << std::endl;
        std::cout << "  block-size: " << block_size << std::endl;
        std::cout << "  threads: " << threads << std::endl;
        std::cout << std::endl;
    }


    //LatticeBasis<double> lB;
    LatticeBasis<long double> lB;


    if (FileExists(input_file)==false) {
        std::cerr << "coundn't find the file: " + input_file << std::endl;
        exit(-3);
    } else {
        LoadLattice(lB, input_file);
    }


    double lll__delta = args["delta"].as<double>();
    if (verbose >= 1) {
        std::cout << "Apply LLL with delta=" << lll__delta << std::endl;
    }
    vl = verbose >= 2 ? 1 : 0;
    if (debug) vl++;
    BigLLL(lB.L, 0, lll__delta, vl);
    if (verbose >= 1) {
        std::cout << std::endl;
    }
    SaveLattice(lB, output_file);


    int PBKZ__beta_start = args["beta-start"].as<int>();
    int PBKZ__beta_shift = args["beta-shift"].as<int>();
    std::string bkz_option = "parallel=" + to_stdstring(threads);
    bkz_option += " temporal=" + output_file;
    bkz_option += " startbeta=" + to_stdstring(PBKZ__beta_start);
    bkz_option += " betashift=" + to_stdstring(PBKZ__beta_shift);
    if (verbose >= 1) {
        std::cout << "Apply ProgressiveBKZ with block_size=" << block_size << std::endl;
    }
    vl = verbose >= 3 ? 1 : 0;
    if (debug) vl++;
    ProgressiveBKZ(lB, 0, block_size, vl, bkz_option);
    if (verbose >= 1) {
        std::cout << std::endl;
    }
    SaveLattice(lB, output_file);


    return 0;
}
