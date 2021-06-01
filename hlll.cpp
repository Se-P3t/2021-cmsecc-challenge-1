#include <cxxopts.hpp>
#include <hplll/hlll.h>

using namespace hplll;



int main(int argc, char *argv[])
{
    cxxopts::Options options(argv[0], "Run hlll for a given matrix");

    options.add_options()
        ("matrix", "input matrix file", cxxopts::value<std::string>())
        ("out", "file to output lll-reduced matrix (overwrite the input file "
          "when not specified)", cxxopts::value<std::string>())
        ("threads", "threads", cxxopts::value<int>()->default_value("1"))
        //("seysen", "use Seysen reduction", cxxopts::value<bool>()
        //  ->implicit_value("true")->default_value("false"))

        ("precision", "bit precision to use", cxxopts::value<unsigned int>())

        ("delta", "LLL: delta", cxxopts::value<double>()
          ->default_value("0.999"))
        ("threshold", "block size threshold for parallelism", cxxopts::value<int>()
          ->default_value("1000000")) // no threads by default

        ("verbose", "print verbose info.", cxxopts::value<bool>()
          ->implicit_value("true")->default_value("false"))

        ("h,help", "Print usage")
    ;
    options.parse_positional({"matrix",});
    options.positional_help("<matrix>");
    //options.show_positional_help();

    auto args = options.parse(argc, argv);

    if (args.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    bool verbose = args["verbose"].as<bool>();
    std::string input_file = args["matrix"].as<std::string>();
    std::string output_file = input_file;
    if (args.count("out")) {
        output_file = args["out"].as<std::string>();;
    }
    int threads = args["threads"].as<int>();
    //bool seysen = args["seysen"].as<bool>();

    unsigned int precision;
    if (args.count("precision")) {
        precision = args["precision"].as<unsigned int>();
        mpfr_set_default_prec(precision);
    }

    double delta = args["delta"].as<double>();
    int threshold = args["threshold"].as<int>();


    filebuf fb;
    iostream os(&fb);
    ofstream outfile;

    typedef mpz_t integer_t;
    typedef matrix<Z_NR<integer_t> > MatrixZT;
    typedef matrix<FP_NR<mpfr_t> > MatrixFT;

    ZZ_mat<integer_t> A;
    ZZ_mat<integer_t> AT;

    fb.open(input_file, ios::in);
    os >> A;
    fb.close();
    if (verbose) {
        cout << "Loaded matrix file: "+input_file+" :: " <<  A.get_rows() << " x "
        << A.get_cols() << endl;
    }
    if (threads > 1 && threshold >= A.get_rows()) {
        cout << "Warning: the threshold is so big that there is no parallelism" << endl;
    }

    AT.resize(A.get_cols(), A.get_rows());
    transpose(AT, A);
    Lattice<integer_t, mpfr_t, MatrixZT, MatrixFT> B(AT, NO_TRANSFORM,
        DEF_REDUCTION);
    //    seysen ? SEYSEN_REDUCTION : DEF_REDUCTION);
    B.set_num_S(threads, threshold);

    B.hlll(delta, verbose);
    transpose(A, B.getbase());

    outfile.open(output_file, ios_base::trunc);
    outfile << A;
    outfile.close();

    return 0;
}



/*

g++ -O3 -Wall -pthread -fopenmp hlll.cpp -o hlll -I/usr/local/include/hplll -I/usr/local/include -L/usr/local/lib -lmpfr -lgmp -lquadmath /usr/local/lib/libfplll.so -lntl

*/
