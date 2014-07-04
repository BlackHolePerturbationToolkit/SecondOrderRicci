#include <boost/multi_array.hpp>

void read_h1(std::vector<double> &r, std::vector<double> &f, std::vector<double> &fp,
             boost::multi_array<std::complex<double>,4> &h,
             boost::multi_array<std::complex<double>,4> &dh,
             boost::multi_array<std::complex<double>,4> &ddh);
