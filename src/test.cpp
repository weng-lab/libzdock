#include "ZDOCK.hpp"
#include <exception>
#include <iostream>

int main() {
  try {
    // /Users/vanderva/git/zdockserver/webservice/test/irad/2MTA.zd.out
    // /Users/vanderva/Desktop/c0d92b1b-0888-4fa3-83a2-0ccc2f7e60af/zdock.out.pruned
    // /Users/vanderva/Desktop/0116d0ab47/job.154074.mzdock_24.out
    std::string file = "/Users/vanderva/Desktop/c0d92b1b-0888-4fa3-83a2-0ccc2f7e60af/zdock.out.pruned";
    zdock::ZDOCK z(file);
    auto &v = z.predictions();
    std::sort(v.begin(), v.end(),
              [](const auto &a, const auto &b) { return a.score > b.score; });
    std::cout << z << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    exit(1);
  }
}

