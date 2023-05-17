#include "dan_tests.hpp"


int main() {
    bool success = true;
    for (auto& group : dan_tests) {
        success &= group.run();
    }
    return success ? 0 : 1;
}

