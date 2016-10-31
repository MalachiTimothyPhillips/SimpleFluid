//
// Created by Malachi Phillips on 9/30/16.
//

#ifndef CFD_HW_RUNTIMEJOKE_H
#define CFD_HW_RUNTIMEJOKE_H

#include <string>

// Generate jokes for runtime

class RunTimeJoke{
public:
    unsigned int rng_; // random number associated with string

    std::string generate_runtime_joke();

private:
    //~virtual RunTimeJoke(){};
};

#endif //CFD_HW_RUNTIMEJOKE_H
