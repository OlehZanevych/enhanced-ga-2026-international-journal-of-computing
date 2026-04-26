#pragma once
#include <vector>
#include <cstdint>

// Frequency type for bi-weekly classes
enum class Freq : uint8_t { Weekly = 0, Numerator = 1, Denominator = 2 };

// A single class (gene) in the schedule
struct Gene {
    int req;   // requirement index
    int day;   // 0..DAYS-1
    int slot;  // 0..SLOTS_PER_DAY-1
    int room;  // 0..ROOMS-1
    Freq freq;
};

// A schedule is a vector of genes
using Schedule = std::vector<Gene>;

// Problem dimensions
struct Problem {
    int V;          // number of classes (genes)
    int L;          // number of lecturers
    int G;          // number of student groups
    int A;          // number of rooms
    int days;       // days per week
    int slots;      // time slots per day

    // req -> lecturer index
    std::vector<int> req_lecturer;
    // req -> list of student group indices
    std::vector<std::vector<int>> req_groups;
    // req -> frequency type (weekly / bi-weekly)
    std::vector<Freq> req_freq;
};
